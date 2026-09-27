#!/usr/bin/env python3
"""Regenerate conf/containers.config from envs/*.yml via Seqera Containers.

Each envs/<name>.yml is a conda environment (channels + pinned
dependencies) and the one place a tool's version is written. For each, this
asks Seqera Containers (the public Wave API behind seqera.io/containers, no
token needed) for its image in every combination the pipeline runs:

    Docker image      linux/amd64, linux/arm64
    Singularity SIF   linux/amd64, linux/arm64

Wave returns the same frozen image for the same package list, so rerunning
this is safe and only builds what is new (a new version, a new
environment). A SIF is written as a direct HTTPS download of its registry
blob rather than its oras:// reference, which the cluster's Singularity
refused ("could not get image manifest", 2026-09-23).

The generated config gives processes labelled env_<name> their image,
chosen when each task runs by container engine and by the host's
architecture (the Nextflow JVM's -- on a local Docker host also the
engine's).

To update a tool: edit its envs/<name>.yml, run this, commit both files.

    tools/update_containers.py            # all environments
    tools/update_containers.py miniprot   # only envs/miniprot.yml

Standard library only.
"""
import argparse
import json
import pathlib
import sys
import time
import urllib.error
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parent.parent
ENVS = ROOT / 'envs'
OUT = ROOT / 'conf' / 'containers.config'

WAVE = 'https://wave.seqera.io'
REGISTRY = 'community.wave.seqera.io'
BLOB_URL = 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/{short}/{digest}/data'
PLATFORMS = ('linux/amd64', 'linux/arm64')
FORMATS = ('docker', 'sif')
BUILD_TIMEOUT_S = 3600


def read_env(path):
    """(channels, dependencies) of a conda environment file -- the flat
    `key:` / `  - item` subset of YAML envs/*.yml is written in."""
    lists, key = {}, None
    for raw in path.read_text().splitlines():
        line = raw.split('#', 1)[0].rstrip()
        if not line.strip():
            continue
        if not line.startswith(' ') and line.endswith(':'):
            key = line[:-1].strip()
            lists[key] = []
        elif line.strip().startswith('- ') and key:
            lists[key].append(line.strip()[2:].strip())
        else:
            sys.exit(f"ERROR: {path}: unsupported line {raw!r} (only channels/dependencies lists)")
    if not lists.get('dependencies'):
        sys.exit(f"ERROR: {path}: no dependencies")
    return lists.get('channels', []), lists['dependencies']


def http_json(url, body=None, headers=None):
    req = urllib.request.Request(url, data=None if body is None else json.dumps(body).encode(),
                                 headers={'Content-Type': 'application/json', **(headers or {})})
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return json.load(r)
    except urllib.error.HTTPError as e:
        sys.exit(f"ERROR: {url}: HTTP {e.code}: {e.read().decode(errors='replace')[:500]}")


def request_image(channels, dependencies, platform, fmt):
    """Frozen Seqera Containers image for this environment, built first if
    Wave doesn't have it yet."""
    d = http_json(f"{WAVE}/v1alpha2/container", {
        'packages': {'type': 'CONDA', 'entries': dependencies, 'channels': channels},
        'freeze': True,
        'containerPlatform': platform,
        'format': fmt,
    })
    image = d.get('targetImage') or d.get('containerImage')
    if not image:
        sys.exit(f"ERROR: Wave returned no image for {dependencies} {platform} {fmt}: {d}")
    if not d.get('cached') and d.get('buildId'):
        wait_for_build(d['buildId'], image)
    return image


def wait_for_build(build_id, image):
    print(f"  building {image} ({build_id}) ...", file=sys.stderr)
    start = time.time()
    while time.time() - start < BUILD_TIMEOUT_S:
        s = http_json(f"{WAVE}/v1alpha1/builds/{build_id}/status")
        if s.get('status') == 'COMPLETED':
            if not s.get('succeeded'):
                sys.exit(f"ERROR: build {build_id} failed -- see {WAVE}/view/builds/{build_id}")
            return
        time.sleep(15)
    sys.exit(f"ERROR: build {build_id} still running after {BUILD_TIMEOUT_S} s")


def sif_blob_url(oras_ref):
    """oras://community.wave.seqera.io/<repo>:<tag> -> the SIF's HTTPS blob URL."""
    repo, tag = oras_ref.removeprefix(f"oras://{REGISTRY}/").rsplit(':', 1)
    token = http_json(f"https://cerbero.seqera.io/auth/token?service={REGISTRY}"
                      f"&scope=repository:{repo}:pull")['token']
    manifest = http_json(f"https://{REGISTRY}/v2/{repo}/manifests/{tag}",
                         headers={'Authorization': f"Bearer {token}",
                                  'Accept': 'application/vnd.oci.image.manifest.v1+json'})
    sifs = [l for l in manifest.get('layers', []) if l.get('mediaType', '').endswith('.sif')]
    if len(sifs) != 1:
        sys.exit(f"ERROR: {oras_ref}: expected one SIF layer, got {manifest.get('layers')}")
    digest = sifs[0]['digest'].removeprefix('sha256:')
    return BLOB_URL.format(short=digest[:2], digest=digest)


def read_existing():
    """label -> (block text) from the current config, kept for the
    environments not being regenerated."""
    blocks = {}
    if OUT.exists():
        current = None
        for line in OUT.read_text().splitlines(keepends=True):
            if line.startswith("    withLabel: 'env_"):
                current = line.split("'")[1]
                blocks[current] = ''
            if current:
                blocks[current] += line
                if line.rstrip() == '    }':
                    current = None
    return blocks


def config_block(name, deps, images):
    arm = "System.getProperty('os.arch') == 'aarch64'"
    return (f"    withLabel: 'env_{name}' {{\n"
            f"        // {' '.join(deps)}\n"
            f"        container = {{ workflow.containerEngine in ['singularity', 'apptainer']\n"
            f"            ? ({arm}\n"
            f"                ? '{images['linux/arm64', 'sif']}'\n"
            f"                : '{images['linux/amd64', 'sif']}')\n"
            f"            : ({arm}\n"
            f"                ? '{images['linux/arm64', 'docker']}'\n"
            f"                : '{images['linux/amd64', 'docker']}') }}\n"
            f"    }}\n")


HEADER = """\
// GENERATED by tools/update_containers.py from envs/*.yml -- do not edit by
// hand: change the environment file, rerun the script, commit both.
//
// Processes labelled env_<name> run envs/<name>.yml's Seqera Containers
// image for the container engine (Docker, or Singularity/Apptainer as an
// HTTPS SIF download) and the host's architecture (linux/arm64 on an
// aarch64 JVM such as Apple Silicon, linux/amd64 otherwise).

"""


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('names', nargs='*', help='environments to regenerate (default: all of envs/)')
    args = parser.parse_args()

    all_names = sorted(p.stem for p in ENVS.glob('*.yml'))
    names = args.names or all_names
    for n in names:
        if n not in all_names:
            sys.exit(f"ERROR: no envs/{n}.yml")

    blocks = read_existing()
    for name in names:
        channels, deps = read_env(ENVS / f"{name}.yml")
        print(f"{name}: {' '.join(deps)}", file=sys.stderr)
        images = {}
        for platform in PLATFORMS:
            for fmt in FORMATS:
                image = request_image(channels, deps, platform, fmt)
                if fmt == 'sif':
                    image = sif_blob_url(image)
                images[platform, fmt] = image
                print(f"  {platform} {fmt}: {image}", file=sys.stderr)
        blocks[f"env_{name}"] = config_block(name, deps, images)

    stale = sorted(set(blocks) - {f"env_{n}" for n in all_names})
    for label in stale:
        del blocks[label]
    OUT.write_text(HEADER + "process {\n" + "\n".join(blocks[k] for k in sorted(blocks)) + "}\n")
    print(f"wrote {OUT.relative_to(ROOT)}", file=sys.stderr)


if __name__ == '__main__':
    main()
