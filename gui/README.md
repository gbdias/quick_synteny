# quick_synteny desktop app

A desktop launcher for the pipeline, written in Electron. You fill in a form,
it runs `nextflow run main.nf` on this machine, shows each step as it runs,
and opens the resulting `*.synteny.interactive.html` in its own window.

It needs nothing preinstalled: no Docker, conda, Java or Nextflow. The app
ships a pinned micromamba, and on first launch uses it to install Nextflow
and Java (under a minute). Each run builds the pipeline's tools as conda
environments from `envs/*.yml`; Docker stays available as an option.

What you get:
- a sidebar of runs, newest first, kept across restarts;
- a **New run** form: your genome, then a reference from your own files, a
  species you pick, or auto-discovery on NCBI. Inputs (taxids, accessions,
  FASTA files, numbers, the output folder, memory) are checked before the
  run starts, against the same rules as `nextflow_schema.json`;
- a run screen with a plain-language step checklist, what was found on NCBI,
  and, when a step fails, Nextflow's error in a readable box;
- the result page in a sandboxed window, with a thumbnail and summary
  numbers in the run's card.

## Platforms

| platform | build | status |
|---|---|---|
| macOS, Apple Silicon | `.dmg` | tested in full by hand; smoke-tested in CI |
| macOS, Intel | `.dmg` | smoke-tested under Rosetta and in CI |
| Linux, x64 and arm64 | AppImage | smoke-tested in CI only |
| Windows | — | not supported |

On recent Ubuntu, the AppImage may need `--no-sandbox`. Without FUSE,
`APPIMAGE_EXTRACT_AND_RUN=1` makes it unpack itself instead of mounting.

Windows isn't supported: Nextflow needs a POSIX system and miniprot has no
Windows build. The Linux build or the command-line pipeline may work inside
WSL2, without support.

## Installing on macOS

The app isn't signed, so macOS blocks it the first time, on every install
**including each update**:

1. Open the `.dmg` and drag **quick_synteny** into Applications.
2. Open it. macOS says it can't verify the app: close that dialog.
3. In System Settings → Privacy & Security, click **Open Anyway** next to
   quick_synteny, and confirm.

Or, in a terminal: `xattr -dr com.apple.quarantine /Applications/quick_synteny.app`.

## Build

```bash
cd gui
npm install
npm test                  # unit tests (test/, node:test)
npm start                 # run from source
npm run dist:mac-arm64    # -> dist/quick_synteny-<version>-mac-arm64.dmg
npm run dist:mac-x64
npm run dist:linux-x64    # -> dist/quick_synteny-<version>-linux-x86_64.AppImage
npm run dist:linux-arm64
```

Each build first puts a pinned, sha256-checked micromamba for that target
into `vendor/` (`scripts/fetch-micromamba.js`). electron-builder then bundles
the pipeline (`main.nf`, `nextflow.config`, `bin/`, `modules/`, `conf/`,
`envs/`), `runtime.yml` and that micromamba, as set in `package.json`'s
`build.extraResources`.

Keep the file name pattern `quick_synteny-<version>-<mac|linux>-<arch>.<ext>`:
the update check looks for it in a release's files. The app is ad-hoc signed
with the hardened runtime, which Apple Silicon requires even without a
Developer ID.

**Smoke test.** `quick_synteny --smoke-test` checks a built app without
opening a window (`lib/smoke.js`): the bundled pipeline is complete, setup
installs Nextflow, Java and the plugins, and `nextflow run main.nf --help`
works. It exits 0 or 1:

```bash
QS_HOME=/tmp/qs dist/mac-arm64/quick_synteny.app/Contents/MacOS/quick_synteny --smoke-test
```

**CI.** `.github/workflows/desktop-app.yml` runs the unit tests, the four
builds (each on a runner of its own platform) and a smoke test of each, on
every change to `gui/` or the pipeline. Builds are kept as run artifacts for
14 days; nothing is published.

**Icons.** macOS uses `build/icon.icns`, built from `build/icon-macos.svg` by
`npm run icons` (macOS only; commit the result). Linux uses `build/icons/`.

## Runtime

Everything the app installs lives under `~/.quick_synteny` (not Electron's
`userData`, since conda environments break in paths with spaces):

| path | what |
|---|---|
| `bin/micromamba` | micromamba, only when running from source (`lib/runtime.js`) |
| `runtime/` | Nextflow and Java, from `runtime.yml` |
| `conda_envs/` | one env per `envs/*.yml`, built on first use |
| `nextflow/` | the app's `NXF_HOME`: the `nf-weblog` and `nf-schema` plugins |
| `mamba/` | micromamba's package cache |

`QS_HOME` moves this folder, and `QS_USER_DATA` the app's settings and run
history, e.g. to test a first launch. `MAMBA_NO_RC=true` keeps the user's own
`~/.condarc` out of every solve. Each run gets an extra config,
`<outdir>/.gui.config`, that enables the weblog (how progress reaches the
app) and sets the conda options. **Settings** can point at your own Nextflow
and Java instead.

## Offline

A connection is needed for the first launch (setup), the first run of each
step (its tool environment), and anything from NCBI. With your own reference
and proteome and the tools already downloaded, a whole run works offline.
Every plugin is pinned and installed by setup so that offline Nextflow
(`NXF_OFFLINE=true`) accepts them; the run config repeats those pins. When
something fails offline, the app says so instead of showing the tools' own
messages (`lib/network.js`).

## Updates

On launch the app asks GitHub's releases API for the newest release with a
file for this platform (`lib/updates.js`) and, if there is one, shows a
banner with **Download**, **What's new** and **Dismiss**. It can't update
itself, since it isn't signed. Releases without app files are ignored; a
release must be tagged `v<version>` from `package.json`. The check is silent
when GitHub can't be reached, can be turned off in **Settings**, and
`QS_RELEASES_URL` points it at another feed for testing.
