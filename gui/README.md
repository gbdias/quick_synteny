# quick_synteny GUI: prototype

A desktop launcher for the pipeline, written in Electron. You fill in a form,
it runs `nextflow run main.nf` on this machine, shows each task as it runs,
and then opens the resulting `*.synteny.interactive.html` in its own window.

It needs nothing preinstalled: no Docker, conda, Java or Nextflow. The app
ships a pinned micromamba, and on first launch it uses it to install
Nextflow and Java, which takes under a minute. Each run then builds the
pipeline's tools as conda environments from `envs/*.yml` (the `conda`
profile). Docker stays available as an option.

It targets macOS and Linux; Windows isn't supported (see
[Platforms](#platforms)). The app isn't signed, and won't be for the
foreseeable future. So macOS asks users to allow it once per install and
per update, and the app can't update itself: it tells you when there's a
new version instead (see [Updates](#updates)).

## Platforms

| platform | build | status |
|---|---|---|
| macOS, Apple Silicon | `.dmg` | tested in full by hand; smoke-tested in CI |
| macOS, Intel | `.dmg` | smoke-tested here under Rosetta and in CI, on an Intel runner |
| Linux, x64 and arm64 | AppImage | smoke-tested in CI only; no full run yet |
| Windows | — | not supported |

On Linux, an AppImage may need `--no-sandbox` on recent Ubuntu versions,
which restrict the sandbox Electron uses (or an AppArmor profile). On
machines without FUSE, `APPIMAGE_EXTRACT_AND_RUN=1` makes it unpack itself
instead of mounting.

Windows would need more than a build. Nextflow needs a POSIX system (bash
for every task, Unix scripts in `bin/`), and miniprot has no Windows conda
build. So the app would have to drive everything inside WSL2, which is a
separate launch path. Windows users can try the Linux build, or the
command-line pipeline, inside WSL2, without support.

## Build

```bash
cd gui
npm install
npm test                  # unit tests (test/, node:test)
npm run dist:mac-arm64    # -> dist/quick_synteny-<version>-mac-arm64.dmg
npm run dist:mac-x64      # -> dist/quick_synteny-<version>-mac-x64.dmg
npm run dist:linux-x64    # -> dist/quick_synteny-<version>-linux-x86_64.AppImage
npm run dist:linux-arm64  # -> dist/quick_synteny-<version>-linux-arm64.AppImage
```

Each build first puts conda-forge's micromamba for that target, pinned and
sha256-checked (`scripts/fetch-micromamba.js`), into `vendor/`. A Mac can
build all four. Linux AppImages cross-build fine, but can only be run on
Linux. electron-builder then
packages the app. Everything it copies into the app's `Resources/` is set
in `package.json` (`build.extraResources`):

- the pipeline itself (`main.nf`, `nextflow.config`, `bin/`, `modules/`,
  `conf/`, `envs/`);
- `runtime.yml`;
- that micromamba.

Files are named `quick_synteny-<version>-<mac|linux>-<arch>.<ext>`
(`build.artifactName`). electron-builder writes an x64 AppImage's arch as
`x86_64`. Keep that pattern: the update notice looks for it, under either
spelling, in a release's files.

The app is ad-hoc signed, with the hardened runtime and
`build/entitlements.mac.plist`. Apple Silicon won't run an app with no
signature at all, so keep this even without a Developer ID. It doesn't
satisfy Gatekeeper, though.

To run from source instead, use `npm start`. The app then downloads
micromamba into the runtime folder on first launch.

### Smoke test

`quick_synteny --smoke-test` checks that a built app works on this machine,
without opening a window (`lib/smoke.js`). It:
1. checks that the bundled pipeline is complete and its scripts executable;
2. runs the real first-launch setup into `QS_HOME`, using the bundled
   micromamba to install Nextflow, Java and the plugins (weblog and
   nf-schema);
3. starts that Nextflow;
4. runs `nextflow run main.nf --help -profile conda` on the bundled
   pipeline.

It prints each step and exits 0 or 1. Point `QS_HOME` at a scratch folder,
or it sets up `~/.quick_synteny`:

```bash
QS_HOME=/tmp/qs dist/mac-arm64/quick_synteny.app/Contents/MacOS/quick_synteny --smoke-test
```

It doesn't run the pipeline itself: that needs NCBI and real genomes.

### CI

`.github/workflows/desktop-app.yml` runs on every push or pull request
that touches `gui/` or the pipeline (the app bundles it), and on demand:
1. the unit tests;
2. the four builds, each on a runner of its own platform and architecture
   (`macos-15`, `macos-15-intel`, `ubuntu-24.04`, `ubuntu-24.04-arm`);
3. a smoke test of each build. On Linux it runs under Xvfb, with
   `--no-sandbox` and `APPIMAGE_EXTRACT_AND_RUN`, since the runners have no
   display and no FUSE.

The builds are kept as run artifacts for 14 days. Nothing is published.

## Installing on macOS

The app isn't signed, so macOS blocks it the first time. This happens on
every install, **including each update**, since every downloaded copy is
new to macOS:

1. Open the `.dmg` and drag **quick_synteny** into Applications, replacing
   the old version if there is one.
2. Open it. macOS says it can't verify the app: close that dialog.
3. Open System Settings → Privacy & Security, and click **Open Anyway**
   next to quick_synteny. Confirm when asked.
4. From then on, that version opens like any other app.

If you're comfortable with a terminal, this does steps 2–3 in one go:
`xattr -dr com.apple.quarantine /Applications/quick_synteny.app`.

## Runtime

Everything the app installs lives under `~/.quick_synteny`. That is not
Electron's `userData` folder, which on macOS is under `Application
Support`: conda environments don't work reliably from a path with a space
in it.

| path | what | from |
|---|---|---|
| `bin/micromamba` | micromamba 2.9.0, only when running from source | conda-forge build, sha256 checked (`lib/runtime.js`) |
| `runtime/` | Nextflow 25.10.4, openjdk 21 | `runtime.yml` |
| `conda_envs/` | one env per `envs/*.yml` | built by Nextflow on first use |
| `nextflow/` | the app's Nextflow home (`NXF_HOME`): the `nf-weblog` and `nf-schema` plugins, run history | setup installs the plugins |
| `mamba/` | micromamba's package cache | shared by all of the above |

A packaged app uses its bundled micromamba instead of `bin/micromamba`.
`QS_HOME` moves the whole folder somewhere else, which is useful for
testing a first launch.

Nextflow runs with the runtime's `bin/` and micromamba first on PATH, and
with `MAMBA_ROOT_PREFIX` pointing into this folder. `MAMBA_NO_RC=true` is
set too. Without it, the user's own `~/.condarc` joins every solve: on this
Mac it added Anaconda's commercial `defaults` channel, which has licence
terms of its own and makes environments differ from one user to the next. The app adds its own
config to every run (`<outdir>/.gui.config`). That config enables the
weblog, and for the conda engine it sets `conda.useMicromamba` and
`conda.cacheDir`. **Settings** can point at your own Nextflow and Java
instead.

## Offline

The app needs a connection for three things:
- **First launch:** setup downloads Nextflow, Java and the plugins
  (`nf-weblog`, `nf-schema`).
- **The first run of each step:** it builds that step's tool environment.
- **Finding a reference on NCBI:** needed on every run that uses it.

Everything else works offline. With your own reference genome and
proteome, and the tools already downloaded, a whole run needs no
connection.

What happens without one:

| situation | what the app does |
|---|---|
| first launch | setup fails with "Setup needs an internet connection" and a **Retry setup** button, instead of micromamba's own message, which blames a corrupted package cache |
| name search, taxid check | "Can't reach NCBI to search by name" under the search box; the taxid check says it couldn't check, without blocking |
| **Run** with NCBI discovery | refused before starting: "Can't reach NCBI…", with the own-files alternative. Without that check, the run would fail only after about 3 minutes of `datasets` retrying |
| a step whose tools aren't downloaded yet | the error box says "Couldn't set up the tools for <step>", and that the download failed because the computer seems to be offline |
| any other failure while offline | the error box adds a hint that it looks like a network problem |
| update check | silent |

To find out whether the computer is online, the app asks NCBI and
conda-forge directly (`lib/network.js`). It does this before a discovery
run, when a run or setup fails, and when a run starts. A working network
connection says nothing about a firewall, a captive portal or a dead
Wi-Fi, so it doesn't rely on that. The tools' own messages don't always
say "network" either: with a cached index, micromamba only reports "Could
not load repodata.json … after retry".

When a run starts offline, Nextflow gets `NXF_OFFLINE=true`, so it doesn't
reach for anything itself. Offline Nextflow refuses an unpinned plugin, so
every plugin is pinned and installed by setup:
- the weblog (`nf-weblog@1.2.0` in `lib/runtime.js`), without which
  progress reporting would need the internet on every run;
- the pipeline's own plugins (`nf-schema`), read from the `plugins` block
  of the bundled `nextflow.config`, so their version is set in one place.

Each run's config repeats all of these pins. A `plugins` block in a `-c`
config replaces the pipeline's own instead of adding to it, and without
the repeat Nextflow loaded nf-schema unpinned, asking the registry for its
latest version. Nextflow's own
version check is always off (`NXF_DISABLE_CHECK_LATEST`): the app pins its
Nextflow, so a "26.x is available" notice would only mislead.

## Updates

On each launch, the app asks GitHub's releases API
(`api.github.com/repos/gbdias/quick_synteny/releases`) for the newest
release that is newer than itself and has a file for this platform and
architecture (`lib/updates.js`). Pipeline-only releases, like
`v0.1.0-alpha`, have no app files and are ignored.

If there is one, a banner offers:
- **Download**: that file, in the browser;
- **What's new**: the release page;
- **Dismiss**: don't mention this version again. A newer one still shows.

On macOS the banner also warns that the new version will need **Open
Anyway** again.

Checking never blocks or delays anything, and when GitHub can't be reached
the app says nothing. The page may only open this repo's GitHub links. A
switch under **Settings** turns the check off. `QS_RELEASES_URL` points it
at another feed, for testing.

The app's version is `package.json`'s `version`. A release tagged
`v<version>`, with the app files attached, is what the notice finds.

## What it does

- **First launch.** Setup starts by itself. A **Retry setup** button
  appears if it fails.
- **Environment check.** The header shows the runtime's Nextflow. Docker
  is only checked while the Docker engine is selected: once each time it's
  picked, and again when a Docker run starts, both on the page (for its
  status pill) and in the main process, which refuses the run if Docker
  isn't up. With Conda, Docker never comes up.
- **Form.** Target FASTA; either NCBI discovery (a species-name search that
  resolves to a taxid through NCBI Datasets' `taxon_suggest`, plus max rank
  and exclude-target) or your own reference and proteome; advanced chaining
  and miniprot options; **Tools** (Conda or Docker); output folder;
  `-resume`. The advanced options, the rank and the same-species switch
  each have a "?" like the result page's. Its text is the parameter's
  `description`, then its `help_text`, then its default, all from the
  bundled `nextflow_schema.json`: the same text as `--help <parameter>`.
  The tip itself is CSS only. Defaults and the rank list come from the
  schema too.
- **Input checks.** These run before a run starts (`lib/inputs.js`, and
  the `taxon:lookup` handler in `main.js`), so a bad input fails here and
  not minutes into a run:
  - The taxid is looked up at NCBI as it's typed. The field then shows the
    taxon's name and rank, "not a taxid" (e.g. `000000`), or "NCBI has no
    taxon"; the last two block **Run**. If NCBI can't be reached, the run
    isn't blocked, since the pipeline checks the taxid again.
  - Each picked file has to start with a `>` line, possibly gzipped.
  - Typed numbers are checked against the schema's types and limits, e.g.
    "Min block (anchors) must be at least 2". These are the rules the
    pipeline applies, so a bad value fails before the run starts rather
    than as Nextflow's error. If one does reach the pipeline, nf-schema's
    message goes in the error box.
  - An output folder that can't be read (macOS privacy settings) or
    written blocks the run. A folder that isn't empty and isn't a previous
    run gets a warning and a **Use a new subfolder instead** button. It
    suggests `<target>_vs_<taxid>`, so a run into `~/Downloads` doesn't
    scatter `work/` and `.nextflow/` there.
- **Memory preflight.** From the FASTA sizes, and README's ~10 GB of RAM per
  Gb of genome, it warns when miniprot won't fit and suggests a
  `--miniprot_chunk_gb` that will. The limit is this computer's RAM for
  Conda and Docker's VM memory for Docker.
- **Progress.** Nextflow's weblog POSTs one JSON event per task state
  change to a localhost port. The weblog doesn't report two things, so
  those come from the log instead:
  - tasks that `-resume` skips (`Cached process` lines);
  - conda environment builds (`Creating env using micromamba`), shown as a
    note above the task list.
- **Errors.** When a run fails, a box above the task list says why. Its
  text comes from Nextflow's own `ERROR ~` report in the console
  (`lib/nferror.js`):
  - the pipeline's message when it rejects its parameters;
  - for a failed step, the step's name and the last line of its
    `Command error:` output, with the full error when there's more than one
    line and a button to show its work folder.

  It also has a button to open `.nextflow.log`. When no report was printed,
  the box shows the log's last lines.
- **Cancel.** It sends SIGINT, like Ctrl-C: Nextflow stops its tasks (and
  containers) itself. Quitting the app cancels a run too.
- **Results.** It opens the page in a sandboxed window. The page's
  SVG/PNG/TSV exports open a save dialog that starts in Downloads. **Open
  result…** takes a page, or any earlier run's outdir.

## Tested

With the parameter schema (nf-schema), on macOS (Apple Silicon):
- The smoke test on an empty runtime folder installed both plugins and
  ran the schema-built `--help` offline.
- An existing runtime gained nf-schema by itself on the next launch.
- The form showed the schema's tips, defaults and rank list. Six
  out-of-range or malformed numbers were each blocked with their own
  message.
- A value forced past the form (`--min_block 1`) showed nf-schema's message
  in the error box.
- A full discovery run completed (23 tasks), and an offline run with your
  own files completed (15 tasks) with both plugins pinned and no registry
  requests.

Earlier, the packaged app was copied out of the `.dmg` and
launched with a Finder-like minimal PATH, an empty runtime folder, and
Docker not running:

- Setup started by itself, used the bundled micromamba, and was done in
  17 s.
- A full *S. cerevisiae* run took 127 s, downloading every package from
  scratch. All 23 tasks completed, and every links TSV, hit table, and
  `stats.json` matched the Docker run byte for byte.
- Input checks: `000000` and `12ab` are flagged as not taxids,
  `999999999` as unknown to NCBI, and `4932` as *S. cerevisiae*. Each bad
  input blocks **Run** with its own message.
- Error box: taxid `0`, forced past the form, is rejected by the pipeline
  with "--taxid must be an NCBI taxid (a positive whole number), got '0'".
  `999999999` fails in `PARSE_LINEAGE` and is reported as "--taxid
  999999999 is not a known NCBI taxid", with a button to its work folder.
- Update notice, against a mock releases feed:
  - it offers the newest release that has a `mac-arm64` file, skipping a
    newer Linux-only release and pipeline-only ones, with the macOS hint;
  - it stays quiet once that version is dismissed (also after a relaunch),
    when it's turned off in Settings, and when the feed can't be reached;
  - it refuses to open a link outside this repo.

  Against the real repo, which only has the pipeline-only `v0.1.0-alpha`,
  it shows nothing.
- Offline, simulated with `sandbox-exec` blocking outgoing connections and
  DNS but allowing localhost:
  - Setup on an empty runtime folder fails with "Setup needs an internet
    connection".
  - A run with your own reference and proteome and the tools already
    downloaded completes: 15 tasks in 21 s. Its links TSVs and hit tables
    match the online discovery run's byte for byte.
  - Missing tools give "Couldn't set up the tools for MINIPROT_INDEX…
    seems to be offline".
  - A discovery run is refused before it starts. Forced past that check,
    it fails after 187 s with "Couldn't reach NCBI" and the network hint.
  - Name search and taxid check say NCBI can't be reached, and the update
    check is silent.

  Two test-only adjustments: Electron ran with `--no-sandbox`, because its
  own sandbox can't start inside `sandbox-exec`, and Java was told to use
  IPv4, because `sandbox-exec` can't allow the `::ffff:127.0.0.1` loopback
  address Java uses for the weblog.
- The `.dmg` is 135 MB. Gatekeeper rejects the app once it's quarantined,
  as expected for an unsigned build.

Earlier, running from source with an empty `~/.quick_synteny`:

- **Set up runtime** took 27 s.
- A full *S. cerevisiae* run (taxid 4932, exclude-target) with Conda took
  66 s, including building all five tool environments. *S. pastorianus* was
  the discovered reference.
- Every output matches the Docker run of the same inputs byte for byte:
  miniprot's raw GFFs, both hit tables, all links TSVs, and `stats.json`.
  In the HTML page, only Bokeh's random element IDs differ.
- `~/.quick_synteny` takes 2.0 GB after that run, most of it micromamba's
  package cache.
- With Docker, as tested earlier: the page's ring PNG, ring SVG, and blocks
  TSV exports save correctly; cancelling mid-run leaves no containers
  behind; and a `-resume` rerun lists all 23 tasks as cached.

One run on one genome pair doesn't show that conda and the containers
always agree. `envs/*.yml` pins each tool's version, but each host solves
the dependencies and picks the builds itself. Per-platform lock files
would close that gap.

## Notes

- **The pipeline barely changed.** The only addition is a `conda` profile
  (`conf/conda.config`). The pipeline's one output is a self-contained,
  offline page that is interactive already, so the app is a thin launcher
  plus a viewer: no framework, no build step.
- **Without a container, macOS doesn't report task runtime metrics** (peak
  RSS and so on); Nextflow warns about this. Durations still show.
- **Icons** live in `build/`. macOS uses `icon-macos.png`, the padded
  rounded square that macOS expects (`mac.icon`). Linux uses the size set
  in `icons/`, and the full-bleed `icon.png` is electron-builder's default
  for anything else. The SVGs are the source artwork. When running from
  source, `main.js` sets the Dock icon itself.
- **`nextflow_schema.json` (nf-core style)** would let the form be generated
  rather than hand-written, and let other launchers (Seqera Platform,
  EPI2ME Desktop) run the pipeline as-is.

## Layout

- `main.js`: windows, IPC, and NCBI taxon lookup.
- `package.json`: scripts and the electron-builder config (`build`).
- `build/`: packaging resources: icons and macOS entitlements.
- `scripts/fetch-micromamba.js`: puts the pinned micromamba in `vendor/`
  for packaging.
- `lib/smoke.js`: `--smoke-test`.
- `test/`: unit tests, with fixtures of Nextflow's real console output.
- `runtime.yml`: the app's own Nextflow and Java.
- `lib/runtime.js`: micromamba download and verification, runtime setup,
  and the environment Nextflow runs with.
- `lib/env.js`: PATH fix-up, and deciding which Nextflow to launch;
  Docker detection.
- `lib/runner.js`: builds the `nextflow run` command line and the per-run
  config, runs the weblog listener, and handles cancel and result lookup.
- `lib/nferror.js`: turns Nextflow's console error report into the error
  box's title, message, and detail.
- `lib/inputs.js`: FASTA sniffing, output-folder inspection, and subfolder
  names.
- `test/schema.test.js`: checks that the schema, `nextflow.config`'s
  params, the form and the runner name the same parameters with the same
  defaults, and that the run config pins every plugin.
- `lib/updates.js`: the new-version check against GitHub Releases.
- `lib/network.js`: whether NCBI and conda-forge can be reached.
- `preload.js`: the IPC bridge exposed to the form as `window.qs`.
- `renderer/`: the form and progress UI (plain HTML/CSS/JS).
