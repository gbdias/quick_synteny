# quick_synteny
Whole-genome synteny plots based on alignment of protein-coding genes

## Overview
This script will use a **proteome** to produce synteny plots between a **target** and a **reference** genome.
- `proteome`: a set of proteins you expect to match both the target and reference genomes reasonably well. These could be proteins derived from either target or reference, or even from a third related species.
- `target`: a genome assembly you want to analyze. Could be contigs, scaffolds, or chromosome-level. Bear in mind that extremelly fragmented assemblies might be hard to visualize.
- `reference`: usually a chromosome-level genome assembly.

Briefly, the script will align the proteome to both genomes using `miniprot`, extract the CDS from each hit using `AGAT`, find orthologous loci using `jcvi.compara.catalog ortholog` (LAST-based), and finally produce a karyotype visualization `jcvi.graphics.karyotype`.

## Requirements
This script must be run on a Linux operating system with Apptainer/Singularity installed.

Apptainer images for each of the individual tools must be obtained prior to execution of the script. The tool versions tested are indicated below.

- miniprot `0.18`
- AGAT `1.5.1`
- ucsc-fasize `482`
- jcvi `1.5.7`

> We reccommend using https://seqera.io/containers/ to obtain the images.

## Execution

The input genome fasta files should be unzipped.

The tool is run using two files, `quick_synteny.config` and `quick_synteny.sh`.
- `quick_synteny.config`: is where the user specifies the paths to all necessary input files, including the apptainer images for the tools.
- `quick_synteny.sh`: running this script will perform the synteny analysis. First it will source the variables defined in `quick_synteny.config`, so make sure this file is linked correctly.

Once all paths are specified correctly, execute the bash script like this:
```
bash quick_synteny.sh
# or
sbatch quick_synteny.sh
```

> If submitting the script to an HPC scheduler (e.g. SLURM or PBS) please add the appropriate headings.

________________________
### Example

Creating synteny plots between the ring-tailed lemur (_Lemur katta_) and human took 21 minutes and a peak RAM of 22.8Gb, while using 20CPUs.

