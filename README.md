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
- BEDOPS `2.4.42`

> We recommend using https://seqera.io/containers/ to obtain the images.

## Execution

The input genome fasta files should be unzipped.

The tool is run using two files, `quick_synteny.config` and `quick_synteny.sh`.
- `quick_synteny.config`: is where the user specifies the paths to all necessary input files, including the apptainer images for the tools.
- `quick_synteny.sh`: running this script will perform the synteny analysis. First it will source the variables defined in the config file, so make sure this file is linked correctly.

Once all paths are specified correctly, execute the bash script like this:
```
bash quick_synteny.sh
# or
sbatch quick_synteny.sh
```

> If submitting the script to an HPC scheduler (e.g. SLURM or PBS) please add the appropriate headings.

________________________
### Example

We can create synteny plots between the ring-tailed lemur (_Lemur katta_, GCF_020740605.2) and human (GCF_009914755.1) using proteins from the Indian Elephan (_Elephas maximus indicus_, GCF_024166365.1). This is to show that even proteins from distant species can be used for this analysis in case the target or reference species are not annotated. 

This analysis took 23 minutes and a peak RAM of 20 GB, while using 48 CPUs. Two PDF files are generated, one dot plot synteny map, and one ideogram synteny map.


<img width="1116" height="1113" alt="Screenshot 2026-03-25 at 16 35 38" src="https://github.com/user-attachments/assets/fd27245c-0e32-40ea-85f3-d6a899030bc5" />

<img width="2494" height="756" alt="Screenshot 2026-03-25 at 16 36 14" src="https://github.com/user-attachments/assets/c9657a6b-8493-4205-94cd-ffa3f8c823af" />

