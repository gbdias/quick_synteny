#!/bin/bash

set -e

source quick_synteny.config

# check if outdir exists
mkdir -p ${OUTDIR}

# rename seqs on reference genome. ! this is finnicky !
# for NCBI genomes where the actual chromosome is at the end of the seqid line
# it's fine to comment it out if you already tidied up chr names beforehand
awk '
BEGIN { scaf = 0 }
/^>/ {
    # Check if "chromosome" exists in the line
    if ($0 ~ /[Cc]hromosome [0-9XYM]+/) {

        # Regex to capture the ID after "chromosome" and before the first comma or space
        match($0, /[Cc]hromosome ([0-9XYM]+)/, m)
        chr = m[1]

        # Check if it is unlocalized
        if ($0 ~ /unlocalized/) {
            u_counts[chr]++
            print ">chr" chr "_u" u_counts[chr]
        } else {
            print ">chr" chr
        }
    }
    # Fallback for scaffolds/organelles
    else {
        scaf++
        print ">scaf" scaf
    }
    next
}
{ print }
' ${REF} > ref_renamed.fa
REF=${OUTDIR}/ref_renamed.fa

# align proteome
singularity exec -B $BASE:$BASE ${MPROT} miniprot -t ${THREADS} -d target.mpi ${TARGET}
singularity exec -B $BASE:$BASE ${MPROT} miniprot -t ${THREADS} -I -u --gff target.mpi ${PROT} > ${TARGETNAME}.gff
singularity exec -B $BASE:$BASE ${MPROT} miniprot -t ${THREADS} -d ref.mpi ${REF}
singularity exec -B $BASE:$BASE ${MPROT} miniprot -t ${THREADS} -I -u --gff ref.mpi ${PROT} > ${REFNAME}.gff

# extract CDS sequences
singularity exec -B $BASE:$BASE ${AGAT} agat_sp_extract_sequences.pl -g ${TARGETNAME}.gff -f ${TARGET} -t cds -o ${TARGETNAME}_cds.fa
singularity exec -B $BASE:$BASE ${AGAT} agat_sp_extract_sequences.pl -g ${REFNAME}.gff -f ${REF} -t cds -o ${REFNAME}_cds.fa

# convert gff to BED
singularity exec -B $BASE:$BASE ${AGAT} agat_convert_sp_gff2bed.pl --gff ${TARGETNAME}.gff -o ${TARGETNAME}.bed
singularity exec -B $BASE:$BASE ${AGAT} agat_convert_sp_gff2bed.pl --gff ${REFNAME}.gff -o ${REFNAME}.bed
#singularity exec -B $BASE:$BASE ${BEDOPS} gff2bed < ${TARGETNAME}.gff > ${TARGETNAME}.bed
#singularity exec -B $BASE:$BASE ${BEDOPS} gff2bed < ${REFNAME}.gff > ${REFNAME}.bed

# select only sequences at least SIZE bp
# add mock features to reflect crom size in bed files
singularity exec -B $BASE:$BASE ${FASIZE} faSize -detailed ${TARGET} | awk -v s=${SIZE} '$2 >= s' | sort -k2,2nr > ${TARGETNAME}.chrom.sizes
cat ${TARGETNAME}.bed <(awk -v FS="\t" '{print $1"\t"($2-1)"\t"$2"\t"$1".chromEnd\t0\t."}' ${TARGETNAME}.chrom.sizes) > ${TARGETNAME}_chromEndMarked.bed
singularity exec -B $BASE:$BASE ${FASIZE} faSize -detailed ${REF} | awk -v s=${SIZE} '$2 >= s' | sort -k2,2nr > ${REFNAME}.chrom.sizes
cat ${REFNAME}.bed <(awk -v FS="\t" '{print $1"\t"($2-1)"\t"$2"\t"$1".chromEnd\t0\t."}' ${REFNAME}.chrom.sizes) > ${REFNAME}_chromEndMarked.bed

# get seqids
cut -f1 ${TARGETNAME}.chrom.sizes | awk '{print}' ORS=',' | sed 's/,$/\n/' > seqids
cut -f1 ${REFNAME}.chrom.sizes | awk '{print}' ORS=',' | sed 's/,$//' >> seqids

# reformat fasta
singularity exec -B $BASE:$BASE ${JCVI} python -m jcvi.formats.fasta format ${TARGETNAME}_cds.fa ${TARGETNAME}.cds
singularity exec -B $BASE:$BASE ${JCVI} python -m jcvi.formats.fasta format ${REFNAME}_cds.fa ${REFNAME}.cds

# pairwise search
singularity exec -B $BASE:$BASE ${JCVI} python -m jcvi.compara.catalog ortholog ${TARGETNAME} ${REFNAME} --no_strip_names

# create layout file
cat << EOF > layout
# y, xstart, xend, rotation, color, label, va,  bed
 .6,     .1,    .95,       0,      , ${TARGETNAME}, top, ${TARGETNAME}_chromEndMarked.bed
 .4,     .1,    .95,       0,      , ${REFNAME}, bottom, ${REFNAME}_chromEndMarked.bed
# edges
e, 0, 1, ${TARGETNAME}.${REFNAME}.anchors.simple
EOF

# create anchors
singularity exec -B $BASE:$BASE ${JCVI} python -m jcvi.compara.synteny screen --minspan=30 --simple ${TARGETNAME}.${REFNAME}.anchors ${TARGETNAME}.${REFNAME}.anchors.simple

# plot
singularity exec -B $BASE:$BASE ${JCVI} python -m jcvi.graphics.karyotype --basepair --keep-chrlabels --diverge Spectral --style white --chrstyle roundrect --font Arial --figsize 20x10 -o ${TARGETNAME}.${REFNAME}.karyotype.pdf seqids layout
