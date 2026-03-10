#!/bin/bash
#SBATCH --account=lipshutzlab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00
#SBATCH --mem=32GB
#SBATCH --output=logs/ragtag.log
#SBATCH --error=logs/ragtag.err
#SBATCH --job-name=ragtag

# Load required modules
module load Minimap2/2.24
module load samtools/1.21

# Activate conda environment (uncomment if needed)
# conda activate genomics

# Define variables
WORKDIR=/hpc/group/lipshutzlab/Leilton/jacana_RCD/scaffold_to_chromosome/chromosome_level
JACANA_SCAFFOLDS=/hpc/group/lipshutzlab/Leilton/jacana_sex_dispersal/SRA/ref_genome/Jacana.2cell.hap1.fa
REFDIR=$WORKDIR/reference_genomes
OUTDIR=$WORKDIR/chromosome_level
LOGDIR=$WORKDIR/logs

echo "[$(date)] Running RagTag scaffold..."
ragtag.py scaffold -t 10 -C \
    -o $WORKDIR/ragtag_output_chicken $REFDIR/Gallus_gallus.GRCg6a.dna.toplevel.fa $JACANA_SCAFFOLDS >> $LOGDIR/ragtag_scaffold_chicken.log 2>&1

# echo "[$(date)] Running RagTag merge with multiple references..."
# ragtag.py merge -t 10 -o $OUTDIR/final_assembly \
#     $WORKDIR/ragtag_output_zebra/ragtag.scaffold.agp \
#     $REFDIR/Gallus_gallus.GRCg6a.dna.toplevel.fa \
#     $REFDIR/Taeniopygia_guttata.bTaeGut1_v1.p.dna.toplevel.fa $JACANA_SCAFFOLDS >> $LOGDIR/ragtag_merge.log 2>&1

echo "[$(date)] Processing final assembly..."
cp $WORKDIR/ragtag_output_chicken/ragtag.scaffold.fasta $OUTDIR/final_assembly/jacana_spinosa.chromosome_level_chicken.fa
samtools faidx $OUTDIR/final_assembly/jacana_spinosa.chromosome_level_chicken.fa

# Generate assembly report
echo "=== Final Report ===" > $OUTDIR/final_assembly/assembly_report_chicken.txt
echo "Date: $(date)" >> $OUTDIR/final_assembly/assembly_report_chicken.txt
echo "Number of chromosomes: $(grep '^>' $OUTDIR/final_assembly/jacana_spinosa.chromosome_level_chicken.fa | wc -l)" >> $OUTDIR/final_assembly/assembly_report_chicken.txt
echo "Total size: $(awk '{sum+=$2} END {print sum}' $OUTDIR/final_assembly/jacana_spinosa.chromosome_level_chicken.fa.fai) bp" >> $OUTDIR/final_assembly/assembly_report_chicken.txt

# Standardize chromosome names
sed 's/>\(.*\)_RagTag/>chr_\1/' $OUTDIR/final_assembly/jacana_spinosa.chromosome_level_chicken.fa > $OUTDIR/final_assembly/jacana_spinosa.chromosome_level_final.fasta
