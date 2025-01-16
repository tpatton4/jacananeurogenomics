
## HiSat2 for Horse Ensembl Assembly
## Configured for GRACE
## FASTER uses module load SAMtools/1.16.1

rm hisat*.job

gendir="/scratch/user/agaricx/EquCab/_EquCab3/"
gen="/scratch/user/agaricx/EquCab/_EquCab3/Equus_caballus.EquCab3.0.dna_sm.toplevel.fa"
fqdir="/scratch/user/agaricx/EquCab/ICSI_RNA/FASTQ/"
out="/scratch/user/agaricx/EquCab/ICSI_RNA/BAM/"

samples="1001-W_S4 1002-W_S5 1003-W_S9 1004-W_S10 1005-W_S25 1006-W_S26 1007-W_S27 1008-W_S28 1009-W_S29 1010-W_S30 1011-W_S31 1012-W_S32 1013-W_S40 1014-W_S41 1015-W_S42 1016-W_S43 1017-W_S50 1018-W_S51 1019-W_S52 1020-W_S53 1021-W_S62 1022-W_S63 1023-W_S64 1026-W_S65 133T-W_S92 134T-W_S87 135T-W_S80 136T-W_S81 137T-W_S82 140T-W_S83 143T-W_S88 144T-W_S89 145T-W_S84 148T-W_S85 150T-W_S86 153T-W_S77 155T-W_S78 156T-W_S79 160T-W_S11 161T-W_S12 162T-W_S13 163T-W_S93 168T-W_S94 169T-W_S6 170T-W_S14 171T-W_S15 172T-W_S2 173T-W_S3 175T-W_S16 177T-W_S95 178T-W_S17 180T-W_S7 181T-W_S18 183V_S8 187V-W_S19 188V-W_S20 189V_S1 191V_S21 192V_S22 193V_S23 194V_S24 198V-W_S37 199T-W_S33 199V-W_S38 1_S72 2_S73 300V-W_S39 303T_S34 304T_S35 305T_S36 306T_S44 307T_S45 308T_S46 309T_S47 310T_S54 313T-W_S55 314T-W_S56 315T-W_S57 316T-W_S58 317T-W_S59 319T-W_S60 320T-W_S61 321T-W_S66 322V-W_S48 323V_S49 324T-W_S67 325T-W_S68 326T-W_S69 327T-W_S70 328T-W_S71 333V-W_S90 344V-W_S91 3_S74 4_S75 A_S76"

for sid in $samples; do

echo "#!/bin/bash 
##NECESSARY JOB SPECIFICATIONS 
#SBATCH --export=NONE
#SBATCH --job-name=rna_${sid}
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=248G 
#SBATCH --output=rna_${sid}.out.%j
#SBATCH --error=rna_${sid}.err.%j
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=bdavis@tamu.edu

module load GCC/11.3.0
module load OpenMPI/4.1.4
module load HISAT2/2.2.1
module load SAMtools/1.17
# module load SAMtools/1.16.1

hisat2 -x ${gendir}Equus_caballus.EquCab3 -1 ${fqdir}${sid}_R1.fq.gz -2 ${fqdir}${sid}_R2.fq.gz -S ${out}${sid}.bam --dta
samtools sort -T ${out}${sid}.sorting -o ${out}${sid}.sort.bam -O BAM ${out}${sid}.bam" > hisat_${sid}.job

done
