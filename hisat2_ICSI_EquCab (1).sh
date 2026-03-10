#HiSat2 for Horse Ensembl Assembly

gendir="/data1/EquCab/_ECA30/"
gen="/data1/EquCab/_ECA30/Equus_caballus.EquCab3.0.dna_sm.toplevel.fa"

indir="/data1/EquCab/FASTQ_RNA_public/"
outdir="/data1/EquCab/BAM_public_RNA_Ensembl/"
gendir="/data1/EquCab/_ECA30/"

rm hisat_sort_all.sh

declare -A samples
samples="1001-W 1002-W 1003-W 1004-W 1005-W 1006-W 1007-W 1008-W 1009-W 1010-W 1011-W 1012-W 1013-W 1014-W 1015-W 1016-W 1017-W 1018-W 1019-W 1020-W 1021-W 1022-W 1023-W 1026-W 133T-W 134T-W 135T-W 136T-W 137T-W 140T-W 143T-W 144T-W 145T-W 148T-W 150T-W 153T-W 155T-W 156T-W 160T-W 161T-W 162T-W 163T-W 168T-W 169T-W 170T-W 171T-W 172T-W 173T-W 175T-W 177T-W 178T-W 180T-W 181T-W 183V 187V-W 188V-W 189V 191V 192V 193V 194V 198V-W 199T-W 199V-W 1 2 300V-W 303T 304T 305T 306T 307T 308T 309T 310T 313T-W 314T-W 315T-W 316T-W 317T-W 319T-W 320T-W 321T-W 322V-W 323V 324T-W 325T-W 326T-W 327T-W 328T-W 333V-W 344V-W 3 4 A"
for sid in $samples; do
echo "(hisat2 -x ${gendir}Equus_caballus.EquCab3 -1 ${indir}${sra}_1.fq.gz -2 ${indir}${sra}_2.fq.gz -S ${dir}${sid}_${sra} --dta; samtools sort -T ${dir}${sid}_${sra}.sorting -o ${dir}${sid}_${sra}.sort.bam -O BAM ${dir}${sid}_${sra}; ) &" >> hisat_sort_all.sh
echo "sleep 600" >> hisat_sort_all.sh
done
unset samples
bash hisat_sort_all.sh



