#peakcalling
#regarde regions tres enrichies, pics utilisés pour faire trucs + fins, detection + fine par la suite
#regarde les fichiers indexés pour l'analyse de peakcalling
#INDEXER LES FICHIERS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
samtools index 2020_372_S2_R_alignment_shortmarked_duplicatesfiltered_duplicates.bam
samtools index 2020_380_S10_R_alignment_shortmarked_duplicatesfiltered_duplicates.bam
samtools index SRR4000473__alignment_shortmarked_duplicatesfiltered_duplicates.bam


#script PEAK CALLING \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
PeakcallingDir=/home/rstudio/mydatalocal/TPNGSepigenomic1CRESSOT/results/Peakcalling
mkdir $PeakcallingDir

bam_suffix=_alignment_shortmarked_duplicatesfiltered_duplicates.bam

workingDir=/home/rstudio/mydatalocal/TPNGSepigenomic1CRESSOT/results/results_trimming/duplicat_marqués
for f in ${workingDir}/*filtered_duplicates.bam
do
IDbase="$(basename -- $f)"
#echo ${IDbase/$bam_suffix/}
macs2 callpeak -f "BAM" -t ${f} -n ${IDbase/$bam_suffix/} --outdir ${PeakcallingDir} -q 0.01 --nomodel --shift -25 --extsize 50 --keep-dup "all" -B --broad --broad-cutoff 0.01 -g 10E7
done



#associer les pics au gènes avec bedtools closest
#ASSOCIER PICS AUX GENES \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#gtf donne coordonnées, nom du gene, orientation pour comparaison, filtre le gtf avant, garde que les gènes dans le fichier 
#on extrait le nom du gène à la fin
gtf=/home/rstudio/mydatalocal/TPNGSepigenomic1CRESSOT/data/Arabidopsis_thaliana.TAIR10.52.gtf
gtf_filtered=${gtf/.gtf/.filtered.gtf}
grep -v "^[MtPt]" ${gtf} | awk '{ if ($3=="gene") print $0 }' |\
awk ' BEGIN { FS=OFS="\t" } { split($9,a,";") ; match(a[1], /AT[0-9]G[0-9]+/) ; id=substr(a[1],RSTART,RLENGTH) ; print $1,$4,$5,id,$7 }' |\
sort -k1,1 -k2,2n > ${gtf_filtered}

for f in ${PeakcallingDir}/*broadPeak
do
#echo $f
bedtools closest -a $f -b ${gtf_filtered} -D ref > ${f/_peaks.broadPeak/}.nearest.genes.txt 
done
#D-ref pour savoir si overlap ou si plus loin




