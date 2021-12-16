#!/bin/bash
#
# Compute basic quality metrics for ATAC-seq data
# 


# >>> change the value of the following variables


#SCRIPT VALEUR 1 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# General variables à définir avant de lancer script
workingDir=/home/rstudio/mydatalocal/TPNGSepigenomic1CRESSOT
scriptDir=/home/rstudio/mydatalocal/TPNGSepigenomic1CRESSOT/src
outputDir=${workingDir}/processed/ATACseq_qualite

ID=2020_372_S2_R # sample ID
#ID=2020_380_S10_R et remplacer ID à chaquoi fois
#ID=SRR4000473
bam_suffix=_alignment_shortmarked_duplicatesfiltered_duplicates.bam


gtf=${workingDir}/data/Arabidopsis_thaliana.TAIR10.52.gtf #gtf fichier avec gènes annotés:chromosome, gene/transcrit,transposon, coordonnées, gene_id, role du gene
selected_regions=${workingDir}/data/TAIR10_selectedRegions.bed
genome=${workingDir}/data/TAIR10_ChrLen.txt #info sur la longueur des chromosomes

#verification de la definition des variables
echo $selected_regions
echo $genome
echo $gtf
head $selected_regions
head $genome
head $gtf

# Variables for TSS enrichment
width=1000 #taille du fragemnt qu'on regarde autour des TSS dans le fichier GTF
flanks=100 #bords des fragemnts qu'on regarde pour voir position du read rapport au TSS pour voir enrichissement en TSS, si distribution centré sur TSS

# Variables for insert size distribution
chrArabido=${workingDir}/data/TAIR10_ChrLen.bed
grep -v -E "Mt|Pt" ${chrArabido} > ${workingDir}/data/TAIR10_ChrLen_1-5.bed #enleve Mt et Pt
chrArabido=${workingDir}/data/TAIR10_ChrLen_1-5.bed

#verification des nouvelles variables
echo $chrArabido
head $chrArabido

#////////////////////// Start of the script

mkdir -p ${outputDir} #crée nouveau fichier de sortie des fichiers

bam=${workingDir}/results/results_trimming/duplicat_marqués/${ID}${bam_suffix}
samtools view ${bam} | head #permet de verifier variable mais compressé donc la on decompresse avant 




# ------------------------------------------------------------------------------------------------------------ #
# --------------------------- Compute TSS enrichment score based on TSS annotation --------------------------- #
# ------------------------------------------------------------------------------------------------------------ #

#1. Define genomic regions of interest
#dit que c'est un fichier tabulé, colonne qui nous interesse à chercher, numero 9, split decoupe la colonne, stocké dans vecteur A, coupe selon ;, dans premier morceau découpé extraire ce qui est AT qlq chose avec split, match renvoit coordonnées dans chaine de caractère, extrait motif entre les slash, et on extrait ça, tout pour extraire le nom du gène 
#extrait nom du gène, regarde colonne qui correspond à orientatio, si en forward affiche tant et manip, sinon fait autre truc
#bedtools: fonction slop: + ou - permet de switcher pour comparer le comparable selon orientation du read (?)
#intersect: donne les intersections regions avec le genome et ce qu' #on a selectionné 

echo "-------------------------- Define genomic regions of interest"
grep -v "^[MtPt]" ${gtf} | awk '{ if ($3=="gene") print $0 }'  |\
grep "protein_coding" |\
awk ' BEGIN { FS=OFS="\t" } { split($9,a,";") ; match(a[1], /AT[0-9]G[0-9]+/) ; id=substr(a[1],RSTART,RLENGTH) ; if ($7=="+") print $1,$4,$4,id,$7 ; else print $1,$5,$5,id,$7 } ' |\
uniq | bedtools slop -i stdin -g ${genome} -b ${width} > ${outputDir}/tss_${width}.bed
bedtools intersect -u -a ${outputDir}/tss_${width}.bed -b ${selected_regions} > ${outputDir}/tmp.tss && mv ${outputDir}/tmp.tss ${outputDir}/tss_${width}.bed
echo `cat ${outputDir}/tss_${width}.bed | wc -l` "roi defined from" ${gtf}

tssFile=${outputDir}/tss_${width}.bed
head ${tssFile}

#attention! TSS defini en tant que bout 5' du gene, donc peut le prendre pour faire la fenetre quand le gene est en forward, mais si en sens inverse c'est le TES que tu prends donc la tu fais la fenetre en 3' pour avoir le TSS


#2. Compute TSS enrichment
echo "-------- Compute per-base coverage around TSS"
#couverture pour chaque base dans l'intervalle --> position dans la sequence + read qui match
bedtools coverage -a ${tssFile} -b ${bam} -d > ${outputDir}/${ID}_tss_depth.txt

#tri données en fonction de postion dans la sequence
awk -v w=${width} ' BEGIN { FS=OFS="\t" } { if ($5=="-") $6=(2*w)-$6+1 ; print $0 } ' ${outputDir}/${ID}_tss_depth.txt > ${outputDir}/${ID}_tss_depth.reoriented.txt
#awk: manipuler les tableaux en bash

#tri les donnees en fonction de leur position dans la sequence
sort -n -k 6 ${outputDir}/${ID}_tss_depth.reoriented.txt > ${outputDir}/${ID}_tss_depth.sorted.txt

bedtools groupby -i ${outputDir}/${ID}_tss_depth.sorted.txt -g 6 -c 7 -o sum > ${outputDir}/${ID}_tss_depth_per_position.sorted.txt

norm_factor=`awk -v w=${width} -v f=${flanks} '{ if ($6<f || $6>(2*w-f)) sum+=$7 } END { print sum/(2*f) } ' ${outputDir}/${ID}_tss_depth.sorted.txt`
echo "Nf: " ${norm_factor}
awk -v w=${width} -v f=${flanks} '{ if ($1>f && $1<(2*w-f)) print $0 }' ${outputDir}/${ID}_tss_depth_per_position.sorted.txt | awk -v nf=${norm_factor} -v w=${width} 'BEGIN { OFS="\t" } { $1=$1-w ; $2=$2/nf ; print $0 }' > ${outputDir}/${ID}_tss_depth_per_position.normalized.txt
Rscript ${scriptDir}/plot_tss_enrich.R -f ${outputDir}/${ID}_tss_depth_per_position.normalized.txt -w ${width} -o ${outputDir}  




# ---------------------------------------------------------------------------------------- #
# ------------------------------- Insert size distribution ------------------------------- #
# ---------------------------------------------------------------------------------------- #

echo "-------- Compute insert size distribution"
samtools view -f 3 -F 16 -L ${chrArabido} ${bam} | awk ' function abs(v){ return v < 0 ? -v : v } { print abs($9) } ' | sort -g | uniq -c | sort -k2 -g > ${outputDir}/${ID}_TLEN_1-5.txt
Rscript ${scriptDir}/plot_tlen.R -f ${outputDir}/${ID}_TLEN_1-5.txt -o ${outputDir}



# End of the script \\\\\\\\\\\\\\\\\\\\\\\\\\\\






# SCRIPT VALEUR 2 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
workingDir=/home/rstudio/mydatalocal/TPNGSepigenomic1CRESSOT
scriptDir=/home/rstudio/mydatalocal/TPNGSepigenomic1CRESSOT/src
outputDir=${workingDir}/processed/ATACseq_qualite

ID=2020_380_S10_R # sample ID
#ID=2020_380_S10_R et remplacer ID à chaquoi fois
#ID=SRR4000473
bam_suffix=_alignment_shortmarked_duplicatesfiltered_duplicates.bam


gtf=${workingDir}/data/Arabidopsis_thaliana.TAIR10.52.gtf #gtf fichier avec gènes annotés:chromosome, gene/transcrit,transposon, coordonnées, gene_id, role du gene
selected_regions=${workingDir}/data/TAIR10_selectedRegions.bed
genome=${workingDir}/data/TAIR10_ChrLen.txt #info sur la longueur des chromosomes

#verification de la definition des variables
echo $selected_regions
echo $genome
echo $gtf
head $selected_regions
head $genome
head $gtf

# Variables for TSS enrichment
width=1000 #taille du fragemnt qu'on regarde autour des TSS dans le fichier GTF
flanks=100 #bords des fragemnts qu'on regarde pour voir position du read rapport au TSS pour voir enrichissement en TSS, si distribution centré sur TSS

# Variables for insert size distribution
chrArabido=${workingDir}/data/TAIR10_ChrLen.bed
grep -v -E "Mt|Pt" ${chrArabido} > ${workingDir}/data/TAIR10_ChrLen_1-5.bed #enleve Mt et Pt
chrArabido=${workingDir}/data/TAIR10_ChrLen_1-5.bed

#verification des nouvelles variables
echo $chrArabido
head $chrArabido

#////////////////////// Start of the script

mkdir -p ${outputDir} #crée nouveau fichier de sortie des fichiers

bam=${workingDir}/results/results_trimming/duplicat_marqués/${ID}${bam_suffix}
samtools view ${bam} | head #permet de verifier variable mais compressé donc la on decompresse avant 



# ------------------------------------------------------------------------------------------------------------ #
# --------------------------- Compute TSS enrichment score based on TSS annotation --------------------------- #
# ------------------------------------------------------------------------------------------------------------ #

#1. Define genomic regions of interest
#dit que c'est un fichier tabulé, colonne qui nous interesse à chercher, numero 9, split decoupe la colonne, stocké dans vecteur A, coupe selon ;, dans premier morceau découpé extraire ce qui est AT qlq chose avec split, match renvoit coordonnées dans chaine de caractère, extrait motif entre les slash, et on extrait ça, tout pour extraire le nom du gène 
#extrait nom du gène, regarde colonne qui correspond à orientatio, si en forward affiche tant et manip, sinon fait autre truc
#bedtools: fonction slop: + ou - permet de switcher pour comparer le comparable selon orientation du read (?)
#intersect: donne les intersections regions avec le genome et ce qu' #on a selectionné 

echo "-------------------------- Define genomic regions of interest"
grep -v "^[MtPt]" ${gtf} | awk '{ if ($3=="gene") print $0 }'  |\
grep "protein_coding" |\
awk ' BEGIN { FS=OFS="\t" } { split($9,a,";") ; match(a[1], /AT[0-9]G[0-9]+/) ; id=substr(a[1],RSTART,RLENGTH) ; if ($7=="+") print $1,$4,$4,id,$7 ; else print $1,$5,$5,id,$7 } ' |\
uniq | bedtools slop -i stdin -g ${genome} -b ${width} > ${outputDir}/tss_${width}.bed
bedtools intersect -u -a ${outputDir}/tss_${width}.bed -b ${selected_regions} > ${outputDir}/tmp.tss && mv ${outputDir}/tmp.tss ${outputDir}/tss_${width}.bed
echo `cat ${outputDir}/tss_${width}.bed | wc -l` "roi defined from" ${gtf}

tssFile=${outputDir}/tss_${width}.bed
head ${tssFile}

#attention! TSS defini en tant que bout 5' du gene, donc peut le prendre pour faire la fenetre quand le gene est en forward, mais si en sens inverse c'est le TES que tu prends donc la tu fais la fenetre en 3' pour avoir le TSS


#2. Compute TSS enrichment
echo "-------- Compute per-base coverage around TSS"
#couverture pour chaque base dans l'intervalle --> position dans la sequence + read qui match
bedtools coverage -a ${tssFile} -b ${bam} -d > ${outputDir}/${ID}_tss_depth.txt

#tri données en fonction de postion dans la sequence
awk -v w=${width} ' BEGIN { FS=OFS="\t" } { if ($5=="-") $6=(2*w)-$6+1 ; print $0 } ' ${outputDir}/${ID}_tss_depth.txt > ${outputDir}/${ID}_tss_depth.reoriented.txt
#awk: manipuler les tableaux en bash

#tri les donnees en fonction de leur position dans la sequence
sort -n -k 6 ${outputDir}/${ID}_tss_depth.reoriented.txt > ${outputDir}/${ID}_tss_depth.sorted.txt

bedtools groupby -i ${outputDir}/${ID}_tss_depth.sorted.txt -g 6 -c 7 -o sum > ${outputDir}/${ID}_tss_depth_per_position.sorted.txt

norm_factor=`awk -v w=${width} -v f=${flanks} '{ if ($6<f || $6>(2*w-f)) sum+=$7 } END { print sum/(2*f) } ' ${outputDir}/${ID}_tss_depth.sorted.txt`
echo "Nf: " ${norm_factor}
awk -v w=${width} -v f=${flanks} '{ if ($1>f && $1<(2*w-f)) print $0 }' ${outputDir}/${ID}_tss_depth_per_position.sorted.txt | awk -v nf=${norm_factor} -v w=${width} 'BEGIN { OFS="\t" } { $1=$1-w ; $2=$2/nf ; print $0 }' > ${outputDir}/${ID}_tss_depth_per_position.normalized.txt
Rscript ${scriptDir}/plot_tss_enrich.R -f ${outputDir}/${ID}_tss_depth_per_position.normalized.txt -w ${width} -o ${outputDir}  

# ---------------------------------------------------------------------------------------- #
# ------------------------------- Insert size distribution ------------------------------- #
# ---------------------------------------------------------------------------------------- #

echo "-------- Compute insert size distribution"
samtools view -f 3 -F 16 -L ${chrArabido} ${bam} | awk ' function abs(v){ return v < 0 ? -v : v } { print abs($9) } ' | sort -g | uniq -c | sort -k2 -g > ${outputDir}/${ID}_TLEN_1-5.txt
Rscript ${scriptDir}/plot_tlen.R -f ${outputDir}/${ID}_TLEN_1-5.txt -o ${outputDir}






#SCRIPT VALEUR 3 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# General variables à définir avant de lancer script
workingDir=/home/rstudio/mydatalocal/TPNGSepigenomic1CRESSOT
scriptDir=/home/rstudio/mydatalocal/TPNGSepigenomic1CRESSOT/src
outputDir=${workingDir}/processed/ATACseq_qualite

ID=SRR4000473 # sample ID
bam_suffix=__alignment_shortmarked_duplicatesfiltered_duplicates.bam


gtf=${workingDir}/data/Arabidopsis_thaliana.TAIR10.52.gtf #gtf fichier avec gènes annotés:chromosome, gene/transcrit,transposon, coordonnées, gene_id, role du gene
selected_regions=${workingDir}/data/TAIR10_selectedRegions.bed
genome=${workingDir}/data/TAIR10_ChrLen.txt #info sur la longueur des chromosomes

#verification de la definition des variables
echo $selected_regions
echo $genome
echo $gtf
head $selected_regions
head $genome
head $gtf

# Variables for TSS enrichment
width=1000 #taille du fragemnt qu'on regarde autour des TSS dans le fichier GTF
flanks=100 #bords des fragemnts qu'on regarde pour voir position du read rapport au TSS pour voir enrichissement en TSS, si distribution centré sur TSS

# Variables for insert size distribution
chrArabido=${workingDir}/data/TAIR10_ChrLen.bed
grep -v -E "Mt|Pt" ${chrArabido} > ${workingDir}/data/TAIR10_ChrLen_1-5.bed #enleve Mt et Pt
chrArabido=${workingDir}/data/TAIR10_ChrLen_1-5.bed

#verification des nouvelles variables
echo $chrArabido
head $chrArabido

#////////////////////// Start of the script

mkdir -p ${outputDir} #crée nouveau fichier de sortie des fichiers

bam=${workingDir}/results/results_trimming/duplicat_marqués/${ID}${bam_suffix}
samtools view ${bam} | head #permet de verifier variable mais compressé donc la on decompresse avant 




# ------------------------------------------------------------------------------------------------------------ #
# --------------------------- Compute TSS enrichment score based on TSS annotation --------------------------- #
# ------------------------------------------------------------------------------------------------------------ #

#1. Define genomic regions of interest
#dit que c'est un fichier tabulé, colonne qui nous interesse à chercher, numero 9, split decoupe la colonne, stocké dans vecteur A, coupe selon ;, dans premier morceau découpé extraire ce qui est AT qlq chose avec split, match renvoit coordonnées dans chaine de caractère, extrait motif entre les slash, et on extrait ça, tout pour extraire le nom du gène 
#extrait nom du gène, regarde colonne qui correspond à orientatio, si en forward affiche tant et manip, sinon fait autre truc
#bedtools: fonction slop: + ou - permet de switcher pour comparer le comparable selon orientation du read (?)
#intersect: donne les intersections regions avec le genome et ce qu' #on a selectionné 

echo "-------------------------- Define genomic regions of interest"
grep -v "^[MtPt]" ${gtf} | awk '{ if ($3=="gene") print $0 }'  |\
grep "protein_coding" |\
awk ' BEGIN { FS=OFS="\t" } { split($9,a,";") ; match(a[1], /AT[0-9]G[0-9]+/) ; id=substr(a[1],RSTART,RLENGTH) ; if ($7=="+") print $1,$4,$4,id,$7 ; else print $1,$5,$5,id,$7 } ' |\
uniq | bedtools slop -i stdin -g ${genome} -b ${width} > ${outputDir}/tss_${width}.bed
bedtools intersect -u -a ${outputDir}/tss_${width}.bed -b ${selected_regions} > ${outputDir}/tmp.tss && mv ${outputDir}/tmp.tss ${outputDir}/tss_${width}.bed
echo `cat ${outputDir}/tss_${width}.bed | wc -l` "roi defined from" ${gtf}

tssFile=${outputDir}/tss_${width}.bed
head ${tssFile}

#attention! TSS defini en tant que bout 5' du gene, donc peut le prendre pour faire la fenetre quand le gene est en forward, mais si en sens inverse c'est le TES que tu prends donc la tu fais la fenetre en 3' pour avoir le TSS


#2. Compute TSS enrichment
echo "-------- Compute per-base coverage around TSS"
#couverture pour chaque base dans l'intervalle --> position dans la sequence + read qui match
#le fichier de l'article est trop gros donc on crée un intermediaire ou stocker pour pouvoir traiter les données
tssFile=${outputDir}/tss_${width}.bed
sort -k1,1 -k2,2n ${tssFile} > ${tssFile/.bed/.test}
tssFile=${outputDir}/tss_${width}.test
bedtools coverage -a ${tssFile} -b ${bam} -d -sorted > ${outputDir}/${ID}_tss_depth.txt

#tri données en fonction de postion dans la sequence
awk -v w=${width} ' BEGIN { FS=OFS="\t" } { if ($5=="-") $6=(2*w)-$6+1 ; print $0 } ' ${outputDir}/${ID}_tss_depth.txt > ${outputDir}/${ID}_tss_depth.reoriented.txt
#awk: manipuler les tableaux en bash

#tri les donnees en fonction de leur position dans la sequence
sort -n -k 6 ${outputDir}/${ID}_tss_depth.reoriented.txt > ${outputDir}/${ID}_tss_depth.sorted.txt

bedtools groupby -i ${outputDir}/${ID}_tss_depth.sorted.txt -g 6 -c 7 -o sum > ${outputDir}/${ID}_tss_depth_per_position.sorted.txt

norm_factor=`awk -v w=${width} -v f=${flanks} '{ if ($6<f || $6>(2*w-f)) sum+=$7 } END { print sum/(2*f) } ' ${outputDir}/${ID}_tss_depth.sorted.txt`
echo "Nf: " ${norm_factor}
awk -v w=${width} -v f=${flanks} '{ if ($1>f && $1<(2*w-f)) print $0 }' ${outputDir}/${ID}_tss_depth_per_position.sorted.txt | awk -v nf=${norm_factor} -v w=${width} 'BEGIN { OFS="\t" } { $1=$1-w ; $2=$2/nf ; print $0 }' > ${outputDir}/${ID}_tss_depth_per_position.normalized.txt
Rscript ${scriptDir}/plot_tss_enrich.R -f ${outputDir}/${ID}_tss_depth_per_position.normalized.txt -w ${width} -o ${outputDir}  




# ---------------------------------------------------------------------------------------- #
# ------------------------------- Insert size distribution ------------------------------- #
# ---------------------------------------------------------------------------------------- #

echo "-------- Compute insert size distribution"
samtools view -f 3 -F 16 -L ${chrArabido} ${bam} | awk ' function abs(v){ return v < 0 ? -v : v } { print abs($9) } ' | sort -g | uniq -c | sort -k2 -g > ${outputDir}/${ID}_TLEN_1-5.txt
Rscript ${scriptDir}/plot_tlen.R -f ${outputDir}/${ID}_TLEN_1-5.txt -o ${outputDir}



# End of the script \\\\\\\\\\\\\\\\\\\\\\\\\\\\
