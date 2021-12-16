#on veut enlever: genome non chromosomique
#read non mappés
#reads mauvaise qualité (473)
#régions blacklistées
#régions dupliquées

cd results
for mapfile in results_trimming/*bam
do 
#samtools view -b #-b permet output in bam
#marquage des duplicats de PCR
#mis un flag associé aux duplicats 
java -jar $PICARD MarkDuplicates \
      I=${mapfile} \
      O=${mapfile/".bam"/"marked_duplicates.bam"} \
      M=marked_dup_metrics.txt
done



#utilisation de https://broadinstitute.github.io/picard/explain-flags.html pour comprendre les SAM flags sur les dossiers bam ouverts
#

#enlever les régions: trier
#enelevr Mt et Pt des séquences de selected regions
cd data
grep -v -E "Mt|Pt" TAIR10_selectedRegions.bed > TAIR10_selectedRegions_nuc.bed #enlève sequences mitochondrie et chloroplastes des données, enlevé avec v, E pour regarder les 2 motifs

for mapfile in /home/rstudio/mydatalocal/TPNGSepigenomic1CRESSOT/results/results_trimming/duplicat_marqués/*
do 
echo $mapfile
samtools view -b ${mapfile} -o ${mapfile/".bam"/"filtered_duplicates.bam"} \
-L /home/rstudio/mydatalocal/TPNGSepigenomic1CRESSOT/data/TAIR10_selectedRegions_nuc.bed \
-F 1024 -f 3 -q 30
done
#L garde que les reeds qui mappent sur le nouveau fichier sans Mt et Pt
#F pour exclure duplicats de PCR avec flag 1024, f3 pour garder ce qui est apparié
#fichier selected regions: régions qu'on veut garder 
#fonction -L: si region correspond au fichier qu'on veut garder, il garde 

