#analyse résultats
##80% alignement sur le genome, considère bon

#Construction d'un index pour analise avec la fonction bowtie2 à partir du génome de référence pour l'alignement avec les reads par la suite
indexdir=/home/rstudio/mydatalocal/TPNGSepigenomic1CRESSOT/data/index
mkdir -p $indexdir
cd $indexdir
bowtie2-build -f ../Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz TAIR10

#déplacement dans fichiers trimmed pour alignement
cd ../../processed/Trimming

#ALIGNEMENT
# Aligning paired reads
alignment=../../results/results_trimming
mkdir -p $alignment

for file in *1_trimmed.fastq # for each sample
do 
  suffixe=${file%%1_trimmed.fastq} # enlever suffixe
  prefixe=${suffixe##"Trimming/"} #enlever prefix et remplacer
 bowtie2 -x ${indexdir}/TAIR10 -1 ${prefixe}1_trimmed.fastq -2 ${prefixe}2_trimmed.fastq -X 2000 \
  --very-sensitive -p 6 | samtools sort -@ 6 --output-fmt bam -o ${prefixe}_alignment_short.bam 
  #fonction bowtie, envoie dans fichier de sortie SAM
#pipe pour ne pas enregistrer ficheirs intermediaires SAM trop lourds, permet de sauver de la place en ne gardant que bam 

done


#473 avec beaucoup de 0 alignés comparé à 472: confirme la théorie de la contamination, on avait deja observé des porblèmes dans les controles qualité 

#réalisation de statistiques sur les données de l'index avec une option de samtools
for bamfile in results/results_trimming/*bam 
do
  samtools index ${bamfile}
  samtools idxstats ${bamfile}
done

#observe 50% presque des reads = reads des Mt et Pt: plus accessibles car pas sous forme chromatine, veut dire que beaucoup dans les échantillons malgré le fait que exprimentalement selectionne les noyaux, donc protocole peut etre amélioré 
#garder les bons mapping: samtools view
#permet de garder ou degager les regions
#téléchargement du fichier qui nous permet de repérer les regions blacklistées (séquences répétées) où les reads ne vont pas être très informatifs et vont embrouiller
wget --user='tp_ngs' --password='Arabido2021!' https://flower.ens-lyon.fr/tp_ngs/Supporting_files/TAIR10_selectedRegions.bed

