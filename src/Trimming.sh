#TRIMMING
##enlève restes d'adapteurs, extremités à pauvre qualité

#chercher la fonction trimmomatic
ls /softwares/Trimmomatic-0.39/

#créer la variable Trimmomatic pour pas avoir à l'appeler à chaque fois
Trimmomatic=/softwares/Trimmomatic-0.39/trimmomatic-0.39.jar 
java -jar $Trimmomatic

##trimming sur toutes les données (article+données TP)
##avec sliding window pour parcourir le read et couper ce qui est en dessous d'une qualité (MINLEN 36)
#!/bin/bash
# arg1: number of threads
# to run: 
# chmod +x trim.sh
# <path>/trim.sh <number of threads>
# Example: ./trim.sh 40
cd processed
mkdir Trimming
cd ../data
output_dir=/home/rstudio/mydatalocal/TPNGSepigenomic1CRESSOT/processed/Trimming #rediriger les nouveaux fichiers dans le nouveau dossier trimming
Nextera=/softwares/Trimmomatic-0.39/adapters/NexteraPE-PE.fa #declare variable Nextera pour aller le chercher plus facilement

for f in data/*1.fastq #boucle pour trimmer tous les fastq
do
n=${f%%1.fastq} #prend le début du nom de fichier pour renomer les nouveaux fichiers
prefixe=${n/"data/"/""} #enelever data pour ne pas rester dans ce fichier sans trouver les autres
java -jar ${Trimmomatic} PE -threads 6 -phred33 ${n}1.fastq  ${n}2.fastq \
${output_dir}/${prefixe}1_trimmed.fastq $output_dir/${prefixe}1_unpaired.fastq ${output_dir}/${prefixe}2_trimmed.fastq \
${output_dir}/${prefixe}2_unpaired.fastq ILLUMINACLIP:${Nextera}:2:30:10 \
SLIDINGWINDOW:4:15 MINLEN:25
#trimming

done