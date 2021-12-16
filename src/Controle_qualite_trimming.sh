#Controle qualité du trimming des reads des séquences de TP 372 et 380, et de l'article 473
cd processed
mkdir Trimming_qualite
cd Trimming
fastqc *trimmed.fastq*
mv *fastqc* ../Trimming_qualite
cd ../Trimming_qualite/
multiqc . #analyse qualité sur tous les fichiers en meme temps
cd ../
mkdir multiqc_trimmed
cd Trimming_qualite/
mv multiqc* ../multiqc_trimmed #dossier multiqc_trimmed avec analyse qualité des read trimmed

#On ne prendra pas en compte les séquences unpaired car aberrantes, pas interessantes pour nous
