#controle qualité des reads depuis le dossier data
fastqc ../data/2020_372_S2_R1.fastq
fastqc ../data/2020_372_S2_R2.fastq
fastqc ../data/2020_380_S10_R1.fastq
fastqc ../data/2020_380_S10_R2.fastq

#mettre données dans le dossier processed
mv data/*fastqc* processed

#executer le script 
bash src/controle_qualite.sh

#analyser tous les fichiers en un fichier
cd processed
multiqc .

#ranger les fichiers
mkdir multiqc_qualite
mv multiqc* multiqc_qualite
mkdir fastqc_qualite
mv *fastqc* fastqc_qualite

#analyse résultats
##éviter d'avoir trop de duplicate reads: 006 pas top, 007 juste derrière
##échantillon de l'article (473): grosse difference %GC, contamination? comparer graphe avec % alignement sur le genome 
##473 et 472 replicats dans article 
##sequence quality: tout est bon
##taille des séquences: dépend de la protéine utilisée pour l'ATACseq
##duplicata levels: pas mal de gènes en plus de 100x quand meme (006 le pire toujours, avec 007)
##sequences surreprésentées dans 006 et 472, 473
##adapter content: adapter standards et communs à tous les reads, parfois séquence adapter, début ou fin des reads pas bien mappés sur génome, peut perdre extrémités de reads à cause des adapteurs. Retirer adapteurs avec logiciel avant le mapping. Ici vert donc pas besoin de les enlever. 
##enrichissement en reads de façon général si gène, à différents starting points pour high complecity // aux memes starting points si low complexity
