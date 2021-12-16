# TP NGS epigenomic1
Analyse des données ATACseq obtenues pour le TP, ici sur des cellules de racine entière d'Arabidopsis thaliana. 
1) D'abord, téléchargement et controle qualité des données par analyse graphique des données obtenues en fastqc. 

Téléchargement des données: script de telechargement.txt dans src
Controle qualité des données téléchargées d'ATACseq du TP: script de controle_qualite.sh dans src
372 et 380: analyses de cellules des racines entières (378 aussi)
fichiers 2019_00X: analyses de cellules quiescentes de la racine

controle qualité: éviter d'avoir trop de duplicates

fichiers dont je m'occupe: 
2020_372_S2_R1
2020_372_S2_R2
2020_380_S10_R1
2020_380_S10_R2

2) Réalisation du trimming pour enlever par swipe window les parties de la séquence de pas assez bonne qualité, ou les extremités ou pourraient se trouver les primers de séquençage. Utilisation du script Trimming.sh dans src.

3) Controle qualité du travail de trimming sur les séquences de l'article et du TP par fastqc. Utilisation du script Controle_qualite_trimming.sh dans src.

4) Il faut ensuite réaliser le mapping: prendre tous les read obtenus dans ces fichiers pour alignement avec la séquence chromosomique obtenue dans le fichier Arabidopsis_thaliana.TAIR10 transformé dans l'index pour reconstruire la séquence, et observer où se trouvent les zones de la chromatine accessibles qui ont donc été séquencées.On utilise le script mapping.sh.

5) Filtrer les données du mapping, en retirant les reads de mauvaise qualité ou non informatifs, pour ne garder que ceux que l'on va analyser. Suivre le script filtering.sh dans src.

6) Effectuer le controle qualité des données ATACseq: utiliser script controle_qualite_ATACseq dans src pour télécharger les autres scripts nécéssaires. Regarder le script atac_qc.sh pour le code.

7) Peakcalling: pour trouver les pics de signal à travers le genome, soit les pics montrant les séquences les plus accessibles dans les differents genomes analysés. Ouvrir le script dans src Peakcalling.sh
puis ouverture du script Intersect.sh dans src pour comparer et ne retenir que les pics différents entre les lignées de cellules souches quiescentes et la lignée de cellules de toute la racine pour voir les différences dans les zones accessibles sur le genome des cellules souches avec le reste des cellules racinaires. On utilise ici la fonction bedtools intersect.

8) absente pour la dernière partie.

penser à préciser de quels fichiers besoin aux étapes et ou les trouver 