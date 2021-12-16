#ce qu'on attend:
##TSS enrichment: zones nucleosome free, donc trouver zones sans nucleosomes donc séquencées: on en attend beaucoup
##insert size distribution: devrait être de la taille des éléments actifs? tout du moins des multiples de choses sachant que coupe à periodicité de nucléosome max. Devrait retrouver périodicité des tailles pour periodicité des coupures des nucléosomes. 

##Récupérer les nouveaux scripts
wget --user='tp_ngs' --password='Arabido2021!' https://flower.ens-lyon.fr/tp_ngs/arabidocontratac/Scripts/atac_qc.sh
wget --user='tp_ngs' --password='Arabido2021!' https://flower.ens-lyon.fr/tp_ngs/arabidocontratac/Scripts/plot_tlen.R
wget --user='tp_ngs' --password='Arabido2021!' https://flower.ens-lyon.fr/tp_ngs/arabidocontratac/Scripts/plot_tss_enrich.R

wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-52/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.52.gtf.gz
#dézipper le fichier
gunzip Arabidopsis_thaliana.TAIR10.52.gtf.gz

wget --user='tp_ngs' --password='Arabido2021!' https://flower.ens-lyon.fr/tp_ngs/arabidocontratac/Supporting_files/TAIR10_ChrLen.bed #format chromosome, start (TSS), stop (TES)
wget --user='tp_ngs' --password='Arabido2021!' https://flower.ens-lyon.fr/tp_ngs/arabidocontratac/Supporting_files/TAIR10_ChrLen.txt