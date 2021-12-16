#comparer avec bedtools intersect chaque fichier WOX5 donc quiescentes souches avec 374 de la racine
annotationsDir=/home/rstudio/mydatalocal/TPNGSepigenomic1CRESSOT/results/annotations
mkdir -p $annotationsDir

#pour 006 quiescent/racine 374
PeakcallingDir=/home/rstudio/mydatalocal/TPNGSepigenomic1CRESSOT/results/Peakcalling

WOX5=2019_006_S6_R.nearest.genes.txt
racine=2020_374_S4.corrected.nearest.genes.txt

bedtools intersect -v -a ${PeakcallingDir}/$WOX5 -b ${PeakcallingDir}/$racine > ${annotationsDir}/${WOX5/.nearest.genes.txt/}_${racine/.nearest.genes.txt/}_difference.txt

#pour 007 quiescent/racine 374
WOX5=2019_007_S7_R.nearest.genes.txt 
racine=2020_374_S4.corrected.nearest.genes.txt

bedtools intersect -v -a ${PeakcallingDir}/$WOX5 -b ${PeakcallingDir}/$racine > ${annotationsDir}/${WOX5/.nearest.genes.txt/}_${racine/.nearest.genes.txt/}_difference.txt

#pour 372 quiescent/racine 374
WOX5=2020_372_S2_R.nearest.genes.txt  
racine=2020_374_S4.corrected.nearest.genes.txt

bedtools intersect -v -a ${PeakcallingDir}/$WOX5 -b ${PeakcallingDir}/$racine > ${annotationsDir}/${WOX5/.nearest.genes.txt/}_${racine/.nearest.genes.txt/}_difference.txt


#comparaison des resultats de cellules quiescentes pour avoir consensus 
bedtools intersect -v -a ${annotationsDir}/2019_006_S6_R_2020_374_S4.corrected_difference.txt -b ${annotationsDir}/2020_372_S2_R_2020_374_S4.corrected_difference.txt ${annotationsDir}/2019_007_S7_R_2020_374_S4.corrected_difference.txt > ${annotationsDir}/${2019_006_S6_R_2020_374_S4.corrected_difference.txt/.corrected_difference.txt/}_${2020_372_S2_R_2020_374_S4.corrected_difference.txt/.corrected_difference.txt/}_${2019_007_S7_R_2020_374_S4.corrected_difference.txt/.corrected_difference.txt/}_difference_quiescentes.txt






#ESSAI 2 AVEC BROADPEAK ET NON NEAREST GENE \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

#pour 006 quiescent/racine 374
PeakcallingDir=/home/rstudio/mydatalocal/TPNGSepigenomic1CRESSOT/results/Peakcalling

WOX5=2019_006_S6_R_peaks.broadPeak  
racine=2020_374_S4.corrected_peaks.broadPeak  

bedtools intersect -v -a ${PeakcallingDir}/$WOX5 -b ${PeakcallingDir}/$racine > ${annotationsDir}/${WOX5/.nearest.genes.txt/}_${racine/.nearest.genes.txt/}_difference.txt

#pour 007 quiescent/racine 374
WOX5=2019_007_S7_R_peaks.broadPeak   
racine=2020_374_S4.corrected_peaks.broadPeak  

bedtools intersect -v -a ${PeakcallingDir}/$WOX5 -b ${PeakcallingDir}/$racine > ${annotationsDir}/${WOX5/.nearest.genes.txt/}_${racine/.nearest.genes.txt/}_difference.txt

#pour 372 quiescent/racine 374
WOX5=2020_372_S2_R_peaks.broadPeak 
racine=2020_374_S4.corrected_peaks.broadPeak  

bedtools intersect -v -a ${PeakcallingDir}/$WOX5 -b ${PeakcallingDir}/$racine > ${annotationsDir}/${WOX5/.nearest.genes.txt/}_${racine/.nearest.genes.txt/}_difference.txt



#comparaison des resultats de cellules quiescentes pour avoir consensus \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
bedtools intersect -v -a ${annotationsDir}/${sample1} -b ${annotationsDir}/${sample2} ${annotationsDir}/${sample3} > ${annotationsDir}/${2019_006_S6_R_2020_374_S4.corrected_difference.txt/.corrected_difference.txt/}_${2020_372_S2_R_2020_374_S4.corrected_difference.txt/.corrected_difference.txt/}_${2019_007_S7_R_2020_374_S4.corrected_difference.txt/.corrected_difference.txt/}_difference_quiescentes.txt


sample1=2019_006_S6_R_peaks.broadPeak_2020_374_S4.corrected_peaks.broadPeak_difference.txt
sample2=2020_372_S2_R_peaks.broadPeak_2020_374_S4.corrected_peaks.broadPeak_difference.txt
sample3=2019_007_S7_R_peaks.broadPeak_2020_374_S4.corrected_peaks.broadPeak_difference.txt

touch ${annotationsDir}/quiescentcells.common.peaks.txt
cat ${sample1} >> ${annotationsDir}/quiescentcells.common.peaks.txt >


