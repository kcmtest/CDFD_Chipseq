for i in ./*.bam; do java -jar -Xms12000m -Xmx18000m /src/picard.jar MarkDuplicates I=${i} O=${i}_rem.bam VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp METRICS_FILE=${i}_rem.txt REMOVE_DUPLICATES=true TMP_DIR=/tmp;done



#for i in ./*_rem.bam; do echo $i ;done

