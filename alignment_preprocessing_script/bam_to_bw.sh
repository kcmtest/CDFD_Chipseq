for i in ./*.bam; do bamCoverage -p 15 -b ${i} -o ${i}.bw ;done


#bamCoverage -b reads.bam -o coverage.bw
