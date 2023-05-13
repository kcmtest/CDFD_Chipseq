for x in `ls *.bam` ; do echo "print current:$x"; bedtools bamtobed -i "$x" > "${x%.bam}.bed"; done
