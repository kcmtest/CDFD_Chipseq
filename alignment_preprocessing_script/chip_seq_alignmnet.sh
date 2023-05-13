#for i in $(ls *.fastq*.gz | cut -f1-4 -d"_" | uniq);

for i in $(ls *.fastq*.gz | cut -f1 -d"_" | uniq)
do 

 bowtie2 -q -p10 -x /run/media/punit/data1/BowtieIndex/hg38  -1 ${i}_R1_001.fastq.gz  -2 ${i}_R2_001.fastq.gz | samtools view -bS - >${i%}.bam  

# bowtie2 -q -p25 -x /run/media/punit/data1/BowtieIndex/hg38  -1 ${i}_R1.fastq.gz  -2 ${i}_R2.fastq.gz | samtools view -bS - >${i%}.bam  
done
