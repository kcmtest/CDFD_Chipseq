library(idr)

url1 <- "https://www.encodeproject.org/files/ENCFF003ZKO/@@download/ENCFF003ZKO.bed.gz"
file1 <- "peaks_rep1.bed.gz"
if (!file.exists(file1)) download.file(url1, file1)
url2 <- "https://www.encodeproject.org/files/ENCFF093VDW/@@download/ENCFF093VDW.bed.gz"
file2 <- "peaks_rep2.bed.gz"
if (!file.exists(file2)) download.file(url2, file2)

f1 <- Ctrl_H2bUb_R2_sorted_peaks

library(readr)
library(GenomicRanges)
df1 <- read_delim("Ctrl_input/Ctrl_input/Ctrl_H2bUb_R1_sorted_peaks.narrowPeak", delim="\t", col_names=FALSE)
df2 <- read_delim("Ctrl_input/Ctrl_input/Ctrl_H2bUb_R2_sorted_peaks.narrowPeak", delim="\t", col_names=FALSE)
peak1 <- GRanges(df1$X1, IRanges(df1$X2, df1$X3), score=df1$X7)
peak2 <- GRanges(df2$X1, IRanges(df2$X2, df2$X3), score=df2$X7)
peak1 <- keepStandardChromosomes(peak1, pruning.mode="coarse")
peak2 <- keepStandardChromosomes(peak2, pruning.mode="coarse")

fo <- findOverlaps(peak1, peak2)
length(peak1)
length(peak2)

length(fo)

table(duplicated(from(fo)))

fo <- as.data.frame(fo)
fo <- fo[!duplicated(fo$queryHits) & !duplicated(fo$subjectHits),]

y1 <- peak1$score[fo[,1]]
y2 <- peak2$score[fo[,2]]
plot(y1, y2, cex=.1)

plot(log10(y1), log10(y2), cex=.1)

plot(rank(-y1), rank(-y2), cex=.1)


library(idr)
dat <- cbind(log10(y1), log10(y2))
dat <- dat[sample(nrow(dat),5000),]
system.time({ 
  res <- est.IDR(dat, mu=3, sigma=1, rho=.9, p=.5)
})


df <- data.frame(rep1=dat[,1],rep2=dat[,2],
                 rank1=rank(-dat[,1]),rank2=rank(-dat[,2]),
                 idr=res$idr)

library(ggplot2)
ggplot(df, aes(rep1,rep2,col=idr)) + geom_point()


ggplot(df, aes(rank1,rank2,col=idr)) + geom_point()


