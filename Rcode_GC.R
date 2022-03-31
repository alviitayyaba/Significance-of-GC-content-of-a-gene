install.packages("seqinr", repos="http://R-Forge.R-project.org")
library(Biostrings)
library(seqinr)

##download following sequences from NCBI
f <- read.fasta("mart_export.txt") ##file containing random genes sequences 
ABO <- read.fasta("ABO.fa.txt")    ## file containing sequence of ABO gene
PITX2 <- read.fasta("PITX2.txt")   ## file containing sequence of PITX2 gene
lenABO <- ABO[[1]]
lenABO <- length(lenABO)
lenABO
lenPITX2 <- PITX2[[1]]
lenPITX2 <- length(lenPITX2)
lenPITX2
##avglen <- (lenABO + lenPITX2)/2
##avglen
summary(f[[1]])  # first sequence
summary(f[[1]])$length  # length of first sequence
f1 <- f[which(getLength(f) < 26000 )]       #change range for a different size of regions. 
f2 <- f1[which(getLength(f1) > 24000 )]
length(f2)   # number of sequences that fall within the range mentioned above
genename <- names(f1)  #name of the genes(region)
result <- vector("list", length(f2))
for (i in 1:length(f2)) {
  result[[i]] <- (summary(f2[[i]])$GC)
}
result[[1]]
data <- as.data.frame(result)
head(data)
##data contains GC content of only genes which have similar length as of genes of interest
Gc_content <- as.matrix(data)
head(Gc_content)
Gc_content <- t(Gc_content)
meanGC <- mean(Gc_content)
#the mean value of all genes GC content
meanGC
length <- vector("list", length(f2))
for (i in 1:length(f2)) {
  length[[i]] <- (summary(f2[[i]])$length)
}
length <- as.data.frame(length)
length <- as.matrix(length)
head(length)
meanlength <- mean(length)
#Mean length of all genes
meanlength
GCABO <- GC(ABO[[1]])
GCABO
GCPITX2 <- GC(PITX2[[1]])
GCPITX2
t_test <- t.test(Gc_content, mu = GCABO, alternative = "two.sided")
t_test
t_test2 <- t.test(Gc_content, mu = GCPITX2, alternative ="two.sided")
t_test2
png(file= "histogramGC.png")
hist(Gc_content, col = "Yellow", border ="blue", )
abline(v= GCABO, col="red")
abline(v= GCPITX2, col="green")
abline(v=meanGC, col="purple")
legend(0.50,27,c("ABO=0.487", "PITX2=0.493", "MeanGC=0.4558"), lwd=c(5,2,3), col=c("red", "green","purple")) #change values accordingly
dev.off()
png(file= "histogramLength.png")
hist(length, col = "Yellow", border ="blue", )
abline(v= lenABO, col="red")
abline(v= lenPITX2, col="green")
abline(v=meanlength, col="purple")
legend(25100,12,c("ABO= 24801", "PITX2 = 24701", "MeanLength=24988"), lwd=c(5,2,3), col=c("red", "green","purple"))
dev.off()

