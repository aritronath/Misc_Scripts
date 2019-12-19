# Read files
rawScores <- read.delim("~/Desktop/Test_data/P15272_scores_RawCounts.txt", header = T, row.names = 1, check.names = F, fill = T)
normScores <- read.delim("~/Desktop/Test_data/P15272_scores_NormelizedCounts.txt", header = T, row.names = 1, check.names = F, fill = T)
annot <- read.delim("~/Desktop/Test_data/P15272_annot.txt", header=T, check.names = F, fill = T)

# Subset on common samples 
normScores <-  normScores[, intersect(colnames(normScores), colnames(rawScores))]
rawScores <-  rawScores[, intersect(colnames(normScores), colnames(rawScores))]

length(intersect(rownames(normScores), rownames(rawScores))) # 3703 of the pathways are in both raw and norm ssGSEA scores

# Subset on common pathways
norm2 <- normScores[intersect(rownames(normScores), rownames(rawScores)), ]
raw2 <- rawScores[intersect(rownames(normScores), rownames(rawScores)), ]

norm2 <- as.matrix(norm2)
raw2 <- as.matrix(raw2)

#-----------------------------------------
# Correlation between the scores 
TheCors <- array(data=NA, dim=nrow(norm2))
for (i in 1:nrow(norm2)) {
  TheCors[i] <- cor(as.numeric(norm2[i,]), as.numeric(raw2[i,]))
  print(i)
}

png("Correl_distribution.png", height=600, width=600, res=100)
hist(TheCors)
abline(v=mean(TheCors), lwd=3, col="red")
dev.off()

# Plot the scores 
png("All_Scores_Correl.png", width=1200, height=1200, res=200)
plot(c(raw2), c(norm2), ylab="Normalized ssGSEA scores", xlab="Raw ssGSEA scores", pch=19, 
     cex=0.1, col="#000FFF50")
dev.off()

cor(c(raw2), c(norm2)) #0.31

#-----------------------------------------
# Relative ranks of pathways
norm.rank <- apply(norm2, 2, rank)
raw.rank <- apply(raw2, 2, rank)

# Correlation between the ranks 
TheRankCors <- array(data=NA, dim=nrow(norm2))
for (i in 1:nrow(norm2)) {
  TheRankCors[i] <- cor(as.numeric(norm.rank[i,]), as.numeric(raw.rank[i,]))
  print(i)
}

png("RankCorrel_distribution.png", height=600, width=600, res=100)
hist(TheRankCors)
abline(v=mean(TheRankCors), lwd=3, col="blue")
dev.off()

# Plot the ranks 
png("All_Ranks_Correl.png", width=1200, height=1200, res=200)
plot(c(raw.rank), c(norm.rank), ylab="Normalized ssGSEA ranks", xlab="Raw ssGSEA ranks", pch=19, 
     cex=0.1, col="#000FFF50")
dev.off()

cor(c(raw.rank), c(norm.rank)) #0.24


#------------------------------------
#Differential scores between inferCNV clusters 
y <- match(colnames(raw2), annot$Cell.ID)
annot <- annot[y,]

TheAOV.raw <- array(data=NA, dim=nrow(raw2))
TheAOV.norm <- array(data=NA, dim=nrow(raw2))

for (i in 1:nrow(raw2)) {
  temp <- summary(aov(raw2[i,] ~ annot$Infer_HMM))
  TheAOV.raw[i] <- temp[[1]][1,5] 
  
  temp <- summary(aov(norm2[i,] ~ annot$Infer_HMM))
  TheAOV.norm[i] <- temp[[1]][1,5] 
  
  print(i)
}

# Plot the p-values 
png("All_AOV_P_Correl.png", width=1200, height=1200, res=200)
plot(-log10(TheAOV.raw), -log10(TheAOV.norm), ylab="Normalized ssGSEA P-values", xlab="Raw ssGSEA P-values", 
     cex=0.25, pch=19, col="#0000FF50")
dev.off()

cor(-log10(TheAOV.raw), -log10(TheAOV.norm)) #0.63

# Differential ranks between inferCNV clusters 
TheKW.raw <- array(data=NA, dim=nrow(raw2))
TheKW.norm <- array(data=NA, dim=nrow(raw2))

for (i in 1:nrow(raw2)) {
  
  temp <- kruskal.test(raw2[i,] ~ annot$Infer_HMM)
  TheKW.raw[i] <- temp$p.value
  
  temp <- kruskal.test(norm2[i,] ~ annot$Infer_HMM)
  TheKW.norm[i] <- temp$p.value
  
  print(i)
}

# Plot the p-values
png("All_KW_P_Correl.png", width=1200, height=1200, res=200)
plot(-log10(TheKW.raw), -log10(TheKW.norm), ylab="Normalized ssGSEA P-values", xlab="Raw ssGSEA P-values", 
     cex=0.25, pch=19, col="#0000FF50")
dev.off()

cor(-log10(TheKW.raw), -log10(TheKW.norm)) #0.5



###-----------------------------------------------
# Which pathways are consistent between the two datasets
library(GSEABase)
c2 <- read.delim("~/Desktop/Test_data/c2.all.v7.0.symbols.txt", header=F, fill=T, sep="\t", 
                 colClasses = "character", na.strings="")
hall <- read.delim("~/Desktop/GeneSets/h.all.v7.0.symbols.gmt", header=F, fill=T, sep="\t", 
                   colClasses = "character", na.strings="")

N.c2 <- apply(c2, 1, function (x) ncol(c2) - sum(is.na(x)))
N.hall <- apply(hall, 1, function (x) ncol(hall) - sum(is.na(x)))

names(N.c2) <- c2$V1
names(N.hall) <- hall$V1
N.gs <- c(N.c2, N.hall)

x <- match(rownames(raw2), names(N.gs))
plot(TheCors, log2(N.gs[x]))

N.gs[match(rownames(raw2)[which(TheCors > 0.9)], names(N.gs))]


rownames(raw2)[which(TheCors > 0.8)]

rownames(raw2)[which(TheAOV.norm < 1e-100 & TheAOV.norm < 1e-100)]

