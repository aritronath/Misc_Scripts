
library(data.table)
library(ggplot2)

# Read ssGSEA scores matrix (pathways in rows, samples in columns)
ssgsea.scores <- fread("ssGSEA_scores.tsv", header=T)

Path_list <- as.character(ssgsea.scores$V1) #Get list of pathways
ssgsea.scores <- as.matrix(ssgsea.scores[,-1]) #Remove pathway names from first column
rownames(ssgsea.scores) <- Path_list #Add pathways to rownames

# Load list of conditions/time points 
Conditions <- c(rep("A", 109), rep("B", 109)) #REPLACE WITH ACTUAL CONDITIONS/TIME POINTS

#### Calculate Kruskal-Wallis chi-square and P-values #####
Results <- data.frame(matrix(data=NA, nrow=length(Path_list), ncol=3))
colnames(Results) <- c("Pathway", "chi_sq", "pvalue")
Results$Pathway <- Path_list

t1 <- proc.time()
for (i in 1:length(Path_list)) {
	temp <- kruskal.test(ssgsea.scores[i,] ~ as.factor(Conditions))
	Results$chi_sq[i] <- temp$statistic
	Results$pvalue[i] <- temp$p.value
	print(i)
}
t2 <- proc.time() - t1

Results$FDR <- p.adjust(Results$pvalue, method='fdr')

write.csv(Results, file="Rank_test_results.csv")

#### Calulate ranks for plotting ####
ranked.scores <- apply(ssgsea.scores, 2, function (x) rank(x, ties.method='average'))
write.csv(ranked.scores, file="Ranked_ssGSEA_scores.csv")

#### Check rank trends with ssGSEA scores across samples ####
plot(ranked.scores[1,], ssgsea.scores[1,])

#### Plot pathway ranked score between conditions ####
# This example plots the first pathway
df <- data.frame("Ranks"=ranked.scores[1,], "Conditions"=Conditions) 
ggplot(df, aes(x=Conditions, y=Ranks, fill=Conditions)) +
	geom_violin() +
	geom_boxplot(width=0.1) + 
	theme_classic()



