##Complete list of genes with p-values and fold change

df1 <- read.csv("~/R-tables/gennimm_baseline_exp.csv")
df2 <- df1[,c(2,72,73)] #Need search function for gene
yscale <- min(df2[df2[,2]>0,2]) #Calculates maximum value for y-scaling
p=10

##Calculate genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off

df2$threshold <- as.factor(abs(df2[,3]) > 2 & df2[,2] < 0.05/ncol(df1))

#Subsets upregulated genes meeting bonferroni cutoff
Cutoff_genes = df2[abs(df2[,3]) > 2 & df2[,2] < 0.05/ncol(df1),]
#Sorts dataframe from highest highest and lowest fold change in expression
Cutoff_genes2 <- Cutoff_genes[order(-Cutoff_genes[,3]),]

#Chooses top ten upregulated and downregulated genes
topten <- Cutoff_genes2[1:p,]
bottomten <- Cutoff_genes2[(nrow(Cutoff_genes2)-p):nrow(Cutoff_genes2),]
dd_text <- rbind(topten, bottomten)

###############################################################################################################################
## Makes Volcano Plot

library(ggplot2)
library(RcolorBrewer)
#Displays colors: display.brewer.all()

g <- ggplot(df2, aes(x=df2[,3], y=-log10(df2[,2])))+ 
  geom_point(alpha=0.7, size=3, aes(colour=threshold)) +
  scale_color_manual(values=c("FALSE"="#666666","TRUE"="#339999")) +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill=NA),
    panel.background = element_rect(fill = "white"),
    axis.text=element_text(size=12),   
    axis.title.x = element_text(size = 16, face="bold"),
    axis.title.y = element_text(size = 16, face="bold", vjust = 1.5))+
  xlim(c(-35, 35)) + ylim(c(0, -log10(yscale)*1.1)) +
  xlab("Log2 Fold Change") + ylab("-Log10 p-value") + 
  geom_vline(xintercept=-2, color="grey", size=1) + geom_vline(xintercept=2, color="grey", size=0)

#Adds labels
g + geom_text(data=dd_text, aes(x=dd_text[,3], y=-log10(dd_text[,2]), label=dd_text$SYMBOL), colour="black")

###################################################################################################################################
## Constructs Histogram for fold expression change (y) versus gene ontology

#Combines top and bottom 100 differentially expressed genes meeting Bonferroni cutoff (Cutoffgenes2)
top100 <- Cutoff_genes2[1:100,]
bottom100 <- Cutoff_genes2[(nrow(Cutoff_genes2)-99):nrow(Cutoff_genes2),]
Histogram <- rbind(top100, bottom100)

#Defines y-axis by corresponding row numbers attaches for 
h_abcissa <- c(1:nrow(Histogram))
Histogram <- cbind(Histogram, h_abcissa)

#Searches rows for gene name
Histogram_36 <- Histogram[which(Histogram[,1]=="CD36"),]

#Creates color flag for CD36
Histogram$Color <- rownames(Histogram) %in% rownames(Histogram_36[1,])

#Constructs histogram
library(RcolorBrewer)
library(ggplot2)

h <- ggplot(Histogram, aes(x=h_abcissa)) + 
  geom_bar(aes(y=Histogram[,3]),  stat="identity", color="#339999") +
  geom_bar(aes(y=Histogram[,3], fill=Histogram$Color),  stat="identity") +
  scale_fill_manual(values=c("FALSE"="#339999","TRUE"="red")) +
  ylim(-35,15) +
  xlab("Gene ID") +
  ylab("Log2 Fold Change") +
  geom_text(aes(label = "CD36", y = Histogram[which(Histogram[,6]=="TRUE"),3]*1.2, x = Histogram[which(Histogram[,6]=="TRUE"),5]), color = "black") +
  theme(
    panel.border = element_rect(colour = "black", fill=NA),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, face="bold"),
    axis.title.y = element_text(size = 16, face="bold", vjust = 1.5),
    axis.ticks = element_blank(),
    legend.position="none",
    legend.title = element_blank(),
    legend.key = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_blank())

###################################################################################################################################

