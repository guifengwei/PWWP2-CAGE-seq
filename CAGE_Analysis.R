
library(dplyr)
library(ggplot2)
library(gridExtra)

data_E14 <- read.table("E14_CAGE_vM22.xls", header=T)
data_A1  <- read.table("A1_CAGE_vM22.xls", header=T)
data_C7  <- read.table("C7_CAGE_vM22.xls", header=T)

n=dim(data_E14)[1]
print("test")
data <- rbind( data_E14 %>% select("seqnames", "start", "end", "width", "strand", "annotation", "genes", "tpm") %>% rename(chr=seqnames),
               data_A1  %>% select("seqnames", "start", "end", "width", "strand", "annotation", "genes", "tpm") %>% rename(chr=seqnames),
			   data_C7  %>% select("seqnames", "start", "end", "width", "strand", "annotation", "genes", "tpm") %>% rename(chr=seqnames) 
			)
data$Sample= factor( rep( c("E14", "A1", "C7"), each=n), levels=c("E14", "A1", "C7")  )

pdf("CAGE.plot.pdf", width=15, height=5)
data1 <- data[data$annotation=="promoter",]
q1 <- ggplot(data1, aes(x=annotation, y=tpm, fill=Sample)) + geom_boxplot(outlier.size=-1)+ylim(0,250)+scale_color_manual('Sample', values=c( "E14"="black", "A1"="#F8766D", "C7"="#619CFF"))
data2 <- data[data$annotation!="promoter",]
q2 <- ggplot(data2, aes(x=annotation, y=tpm, fill=Sample)) + geom_boxplot(outlier.size=-1)+ylim(0,20)+scale_color_manual('Sample', values=c( "E14"="black", "A1"="#F8766D", "C7"="#619CFF"))
grid.arrange(q1, q2, nrow=1, ncol=2)
dev.off()

############################
data_E14 <- data_E14 %>% select("seqnames", "start", "end", "width", "strand", "annotation", "genes", "tpm") %>% rename(chr=seqnames, E14_tpm=tpm)
data_A1  <- data_A1  %>% select("seqnames", "start", "end", "width", "strand", "annotation", "genes", "tpm") %>% rename(chr=seqnames, A1_tpm =tpm)
data_C7  <- data_C7  %>% select("seqnames", "start", "end", "width", "strand", "annotation", "genes", "tpm") %>% rename(chr=seqnames, C7_tpm =tpm)
Data <- left_join(data_E14, data_A1)
Data <- left_join(Data,    data_C7)

write.table(Data, file="CAGE_Table.xls", sep="\t", quote=F, col.names=T, row.names=F)

################# On, Up, Off, Down
print("###################################")
print("Comparing A1 and E14: Promoter")
print("On:")
dim(Data%>%filter(annotation=="promoter") %>% filter(A1_tpm > E14_tpm) %>% filter(E14_tpm==0))[1]
print("Up:")
dim(Data%>%filter(annotation=="promoter") %>% filter(A1_tpm > E14_tpm) %>% filter(E14_tpm>0))[1]
print("Off:")
dim(Data%>%filter(annotation=="promoter") %>% filter(A1_tpm < E14_tpm) %>% filter(A1_tpm==0))[1]
print("Down:")
dim(Data%>%filter(annotation=="promoter") %>% filter(A1_tpm < E14_tpm) %>% filter(A1_tpm>0))[1]
print("###################################")
print("Comparing C7 and E14: Promoter")
print("On:")
dim(Data%>%filter(annotation=="promoter") %>% filter(C7_tpm > E14_tpm) %>% filter(E14_tpm==0))[1]
print("Up:")
dim(Data%>%filter(annotation=="promoter") %>% filter(C7_tpm > E14_tpm) %>% filter(E14_tpm>0))[1]
print("Off:")
dim(Data%>%filter(annotation=="promoter") %>% filter(C7_tpm < E14_tpm) %>% filter(C7_tpm==0))[1]
print("Down:")
dim(Data%>%filter(annotation=="promoter") %>% filter(C7_tpm < E14_tpm) %>% filter(C7_tpm>0))[1]

print("###################################")
print("Comparing A1 and E14: exon")
print("On:")
dim(Data%>%filter(annotation=="exon") %>% filter(A1_tpm > E14_tpm) %>% filter(E14_tpm==0))[1]
print("Up:")
dim(Data%>%filter(annotation=="exon") %>% filter(A1_tpm > E14_tpm) %>% filter(E14_tpm>0))[1]
print("Off:")
dim(Data%>%filter(annotation=="exon") %>% filter(A1_tpm < E14_tpm) %>% filter(A1_tpm==0))[1]
print("Down:")
dim(Data%>%filter(annotation=="exon") %>% filter(A1_tpm < E14_tpm) %>% filter(A1_tpm>0))[1]
print("###################################")
print("Comparing C7 and E14: exon")
print("On:")
dim(Data%>%filter(annotation=="exon") %>% filter(C7_tpm > E14_tpm) %>% filter(E14_tpm==0))[1]
print("Up:")
dim(Data%>%filter(annotation=="exon") %>% filter(C7_tpm > E14_tpm) %>% filter(E14_tpm>0))[1]
print("Off:")
dim(Data%>%filter(annotation=="exon") %>% filter(C7_tpm < E14_tpm) %>% filter(C7_tpm==0))[1]
print("Down:")
dim(Data%>%filter(annotation=="exon") %>% filter(C7_tpm < E14_tpm) %>% filter(C7_tpm>0))[1]



print("###################################")
print("Comparing A1 and E14: intron")
print("On:")
dim(Data%>%filter(annotation=="intron") %>% filter(A1_tpm > E14_tpm) %>% filter(E14_tpm==0))[1]
print("Up:")
dim(Data%>%filter(annotation=="intron") %>% filter(A1_tpm > E14_tpm) %>% filter(E14_tpm>0))[1]
print("Off:")
dim(Data%>%filter(annotation=="intron") %>% filter(A1_tpm < E14_tpm) %>% filter(A1_tpm==0))[1]
print("Down:")
dim(Data%>%filter(annotation=="intron") %>% filter(A1_tpm < E14_tpm) %>% filter(A1_tpm>0))[1]
print("###################################")
print("Comparing C7 and E14: intron")
print("On:")
dim(Data%>%filter(annotation=="intron") %>% filter(C7_tpm > E14_tpm) %>% filter(E14_tpm==0))[1]
print("Up:")
dim(Data%>%filter(annotation=="intron") %>% filter(C7_tpm > E14_tpm) %>% filter(E14_tpm>0))[1]
print("Off:")
dim(Data%>%filter(annotation=="intron") %>% filter(C7_tpm < E14_tpm) %>% filter(C7_tpm==0))[1]
print("Down:")
dim(Data%>%filter(annotation=="intron") %>% filter(C7_tpm < E14_tpm) %>% filter(C7_tpm>0))[1]


print("###################################")
print("Comparing A1 and E14: unknown")
print("On:")
dim(Data%>%filter(annotation=="unknown") %>% filter(A1_tpm > E14_tpm) %>% filter(E14_tpm==0))[1]
print("Up:")
dim(Data%>%filter(annotation=="unknown") %>% filter(A1_tpm > E14_tpm) %>% filter(E14_tpm>0))[1]
print("Off:")
dim(Data%>%filter(annotation=="unknown") %>% filter(A1_tpm < E14_tpm) %>% filter(A1_tpm==0))[1]
print("Down:")
dim(Data%>%filter(annotation=="unknown") %>% filter(A1_tpm < E14_tpm) %>% filter(A1_tpm>0))[1]
print("###################################")
print("Comparing C7 and E14: unknown")
print("On:")
dim(Data%>%filter(annotation=="unknown") %>% filter(C7_tpm > E14_tpm) %>% filter(E14_tpm==0))[1]
print("Up:")
dim(Data%>%filter(annotation=="unknown") %>% filter(C7_tpm > E14_tpm) %>% filter(E14_tpm>0))[1]
print("Off:")
dim(Data%>%filter(annotation=="unknown") %>% filter(C7_tpm < E14_tpm) %>% filter(C7_tpm==0))[1]
print("Down:")
dim(Data%>%filter(annotation=="unknown") %>% filter(C7_tpm < E14_tpm) %>% filter(C7_tpm>0))[1]

print("###################################")
print("Comparing A1, C7, and E14: promoter")
print("On:")
dim(Data%>%filter(annotation=="promoter") %>% filter(A1_tpm > E14_tpm & C7_tpm > E14_tpm) %>% filter(E14_tpm==0))[1]
print("Up:")
dim(Data%>%filter(annotation=="promoter") %>% filter(A1_tpm > E14_tpm & C7_tpm > E14_tpm) %>% filter(E14_tpm>0))[1]
print("Off:")
dim(Data%>%filter(annotation=="promoter") %>% filter(A1_tpm < E14_tpm & C7_tpm < E14_tpm) %>% filter(A1_tpm==0 & C7_tpm==0))[1]
print("Down:")
dim(Data%>%filter(annotation=="promoter") %>% filter(A1_tpm < E14_tpm & C7_tpm < E14_tpm) %>% filter(A1_tpm>0 & C7_tpm>0))[1]

print("###################################")
print("Comparing A1, C7, and E14: exon")
print("On:")
dim(Data%>%filter(annotation=="exon") %>% filter(A1_tpm > E14_tpm & C7_tpm > E14_tpm) %>% filter(E14_tpm==0))[1]
print("Up:")
dim(Data%>%filter(annotation=="exon") %>% filter(A1_tpm > E14_tpm & C7_tpm > E14_tpm) %>% filter(E14_tpm>0))[1]
print("Off:")
dim(Data%>%filter(annotation=="exon") %>% filter(A1_tpm < E14_tpm & C7_tpm < E14_tpm) %>% filter(A1_tpm==0 & C7_tpm==0))[1]
print("Down:")
dim(Data%>%filter(annotation=="exon") %>% filter(A1_tpm < E14_tpm & C7_tpm < E14_tpm) %>% filter(A1_tpm>0 & C7_tpm>0))[1]

print("###################################")
print("Comparing A1, C7, and E14: intron")
print("On:")
dim(Data%>%filter(annotation=="intron") %>% filter(A1_tpm > E14_tpm & C7_tpm > E14_tpm) %>% filter(E14_tpm==0))[1]
print("Up:")
dim(Data%>%filter(annotation=="intron") %>% filter(A1_tpm > E14_tpm & C7_tpm > E14_tpm) %>% filter(E14_tpm>0))[1]
print("Off:")
dim(Data%>%filter(annotation=="intron") %>% filter(A1_tpm < E14_tpm & C7_tpm < E14_tpm) %>% filter(A1_tpm==0 & C7_tpm==0))[1]
print("Down:")
dim(Data%>%filter(annotation=="intron") %>% filter(A1_tpm < E14_tpm & C7_tpm < E14_tpm) %>% filter(A1_tpm>0 & C7_tpm>0))[1]

print("###################################")
print("Comparing A1, C7, and E14: unknown")
print("On:")
dim(Data%>%filter(annotation=="unknown") %>% filter(A1_tpm > E14_tpm & C7_tpm > E14_tpm) %>% filter(E14_tpm==0))[1]
print("Up:")
dim(Data%>%filter(annotation=="unknown") %>% filter(A1_tpm > E14_tpm & C7_tpm > E14_tpm) %>% filter(E14_tpm>0))[1]
print("Off:")
dim(Data%>%filter(annotation=="unknown") %>% filter(A1_tpm < E14_tpm & C7_tpm < E14_tpm) %>% filter(A1_tpm==0 & C7_tpm==0))[1]
print("Down:")
dim(Data%>%filter(annotation=="unknown") %>% filter(A1_tpm < E14_tpm & C7_tpm < E14_tpm) %>% filter(A1_tpm>0 & C7_tpm>0))[1]

##################################################################
data_exon_up   <- Data%>%filter(annotation=="exon") %>% filter(A1_tpm > E14_tpm & C7_tpm > E14_tpm) %>% filter(E14_tpm>0)
data_intron_up <- Data%>%filter(annotation=="intron") %>% filter(A1_tpm > E14_tpm & C7_tpm > E14_tpm) %>% filter(E14_tpm>0)

data_up <- rbind(data_exon_up, data_intron_up) %>% mutate(log2FC=0.5*(log2(A1_tpm/E14_tpm) + log2(C7_tpm/E14_tpm)))  %>% arrange(log2FC)

#write.table(data_up, file="Spurious_TSS.bed", sep="\t", quote=F, row.names=F, col.names=T)

data_exon_down   <- Data%>%filter(annotation=="exon") %>% filter(A1_tpm < E14_tpm & C7_tpm < E14_tpm) %>% filter(A1_tpm>0)
data_intron_down <- Data%>%filter(annotation=="intron") %>% filter(A1_tpm < E14_tpm & C7_tpm < E14_tpm)

data_down <- rbind(data_exon_down, data_intron_down) %>% mutate(log2FC=0.5*(log2(A1_tpm/E14_tpm) + log2(C7_tpm/E14_tpm))) %>% arrange(desc(log2FC))

pdf("up_and_down.pdf", width=6, height=11)
q1 <- qplot(1:2232, data_up$log2FC)
q2 <- qplot(1:377,  data_down$log2FC)
grid.arrange(q1, q2, nrow=2, ncol=1)
dev.off()

UP <- data.frame( "Expression"=c( data_up$E14_tpm, data_up$A1_tpm, data_up$C7_tpm), "Sample"=factor(rep(c("E14", "A1", "C7"),each=2232), levels=c("E14", "A1", "C7")))
DOWN <- data.frame( "Expression"=c( data_down$E14_tpm, data_down$A1_tpm, data_down$C7_tpm), "Sample"=factor(rep(c("E14", "A1", "C7"),each=377),levels=c("E14", "A1", "C7")))

ggplot(UP, aes(x=Sample, y=Expression, col=Sample))+geom_boxplot(notch=TRUE) + ylim(0,25)
ggplot(DOWN, aes(x=Sample, y=Expression, col=Sample))+geom_boxplot(notch=TRUE) + ylim(0,45)

library(beanplot)
#beanplot(Expression ~ Sample, data=UP, col = list("grey50", "red", "blue"), beanlines="mean", maxstripline=0.03, yaxs="i")
#beanplot(Expression ~ Sample, data=DOWN, col = list("grey50", "red", "blue"), beanlines="mean", maxstripline=0.03, yaxs="i")

######
source("/usr/people/bioc1387/Scripts/LoadingExpressionGeneGroup.R")

##### PWWP2A taragets
PWWP2A_taragets <- read.table("/usr/people/bioc1387/Project/Pwwp2a_Tianyi/Bam/Peak2/PWWP2A_targets2_GeneSymbol", header=F)
PWWP2A_non_taragets <- read.table("/usr/people/bioc1387/Project/Pwwp2a_Tianyi/Bam/Peak2/PWWP2A_non_targets2.GeneSymbol", header=F)

TakeGenelist <- function(data){
	data$genes2=""
	N=dim(data)
	for(i in 1:N[1]){
		data[i, N[2] ] = unlist(strsplit( toString( data[i,7]), ";"))[1]
		}
	return(data$genes2)
	}

write.table(unique(TakeGenelist(data_up)), file="Spurious_TSS_Gene.txt", sep="\t", quote=F, col.names=F, row.names=F)

