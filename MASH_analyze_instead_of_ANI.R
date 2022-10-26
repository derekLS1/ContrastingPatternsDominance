	
library(RColorBrewer)
library(pheatmap)
library(viridis)
library(gplots)

#1) run WGS, otu99=FALSE
#2) run core, otu99=FALSE
#3) run core, otu99=TRUE
#1) run WGS, otu99=FALSE
#must run WGS first
matrixType="WGS"    #core or WGS 
otu99=FALSE
fake_core_pseudomonas=FALSE


date=format(Sys.Date(), format="%Y%m%d")

#Mash_Sphingomonas
#LOAD dist matrix
if(matrixType=="WGS"){ fullmatrix_Sph<-read.table(file="/Users/dlundberg/Documents/abt6/Pratchaya/PanX/MASH/fullmatrix.txt", sep="\t") }

as.matrix(fullmatrix_Sph)->fullmatrix_Sph
fullmatrix_Sph[1,]->samplenames_Sph
fullmatrix_Sph[2:nrow(fullmatrix_Sph),]->fullmatrix_Sph
apply(fullmatrix_Sph, 2, as.numeric)->fullmatrix_Sph
	#gsub("GCF_000259955.1_ASM25995v1_genomic.fna", "GCF_000259955.1_ASM25995v1_genomic_Sphingomonas_melonis_FR1", samplenames_Sph)->samplenames_Sph
	#gsub("GCF_000419565.1_ASM41956v1_genomic.fna", "GCF_000419565.1_ASM41956v1_genomic_Sphingomonas_wittichii_DP58", samplenames_Sph)->samplenames_Sph
gsub("_coreAln.fasta", "", samplenames_Sph)->samplenames_Sph
samplenames_Sph->colnames(fullmatrix_Sph)
samplenames_Sph->rownames(fullmatrix_Sph)
fullmatrix_Sph->reads



#USES A PREVIOUSLY DEFINED CORE TREE TO SELECT WHICH ONES TO SHOW
#Here, must subset for ASV1. Run "vis_panXoutput_core_genome_tree.R" with SphASV1 as base folder
gsub("GCF_", "GCF@", colnames(reads))->readsMatcher
gsub("_S", "@S", readsMatcher)->readsMatcher
gsub("_.*", "", readsMatcher)->readsMatcher
gsub("-S", "@S", core$tip.label)->coretipMatcher
gsub("GCF_", "GCF@", coretipMatcher)->coretipMatcher
gsub("_.*", "", coretipMatcher)->coretipMatcher
gsub("GCF@000971055.*", "GCF@000971055", coretipMatcher)->coretipMatcher
gsub("GCF@000971055.*", "GCF@000971055", readsMatcher)->readsMatcher
match(coretipMatcher, readsMatcher)->reads_subset
reads[reads_subset,reads_subset]->reads

#convert to ANI-like similarity
100*(1-reads)->reads

		
pal1 <- c(colorRampPalette(c("black","#27004E", "#340069","#6A1B9A", "#CB018E", "#EB6E99", "#FFBC7A", "#FFFFFF"))(n=9))
	#pal1.1[length(pal1.1):1]->pal1

if(genus=="Sphingomonas"){
pdf(file=paste("/Users/dlundberg/Documents/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS/MASH_ASV1_heatmap_LEGEND", date, ".pdf", sep="", collapse=""), width = 3, height = 3, useDingbats=FALSE)
image(1:length(pal1), 1, as.matrix(1:length(pal1)), 
      col=pal1, xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
dev.off()
pdf(file=paste("/Users/dlundberg/Documents/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS/MASH_ASV1_heatmap", date, ".pdf", sep="", collapse=""), width = 3, height = 3, useDingbats=FALSE)
heatmap.2(reads, density.info = "none",
                    trace = "none", dendrogram="none",
                    col = pal1, Colv=FALSE, Rowv=FALSE, key=FALSE, 
                    key.title = NA, key.xlab = NA, sepcolor=NA,
                    breaks = c(80, 82.5, 85, 87.5, 90, 92.5, 95, 97.5, 99.9, 100),
                    labCol = F, labRow=F, sepwidth=c(0,0))
dev.off()
}    
 
if(genus=="Pseudomonas"){   
pdf(file=paste("/Users/dlundberg/Documents/abt6/Pratchaya/PanX/Pseudo_runs/Pseudo_ASV1/MASH_ASV1_heatmap", date, ".pdf", sep="", collapse=""), width = 3, height = 3, useDingbats=FALSE)
heatmap.2(reads, density.info = "none",
                    trace = "none", dendrogram="none",
                    col = pal1, Colv=FALSE, Rowv=FALSE, key=FALSE, 
                    key.title = NA, key.xlab = NA, sepcolor=NA,
                    breaks = c(80, 82.5, 85, 87.5, 90, 92.5, 95, 97.5, 99.9, 100),
                    labCol = F, labRow=F, sepwidth=c(0,0))
dev.off() 
}         
                
    
    
#compare intra-vs-inter species ANI
#LOAD METADATA
metadata=read.table(file="/Users/dlundberg/Documents/abt6/Pratchaya/Single_genomes/core_genome_tree_metadata.txt", comment.char="", quote = "")
as.matrix(metadata[1,])->names(metadata)
metadata[2:nrow(metadata),]->metadata


gsub("GCF_", "GCF@", colnames(reads))->readsMatcher
gsub(".*_S", "S", readsMatcher)->readsMatcher
gsub("_.*", "", readsMatcher)->readsMatcher
gsub("GCF@000971055.*", "GCF@000971055", readsMatcher)->readsMatcher
gsub("GCF_", "GCF@", metadata$SeqID)->metadataMatcher
gsub("_.*", "", metadataMatcher)->metadataMatcher
gsub("GCF@000971055.*", "GCF@000971055", metadataMatcher)->metadataMatcher
as.matrix(metadata$host_species[match(readsMatcher,metadataMatcher)])->reads_plantlabels
as.matrix(metadata$plantID[match(readsMatcher,metadataMatcher)])->reads_samplelabels
cbind(readsMatcher, reads_plantlabels, reads_samplelabels)->reads_key
	AthalianaANI=reads[which(reads_plantlabels=="Athaliana"),which(reads_plantlabels=="Athaliana")]
    	#the below gets rid of duplicate calcualtions because its a square matrix against itself
    	for(r in 1:nrow(AthalianaANI)){
    		AthalianaANI[r,r:ncol(AthalianaANI)]<-NA
    	}
    NotAthalianaANI=reads[which(reads_plantlabels=="Athaliana"),which(reads_plantlabels!="Athaliana")]                
    as.data.frame(rbind(
				cbind(rep("AthalianaANI", length(as.vector(AthalianaANI))), as.vector(AthalianaANI)),
				cbind(rep("NotAthalianaANI", length(as.vector(NotAthalianaANI))), as.vector(NotAthalianaANI))
				))->histANI
				
				as.numeric(as.character(histANI[,2]))->histANI[,2]
				histANI[,1]=factor(histANI[,1], levels=c("AthalianaANI", "NotAthalianaANI"))	
#75 athalianas from H133 vs. 71 non athalianas from H133
#de2e26	#000000	#000000
DerekColors=c("#A9CF38","#6050A1")
#WG MASH

mean(NotAthalianaANI)
mean(sort(AthalianaANI))
median(NotAthalianaANI)
median(sort(AthalianaANI))
t.test(sort(AthalianaANI), sort(NotAthalianaANI), alternative="greater")
wilcox.test(sort(AthalianaANI), sort(NotAthalianaANI), alternative="greater")

if(genus=="Pseudomonas"){
	pdf(file=paste("/Users/dlundberg/Documents/abt6/Pratchaya/PanX/Pseudo_runs/Pseudo_ASV1/MASH_ASV1_Athaliana_V_Not_density_", date, ".pdf", sep="", collapse=""), width = 5, height = 2.5, useDingbats=FALSE)
	limits=c(96, 100)
}
if(genus=="Sphingomonas"){
	pdf(file=paste("/Users/dlundberg/Documents/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS/MASH_ASV1_Athaliana_V_Not_density", date, ".pdf", sep="", collapse=""), width = 5, height = 2.5, useDingbats=FALSE)
	limits=c(80, 100)
}
ggplot(histANI, aes(x=histANI[,2], fill=histANI[,1], color=histANI[,1])) +  
 	geom_density(alpha=.2, adjust=.25) + 
 	theme_classic() + 
	theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15)) + 
	theme(axis.line = element_line(colour = "black", size = .15), axis.ticks = element_line(colour = "black", size = .15))+
	scale_color_manual(values=DerekColors) + scale_fill_manual(values=DerekColors) + 
	#scale_y_continuous(breaks = c(0, 2.5, 5), labels = c(0, 2.5, 5), limits = c(0,7)) + 
	scale_x_continuous(limits = limits) + 
	theme(axis.text.x = element_text(size=15, color="black"),axis.text.y = element_text(size=15, color="black")) + 
	theme(axis.line = element_line(color = "black", size = .15), axis.ticks = element_line(color = "black", size = .15))
dev.off()



