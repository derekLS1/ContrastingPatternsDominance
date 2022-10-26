library(gplots)

#Pseudomonas
genus="Pseudomonas"
if(genus=="Pseudomonas"){
	fullmatrixP=read.table("/Users/dlundberg/Documents/abt6/Pratchaya/Bulk_metagenomes/NonRedundantMapping/Pseudomonas/fullmatrix.txt", sep="\t")
}
if(genus=="Sphingomonas"){
	fullmatrixP=read.table("/Users/dlundberg/Documents/abt6/Pratchaya/Bulk_metagenomes/NonRedundantMapping/Sphingomonas/fullmatrix.txt", sep="\t")
}

as.matrix(fullmatrixP)->fullmatrixP
fullmatrixP[1,]->samplenamesP
fullmatrixP[2:nrow(fullmatrixP),]->fullmatrixP
apply(fullmatrixP, 2, as.numeric)->fullmatrixP
samplenamesP->colnames(fullmatrixP)
samplenamesP->rownames(fullmatrixP)
		
100*(1-fullmatrixP)->fullmatrixP


#bad Pseudomonas genomes by N50, or not in the core genome tree
if(genus=="Pseudomonas"){ keep=grep("p22.D2|p26.F7|p26.F6|p22.D5|p22.D7|p23.D7|p25.C11|S64H133|p26.D2", colnames(fullmatrixP), invert=TRUE) }
if(genus=="Sphingomonas"){ keep=grep("S208H133|S136H133|S101H133|S99H133|S127H133|S105H133|S104H133|S97H133|S147H133|S98H133|S216H133|S112H133|S279H133|S275H133|S285H133|S201H133|S107H133|S100H133|S74H113|S167H113|S372H113|S168H113|S170H113|S165H113|S169H113|S172H113|S171H113|S67H113", colnames(fullmatrixP), invert=TRUE) }






fullmatrixP[keep,keep]->fullmatrixP

print(dim(fullmatrixP))
#try deleting anything that is more than 99% similar to something else
s=1
while(s < nrow(fullmatrixP)){
	#find all non-self matches, and delete those from the matrix. 
	kill=setdiff(which(fullmatrixP[s,]>97.5), s)
	if(length(kill)>0){
		keep=setdiff(1:nrow(fullmatrixP), kill)
		fullmatrixP[keep, keep]->fullmatrixP
	}
	s+1->s
}

print(dim(fullmatrixP))


pal1 <- c(colorRampPalette(c("black","#27004E", "#340069","#6A1B9A", "#CB018E", "#EB6E99", "#FFBC7A", "#FFFFFF"))(n=9))
image(1:length(pal1), 1, as.matrix(1:length(pal1)), 
      col=pal1, xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
if(genus=="Pseudomonas"){
	pdf(file="/Users/dlundberg/Documents/abt6/Pratchaya/Bulk_metagenomes/NonRedundantMapping/Pseudomonas/fullmatrixP.pdf", width=10, height=10, useDingbats=FALSE)
	write.table(rownames(fullmatrixP), "/Users/dlundberg/Documents/abt6/Pratchaya/Bulk_metagenomes/NonRedundantMapping/Pseudomonas/condensed_Pref.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

}
if(genus=="Sphingomonas"){
	pdf(file="/Users/dlundberg/Documents/abt6/Pratchaya/Bulk_metagenomes/NonRedundantMapping/Sphingomonas/fullmatrixS.pdf", width=10, height=10, useDingbats=FALSE)
	write.table(rownames(fullmatrixP), "/Users/dlundberg/Documents/abt6/Pratchaya/Bulk_metagenomes/NonRedundantMapping/Sphingomonas/condensed_Sref.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
}
heatmap.2(fullmatrixP, density.info = "none",
                    trace = "none", dendrogram="none",
                    col = pal1, Colv=TRUE, Rowv=TRUE, key=FALSE, 
                    key.title = NA, key.xlab = NA, sepcolor=NA,
                    breaks = c(80, 82.5, 85, 87.5, 90, 92.5, 95, 97.5, 99, 100),
                    labCol = F, labRow=rownames(fullmatrixP), cexRow=0.1, sepwidth=c(0,0))
dev.off()
