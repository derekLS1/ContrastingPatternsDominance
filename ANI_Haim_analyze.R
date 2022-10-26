#transfer to MAC
scp dlundberg@burrito.eb.local:/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/ANI/Sphingo_noIBVSS_ANI.ani.matrix /Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS

#here is a list of the 174 ASV1 Sphingomonas genomes:
# S282H133 S49H113 S128H113 S220H113 S227H113 S231H113 S265H133 S212H113 S363H113 S301H113 S300H113 S354H113 S266H113 S377H113 S222H113 S249H113 S219H113 S283H113 S35H113 S159H133 S105H113 S308H113 S225H113 S226H113 S130H113 S233H113 S241H133 S373H113 S365H113 S237H113 S314H113 S243H113 S310H113 S234H113 S232H113 S328H113 S257H113 S258H133 S253H133 S271H133 S61H113 S247H133 S235H133 S133H113 S250H133 S197H133 S236H133 S254H133 S351H113 S204H133 S151H133 S364H113 S358H113 S250H113 S244H113 S241H113 S199H133 S361H113 S216H113 S260H133 S273H133 S160H133 S148H133 S152H133 S257H133 S161H133 S161H113 S266H133 S239H133 S193H133 S245H133 S206H133 S249H133 S208H113 S228H133 S214H133 S156H133 S108H113 S153H133 S166H133 S371H113 S200H113 S204H113 S222H133 S231H133 S110H113 S116H113 S166H113 S164H113 S158H113 S160H113 S27H113 S21H113 S278H113 S304H113 S254H113 S59H113 S299H113 S298H113 S274H113 S76H113 S293H113 S119H113 S18H113 S280H113 S276H113 S245H113 S246H113 S235H113 S236H113 S238H113 S218H113 S279H113 S127H113 S356H113 S242H113 S284H133 S6H113 S353H113 S275H113 S281H113 S282H113 S284H113 S163H133 S270H133 S149H133 S111H113 S155H133 S263H133 S218H133 S33H113 S322H113 S117H113 S229H113 S268H113 S261H113 S262H113 S269H113 S271H113 S265H113 S264H113 S267H113 S258H113 S374H113 S259H113 S230H113 S215H113 S213H113 S309H113 S228H113 S252H113 S223H113 S247H113 S253H113 S251H113 S292H113 S294H113 S357H113 S289H113 S285H113 S287H113 S256H113 S286H113 S255H113 S239H113 S240H113 S232H133 S355H113 S44H113 S57H113 S118H113 S25H113 S273H113 S277H113


date=format(Sys.Date(), format="%Y%m%d")
library(reshape)


Sani<-scan("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS/Sphingo_noIBVSS_ANI.ani.matrix", sep="\t", quote="", fill=TRUE, what="character")
#Sani<-scan("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS_ASV1/Sphingo_ASV1_ANI.ani.matrix", sep="\t", quote="", fill=TRUE, what="character")
c(grep("S.*H.*", Sani, value=TRUE), grep("GCF.*", Sani, value=TRUE))->allgenomes
matrix(ncol=length(allgenomes)+1, nrow=length(allgenomes)+1)->Sanimat
r=1
c=2
for(i in Sani){
	if(length(which(i==allgenomes)>0)>0){
		#print("newrow")
		r+1->r
		c=1
	}	
	Sanimat[r, c]=Sanimat[c, r]=i
	c+1->c
}
Sanimat[1,2:ncol(Sanimat)]->ScolNames
Sanimat[2:nrow(Sanimat),1]->SrowNames
apply(Sanimat[2:nrow(Sanimat), 2:ncol(Sanimat)], 2, as.numeric)->Sanimat
gsub(".*\\/", "", ScolNames)->ScolNames
gsub(".*\\/", "", SrowNames)->SrowNames
ScolNames->colnames(Sanimat)
SrowNames->rownames(Sanimat)
Sanimat[is.na(Sanimat)]<-100
#USES A PREVIOUSLY DEFINED CORE TREE TO SELECT WHICH ONES TO SHOW
#Here, must subset for ASV1. Run "vis_panXoutput_core_genome_tree.R" with SphASV1 as base folder
gsub("GCF_", "GCF@", colnames(Sanimat))->SanimatMatcher
gsub("_S", "@S", SanimatMatcher)->SanimatMatcher
gsub("_.*", "", SanimatMatcher)->SanimatMatcher
gsub("-S", "@S", core$tip.label)->coretipMatcher
gsub("GCF_", "GCF@", coretipMatcher)->coretipMatcher
gsub("_.*", "", coretipMatcher)->coretipMatcher
gsub("GCF@000971055.*", "GCF@000971055", coretipMatcher)->coretipMatcher
gsub("GCF@000971055.*", "GCF@000971055", SanimatMatcher)->SanimatMatcher
match(coretipMatcher, SanimatMatcher)->Sanimat_subset
Sanimat[Sanimat_subset,Sanimat_subset]->Sanimat

Pani<-scan("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Pseudo_runs/Pseudo_ASV1/Pseudo_ASV1_ANI.ani.matrix", sep="\t", quote="", fill=TRUE, what="character")
c(grep("S.*H.*", Pani, value=TRUE), grep("GCF.*", Pani, value=TRUE), grep("p.*", Pani, value=TRUE))->allgenomes
matrix(ncol=length(allgenomes)+1, nrow=length(allgenomes)+1)->Panimat
r=1
c=2
for(i in Pani){
	if(length(which(i==allgenomes)>0)>0){
		print("newrow")
		r+1->r
		c=1
	}	
	Panimat[r, c]=Panimat[c, r]=i
	c+1->c
}
Panimat[1,2:ncol(Panimat)]->PcolNames
Panimat[2:nrow(Panimat),1]->ProwNames
apply(Panimat[2:nrow(Panimat), 2:ncol(Panimat)], 2, as.numeric)->Panimat
PcolNames->colnames(Panimat)
ProwNames->rownames(Panimat)
Panimat[is.na(Panimat)]<-100


#USES A PREVIOUSLY DEFINED CORE TREE TO SELECT WHICH ONES TO SHOW
#Here, must subset for ASV1. Run "vis_panXoutput_core_genome_tree.R" with SphASV1 as base folder
gsub("GCF_", "GCF@", colnames(Panimat))->PanimatMatcher
gsub("_S", "@S", PanimatMatcher)->PanimatMatcher
gsub("_.*", "", PanimatMatcher)->PanimatMatcher
gsub("-S", "@S", core$tip.label)->coretipMatcher
gsub("GCF_", "GCF@", coretipMatcher)->coretipMatcher
gsub("_.*", "", coretipMatcher)->coretipMatcher
gsub("GCF@000971055.*", "GCF@000971055", coretipMatcher)->coretipMatcher
gsub("GCF@000971055.*", "GCF@000971055", PanimatMatcher)->PanimatMatcher
match(coretipMatcher, PanimatMatcher)->Panimat_subset
Panimat[Panimat_subset,Panimat_subset]->Panimat

#SphOrder<-as.matrix(read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS_ASV1/strain_tree_order.txt", quote="", sep="\t"))
#PsOrder<-as.matrix(read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Pseudo_runs/Pseudo_ASV1/strain_tree_order.txt", quote="", sep="\t"))
#gsub("-S", "_S", SphOrder)->SphOrder
#gsub("GCF_", "GCF", SphOrder)->SphOrder
#gsub("_.*", "", SphOrder)->SphOrder
#gsub("GCF", "GCF_", SphOrder)->SphOrder

#gsub("GCF_", "GCF", ScolNames)->ScolNames.tmp
#gsub("_.*", "", ScolNames.tmp)->ScolNames.tmp
#gsub("GCF", "GCF_", ScolNames.tmp)->ScolNames.tmp
#Sanimat[match(SphOrder, ScolNames.tmp),match(SphOrder, ScolNames.tmp)]->Sanimat
		
#PcolNames->PcolNames.tmp
#Panimat[match(PsOrder, PcolNames.tmp),match(PsOrder, PcolNames.tmp)]->Panimat

				




		
pal1 <- c(colorRampPalette(c("black","#27004E", "#340069","#6A1B9A", "#CB018E", "#EB6E99", "#FFBC7A", "#FFFFFF"))(n=9))
	#pal1.1[length(pal1.1):1]->pal1

pdf(file=paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS_ASV1/ANI_heatmap_LEGEND_", date, ".pdf", sep="", collapse=""), width = 3, height = 3, useDingbats=FALSE)
image(1:length(pal1), 1, as.matrix(1:length(pal1)), 
      col=pal1, xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
dev.off()  
pdf(file=paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS_ASV1/ANI_heatmap_", date, ".pdf", sep="", collapse=""), width = 3, height = 3, useDingbats=FALSE)
heatmap.2(Sanimat, density.info = "none",
                    trace = "none", dendrogram="none",
                    col = pal1, Colv=FALSE, Rowv=FALSE, key=FALSE, 
                    key.title = NA, key.xlab = NA, sepcolor=NA,
                    breaks = c(80, 82.5, 85, 87.5, 90, 92.5, 95, 97.5, 99.9, 100),
                    labCol = F, labRow=F, sepwidth=c(0,0))
dev.off()
    
    
pdf(file=paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Pseudo_runs/Pseudo_ASV1/ANI_heatmap_", date, ".pdf", sep="", collapse=""), width = 3, height = 3, useDingbats=FALSE)
heatmap.2(Panimat, density.info = "none",
                    trace = "none", dendrogram="none",
                    col = pal1, Colv=FALSE, Rowv=FALSE, key=FALSE, 
                    key.title = NA, key.xlab = NA, sepcolor=NA,
                    breaks = c(80, 82.5, 85, 87.5, 90, 92.5, 95, 97.5, 99.9, 100),
                    labCol = F, labRow=F, sepwidth=c(0,0))
dev.off() 
                
                
    
    
#compare intra-vs-inter species ANI
#LOAD METADATA
metadata=read.table(file="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Single_genomes/core_genome_tree_metadata.txt", comment.char="", quote = "")
as.matrix(metadata[1,])->names(metadata)
metadata[2:nrow(metadata),]->metadata


#SPHINGOMONAS
gsub("GCF_", "GCF@", colnames(Sanimat))->SanimatMatcher
gsub(".*_S", "S", SanimatMatcher)->SanimatMatcher
gsub("_.*", "", SanimatMatcher)->SanimatMatcher
gsub("GCF@000971055.*", "GCF@000971055", SanimatMatcher)->SanimatMatcher
gsub("GCF_", "GCF@", metadata$SeqID)->metadataMatcher
gsub("_.*", "", metadataMatcher)->metadataMatcher
gsub("GCF@000971055.*", "GCF@000971055", metadataMatcher)->metadataMatcher
as.matrix(metadata$host_species[match(SanimatMatcher,metadataMatcher)])->sanimat_plantlabels
	AthalianaANI=Sanimat[which(sanimat_plantlabels=="Athaliana"),which(sanimat_plantlabels=="Athaliana")]
    	#the below gets rid of duplicate calcualtions because its a square matrix against itself
    	for(r in 1:nrow(AthalianaANI)){
    		AthalianaANI[r,r:ncol(AthalianaANI)]<-NA
    	}
    NotAthalianaANI=Sanimat[which(sanimat_plantlabels=="Athaliana"),which(sanimat_plantlabels!="Athaliana")]                
        #for(r in 1:nrow(NotAthalianaANI)){
        #	if(r<=ncol(NotAthalianaANI)){
    	#		NotAthalianaANI[r,r:ncol(NotAthalianaANI)]<-NA
    	#	}
    	#}
    as.data.frame(rbind(
				cbind(rep("AthalianaANI", length(as.vector(AthalianaANI))), as.vector(AthalianaANI)),
				cbind(rep("NotAthalianaANI", length(as.vector(NotAthalianaANI))), as.vector(NotAthalianaANI))
				))->histANI
				
				as.numeric(as.character(histANI[,2]))->histANI[,2]
				histANI[,1]=factor(histANI[,1], levels=c("AthalianaANI", "NotAthalianaANI"))	
#75 athalianas from H133 vs. 71 non athalianas from H133
#de2e26	#000000	#000000
DerekColors=c("#A9CF38","#6050A1")
pdf(file=paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS/Athaliana_V_Not_density_", date, ".pdf", sep="", collapse=""), width = 5, height = 2.5, useDingbats=FALSE)
ggplot(histANI, aes(x=histANI[,2], fill=histANI[,1], color=histANI[,1])) +  #geom_histogram(bins=50, position="dodge")
	geom_density(alpha=.2, adjust=.25) + theme_classic() + 
	theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15)) + 
	theme(axis.line = element_line(colour = "black", size = .15), axis.ticks = element_line(colour = "black", size = .15))+
	scale_color_manual(values=DerekColors) + scale_fill_manual(values=DerekColors) + 
	#scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), limits = c(0,.5)) + 
	scale_x_continuous(limits = c(80,100)) + 
	theme(axis.text.x = element_text(size=15, color="black"),axis.text.y = element_text(size=15, color="black")) + 
	theme(axis.line = element_line(color = "black", size = .15), axis.ticks = element_line(color = "black", size = .15))
dev.off()

AthalianaANIhist<-hist(AthalianaANI, probability=T, breaks=seq(82, 100, by=.02))
	AthalianaANIhist$density
NotAthalianaANIhist<-hist(NotAthalianaANI, probability=T, breaks=seq(82, 100, by=.02))
	NotAthalianaANIhist$density
	
	mean(
	wilcox.test(NotAthalianaANIhist$density, AthalianaANIhist$density)

par(mfrow=c(2,1))
hist(AthalianaANI, probability=T, breaks=seq(82, 100, by=.02))
hist(NotAthalianaANI, probability=T, breaks=seq(82, 100, by=.02))


t.test(sort(AthalianaANI),sort(NotAthalianaANI))

mean(NotAthalianaANI)
mean(sort(AthalianaANI))


t.test(sort(AthalianaANI), sort(NotAthalianaANI), alternative="greater")





boxplot(sort(AthalianaANI), sort(NotAthalianaANI), outline=FALSE, ylim=c(80, 100))
	points(jitter(rep(1, length(sort(AthalianaANI))), amount=.2), sort(AthalianaANI), pch=".")
	points(jitter(rep(2, length(sort(NotAthalianaANI))), amount=.2), sort(NotAthalianaANI), pch=".")



ks.test(AthalianaANIhist$density, NotAthalianaANIhist$density)

ks.test(AthalianaANIhist$density, test3)

a<-qqplot(AthalianaANI, NotAthalianaANI)
ks.test(a$x, a$y)
ks.test(a$x, jitter(a$x, amount=.6))
par(mfrow=c(1,1))
qqplot(a$x, a$y)
qqplot(a$x, jitter(a$x, amount=.6))
plot(a$x)
plot(a$y)

#PSEUDOMONAS
gsub("GCF_", "GCF@", colnames(Panimat))->PanimatMatcher
gsub(".*_S", "S", PanimatMatcher)->PanimatMatcher
gsub("_.*", "", PanimatMatcher)->PanimatMatcher
gsub("GCF@000971055.*", "GCF@000971055", PanimatMatcher)->PanimatMatcher
gsub("GCF_", "GCF@", metadata$SeqID)->metadataMatcher
gsub("_.*", "", metadataMatcher)->metadataMatcher
gsub("GCF@000971055.*", "GCF@000971055", metadataMatcher)->metadataMatcher
as.matrix(metadata$host_species[match(PanimatMatcher,metadataMatcher)])->Panimat_plantlabels
as.matrix(metadata$plantID[match(PanimatMatcher,metadataMatcher)])->Panimat_samplelabels
cbind(PanimatMatcher, Panimat_plantlabels, Panimat_samplelabels)->Panimat_key
	AthalianaANI=Panimat[which(Panimat_plantlabels=="Athaliana"),which(Panimat_plantlabels=="Athaliana")]
    	#the below gets rid of duplicate calcualtions because its a square matrix against itself
    	for(r in 1:nrow(AthalianaANI)){
    		AthalianaANI[r,r:ncol(AthalianaANI)]<-NA
    	}
    NotAthalianaANI=Panimat[which(Panimat_plantlabels=="Athaliana"),which(Panimat_plantlabels!="Athaliana")]                
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


dx=Panimat
dx[which(dx<99.9)]=0
write.table(dx, "/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Pseudo_runs/Pseudo_ASV1/Panimat.txt")
dx=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Pseudo_runs/Pseudo_ASV1/Panimat2.txt")
apply(dx, 1, max)->Pse_strainIDs

#compare all to all, no duplicated samples
ANI=as.matrix(dx)
ANI=melt(ANI)
ANI=ANI[ANI[,3]>0,]
	ANI=cbind(ANI, Panimat_key[match(ANI[,1], Panimat_key[,1]),2:3])
	ANI=cbind(ANI, Panimat_key[match(ANI[,2], Panimat_key[,1]),2:3])

#make row strain by sample combos
keeps=vector(length=0)
sample_strain=paste(as.matrix(ANI[,5]), as.matrix(ANI[,3]), sep="_")
for(u in unique(sample_strain)){
	print("new")
	repr=names(sort(table(matrix(as.matrix(ANI[which(sample_strain==u),1:2]))),decreasing=TRUE)[1])
	print(repr)
	keeps=c(keeps, repr)
}
ANI=ANI[match(as.matrix(ANI[,1]), keeps, nomatch=0)>0,]
ANI=ANI[match(as.matrix(ANI[,2]), keeps, nomatch=0)>0,]
	
Panimat_test=Panimat
Panimat_test=Panimat_test[match(rownames(Panimat_test), keeps, nomatch=0)>0, match(colnames(Panimat_test), keeps, nomatch=0)>0]
	#removes redundancy of duplicate comparisons and 100% comparisons
	for(r in 1:nrow(Panimat_test)){
    		Panimat_test[r,r:ncol(Panimat_test)]<-NA
    }
    heatmap.2(Panimat_test, trace="none")
gsub("GCF_", "GCF@", colnames(Panimat_test))->PanimatMatcher
gsub(".*_S", "S", PanimatMatcher)->PanimatMatcher
gsub("_.*", "", PanimatMatcher)->PanimatMatcher
gsub("GCF@000971055.*", "GCF@000971055", PanimatMatcher)->PanimatMatcher
gsub("GCF_", "GCF@", metadata$SeqID)->metadataMatcher
gsub("_.*", "", metadataMatcher)->metadataMatcher
gsub("GCF@000971055.*", "GCF@000971055", metadataMatcher)->metadataMatcher
as.matrix(metadata$host_species[match(PanimatMatcher,metadataMatcher)])->Panimat_plantlabels
as.matrix(metadata$plantID[match(PanimatMatcher,metadataMatcher)])->Panimat_samplelabels
cbind(PanimatMatcher, Panimat_plantlabels, Panimat_samplelabels)->Panimat_key
	AthalianaANI=Panimat[which(Panimat_plantlabels=="Athaliana"),which(Panimat_plantlabels=="Athaliana")]
    	#the below gets rid of duplicate calcualtions because its a square matrix against itself
    	for(r in 1:nrow(AthalianaANI)){
    		AthalianaANI[r,r:ncol(AthalianaANI)]<-NA
    	}
    NotAthalianaANI=Panimat[which(Panimat_plantlabels=="Athaliana"),which(Panimat_plantlabels!="Athaliana")]                
    as.data.frame(rbind(
				cbind(rep("AthalianaANI", length(as.vector(AthalianaANI))), as.vector(AthalianaANI)),
				cbind(rep("NotAthalianaANI", length(as.vector(NotAthalianaANI))), as.vector(NotAthalianaANI))
				))->histANI
				
				as.numeric(as.character(histANI[,2]))->histANI[,2]
				histANI[,1]=factor(histANI[,1], levels=c("AthalianaANI", "NotAthalianaANI"))	




	
cbind(ANI, duplicated(ANI[as.character(ANI[,4])==as.character(ANI[,5]),4]))

AthalianaANIhist<-hist(AthalianaANI, probability=T, breaks=seq(82, 100, by=2))
	AthalianaANIhist$density
NotAthalianaANIhist<-hist(NotAthalianaANI, probability=T, breaks=seq(82, 100, by=2))
	NotAthalianaANIhist$density

par(mfrow=c(2,1))
hist(AthalianaANI, probability=T, breaks=seq(82, 100, by=.1))
hist(NotAthalianaANI, probability=T, breaks=seq(82, 100, by=.1))

wilcox.test(sort(AthalianaANI), sort(NotAthalianaANI))
wilcox.test(sort(AthalianaANI), sort(AthalianaANI))

median(sort(AthalianaANI))
median(NotAthalianaANI)
par(mfrow=c(1,1))
boxplot(sort(AthalianaANI), sort(NotAthalianaANI), outline=FALSE, ylim=c(96, 100))
	points(jitter(rep(1, length(sort(AthalianaANI))), amount=.2), sort(AthalianaANI), pch=".")
	points(jitter(rep(2, length(sort(NotAthalianaANI))), amount=.2), sort(NotAthalianaANI), pch=".")



pdf(file=paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Pseudo_runs/Pseudo_ASV1/Athaliana_V_Not_density_", date, ".pdf", sep="", collapse=""), width = 5, height = 2.5, useDingbats=FALSE)
ggplot(histANI, aes(x=histANI[,2], fill=histANI[,1], color=histANI[,1])) +  
 	geom_density(alpha=.2, adjust=.25) + 
 	theme_classic() + 
	theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15)) + 
	theme(axis.line = element_line(colour = "black", size = .15), axis.ticks = element_line(colour = "black", size = .15))+
	scale_color_manual(values=DerekColors) + scale_fill_manual(values=DerekColors) + 
	#scale_y_continuous(breaks = c(0, 2.5, 5), labels = c(0, 2.5, 5), limits = c(0,7)) + 
	scale_x_continuous(limits = c(96,100)) + 
	theme(axis.text.x = element_text(size=15, color="black"),axis.text.y = element_text(size=15, color="black")) + 
	theme(axis.line = element_line(color = "black", size = .15), axis.ticks = element_line(color = "black", size = .15))
dev.off()


####### finds how many are above 99.99

countPanimat=matrix(ncol=ncol(Panimat), nrow=nrow(Panimat), data=0)
rownames(countPanimat)<-rownames(Panimat)
colnames(countPanimat)<-colnames(Panimat)
#number of isolates >= 99.9 to another
	countPanimat[which(Panimat<100 & Panimat>=99.9)]<-1
	length(which(apply(countPanimat, 1, max)==1))
	isolates99.9=rownames(countPanimat)[which(apply(countPanimat, 1, max)==1)]




length(which(Sanimat<100 & Sanimat>=99.9))/2



