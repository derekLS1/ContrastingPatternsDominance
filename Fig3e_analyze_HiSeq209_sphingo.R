source("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/scripts/R/functions/microbiome_custom_functions.R")
#install.packages("vegan")
#install.packages("reshape")
library(vegan)
library(reshape)
library(ggplot2)
library(RColorBrewer)

date=format(Sys.Date(), format="%Y%m%d")





	#otutab="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/soil_protection_photos_202107/16S_hamPCR/HiSeq0209_SphingoSoil_20210823.txt"
	otutab="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/soil_protection_photos_202107/16S_hamPCR/HiSeq0209_SphingoSoil_concat_20210917.txt"
	
	#tax="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/soil_protection_photos_202107/16S_hamPCR/HiSeq0209_SphingoSoil.tax"
	tax="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/soil_protection_photos_202107/16S_hamPCR/HiSeq0209_SphingoSoil_concat.tax"
	
	metadataT="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/soil_protection_photos_202107/16S_hamPCR/HiSeq0209_sphingo_metadata.txt"

# name	well	barcode	primer	sample


#load OTU table
	read.countTable(file=otutab)->otureads    #CUSTOM FUNCTION
	
#load metadata
	as.matrix(read.table(metadataT, sep="\t"))->metadata
	metadata[1,]->colnames(metadata)
	metadata[2:nrow(metadata),]->metadata

	#colnames(otureads)=metadata[match(colnames(otureads), metadata[,"sequencefile"]),"sample"]

#load taxonomy
	as.matrix(read.table(file=tax, sep="\t"))->taxonomy
	
#remove contaminants
	otureads=remove_contaminant_taxa(otureads, taxonomy, keywords=c("Chloroplast", "Mitochondria"))   # CUSTOM FUNCTION
			
#if doing concatenated R1 and R2, first remove the "concat" from the sample names in the OTU table
	colnames(otureads)=gsub("concat", "", colnames(otureads))
			
#keep only samples intended to be included
	otureads=otureads[,match(colnames(otureads), metadata[,"sample"], nomatch=0)>0]

#keep samples with more than X total reads
	otureads=otureads[,which(colSums(otureads)>=100)]  # This is a comment
	
#find low abundance OTUs to later remove from whole dataset.
	#minimum_readcount=25
	#minimum_samples=5
	#quantitative_OTUs=rownames(otureads)[which(apply(otureads, 1, max)>=minimum_readcount)]
	#quantitative_OTUs=c("HOST", quantitative_OTUs)

#remove low abundance rows 0.5% used 20200229 for DISPLAY ONLY. 0.05 good for 16S comparison
	#otureads=low_abundance(otureads, percent=0.05)
	
#order full table by rowsums
	order(rowSums(otureads), decreasing=FALSE)->ordering_vector
		otureads[ordering_vector,]->otureads

#remove top Pseudo and Sphingo otus and add them back after family classification
CONCAT=TRUE
if(CONCAT==FALSE){
	HOST=otureads[which(rownames(otureads)=="Otu3"),]
	S63D7=otureads[which(rownames(otureads)=="Otu5"),]
	DC3000=otureads[which(rownames(otureads)=="Otu2"),]
	FR1=otureads[which(rownames(otureads)=="Otu4"),]
	p25c2=otureads[which(rownames(otureads)=="Otu1"),]
	otureads2=otureads[match(rownames(otureads), c("Otu1", "Otu2", "Otu3", "Otu4", "Otu5"), nomatch=0)==0,]
}
if(CONCAT==TRUE){
	HOST=otureads[which(rownames(otureads)=="Otu4"),]
	S63D7=10*otureads[which(rownames(otureads)=="Otu6"),]
	DC3000=otureads[which(rownames(otureads)=="Otu2"),]
	FR1=10*otureads[which(rownames(otureads)=="Otu5"),]
	p25c2=otureads[which(rownames(otureads)=="Otu1"),]
	otureads2=otureads[match(rownames(otureads), c("Otu1", "Otu2", "Otu4", "Otu5", "Otu6"), nomatch=0)==0,]
}

	
#include only phylum-classified to remove mitochondria
	otureads2<-filter_by_taxonomyic_level(countTable=otureads2, taxonomy=taxonomy, keywords=c("p"), last_filter=FALSE)

#keep only those classified to this taxonomy
	#otureads2<-filter_by_taxonomyic_level(countTable=otureads2, taxonomy=taxonomy, keywords=c("f"), last_filter=TRUE)
	otureads2<-filter_by_taxonomyic_level(countTable=otureads2, taxonomy=taxonomy, keywords=c("g"), last_filter=TRUE)

#convert all non-added Sphingomonadaceae to "otherSphingo"
#Sphingomonadaceae=otureads2["Sphingomonadaceae",]
#Pseudomonadaceae=otureads2["Pseudomonadaceae",]
Sphingomonas=otureads2["Sphingomonas",]
Pseudomonas=otureads2["Pseudomonas",]

#pull out those families
#otureads2=otureads2[match(rownames(otureads2), c("Sphingomonadaceae", "Pseudomonadaceae"), nomatch=0)==0,]
otureads2=otureads2[match(rownames(otureads2), c("Sphingomonas", "Pseudomonas"), nomatch=0)==0,]
	
#convert all non-added organisms to "remainder"
remainder=colSums(otureads2)


#add back known OTUs

	#otureads2=rbind(otureads2, HOST, S63D7, DC3000, FR1, p25c2)
	#otureads2=rbind(remainder, HOST, S63D7, DC3000, FR1, p25c2, Sphingomonadaceae, Pseudomonadaceae)
	otureads2=rbind(remainder, HOST, S63D7, DC3000, FR1, p25c2, Sphingomonas, Pseudomonas)
	
	
	order(apply(otureads2, 1, median), decreasing=FALSE)->ordering_vector
	otureads2[ordering_vector,]->otureads2

	
	toporder=c("empty", "low_abundance")
	bottomorder=c("Sphingomonadaceae", "SphASV1", "Pseudomonadaceae", "PseASV1")	
	bottomorder=c("DC3000", "p25c2", "Pseudomonas", "S63D7", "FR1", "Sphingomonas", "remainder")
		otureads2=topOrder(otureads2, toporder, bottomorder)

	abundant_families=names(which(apply(otureads2, 1, median)>=1))


#calculate samples with high noise because of low host amplicon.
#Do this before normalizing
	T=colSums(otureads2)
#acceptable percentages
	N=0.22
	minimum_percent_host=100*sqrt(-N*((-2*T)-N))/(N*T)

#normalize to %
	otureads2=normalize100(otureads2)

colnames(otureads2[,which(otureads2["HOST",]-minimum_percent_host < 0)])

actual_minus_acceptable=otureads2["HOST",]-minimum_percent_host

exclude_inaccurate_load=names(actual_minus_acceptable[which(actual_minus_acceptable<0)])

keep_accurate_load=names(actual_minus_acceptable[which(actual_minus_acceptable>=0)])


#convert to LOAD 
divide_by_host=function(tab){
	tab["HOST", which(tab["HOST",]==0)]<-0.00000001
	for(c in 1:ncol(tab)){
		tab[,c]/tab["HOST",c]->tab[,c]
	}
	return(tab)
}
otureads2_load=divide_by_host(otureads2)
otureads2_load=otureads2_load[order(rowSums(otureads2_load)),]

#remove host amplicon
otureads2_load=otureads2_load[which(rownames(otureads2_load)!="HOST"),]



#find maximum reliable load in this dataset
max_believeble_load=max(colSums(otureads2_load)[keep_accurate_load])

#adjust_unbelievable_samples_to_max
adjust_to_max=function(tab){
	for(c in 1:ncol(tab)){
		if(sum(tab[,c])> max_believeble_load){
			tab[,c]=tab[,c] * ( max_believeble_load / sum(tab[,c]) ) 
		}
	}
	return(tab)
}

otureads2_load=adjust_to_max(otureads2_load)

#otureads2_load=otureads2_load[,keep_accurate_load]



melt(otureads2_load)->histDL
c("organism", "sample_name", "abundance")->names(histDL)

histDL$genotype=factor(metadata[match(histDL$sample_name, metadata[,"sample"]),"genotype"], levels=c("WT", "eds", "coi"))
histDL$sphingo=factor(metadata[match(histDL$sample_name, metadata[,"sample"]),"sphingo"], levels=c("M", "B", "F", "S"))
histDL$pseudo=factor(metadata[match(histDL$sample_name, metadata[,"sample"]),"pseudo"], levels=c("DC", "p25c2", "MgCl"))
histDL$display_order=factor(metadata[match(histDL$sample_name, metadata[,"sample"]),"display_order"], levels=c(1:191))
histDL$combined_name=factor(metadata[match(histDL$sample_name, metadata[,"sample"]),"combined_name"])

histDL$pointcolor=rep("black", nrow(histDL))
histDL$pointcolor[match(histDL$sample_name, keep_accurate_load, nomatch=0)==0]<-"red"


histDL$organism=factor(histDL$organism, levels=unique(as.matrix(histDL$organism)))

#order by display group
#histDL$display_group=factor(metadata[match(histDL$sample_name, metadata[,"sequencefile"]),"display_group"], levels=1:18)
	#histDL=histDL[order(as.numeric(as.matrix(histDL$display_group))),]
#histDL$display_order=factor(metadata[match(histDL$sample_name, metadata[,"sequencefile"]),"display_order"], levels=1:125)
	
histDL=histDL[order(as.numeric(as.matrix(histDL$display_order))),]

histDLexclude=which(match(histDL$sample_name, exclude_inaccurate_load, nomatch="0")>0)
histDL$combined_name=as.matrix(histDL$combined_name)
histDL$combined_name[histDLexclude]<-paste("EXCLUDEEXCLUDE_",histDL$combined_name[histDLexclude], sep="")
histDL$combined_name=factor(histDL$combined_name, levels=unique(as.matrix(histDL$combined_name)))


#make more informative name column
histDL$gen_dpi_lys=apply(cbind(as.matrix(histDL$genotype), as.matrix(histDL$sphingo), as.matrix(histDL$pseudo)), 1, paste, collapse="_")

histDL$gen_dpi_lys=factor(histDL$gen_dpi_lys, levels=unique(as.matrix(histDL$gen_dpi_lys)))
histDL$sample_name=factor(histDL$sample_name, levels=unique(as.matrix(histDL$sample_name)))


levels(histDL$organism)->taxa_color_pairs
	COLUMN_IN_COLORSLIST=3    # 3 is colors, 4 is black / white / green / purple
	#ALPHA_COLUMN_IN_COLORSLIST=7
	source("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Regalado2017/family_colorscheme.R")
	"#C577B1"->taxa_color_pairs[which(taxa_color_pairs[,1]=="DC3000"),2]
	"#85D4EC"->taxa_color_pairs[which(taxa_color_pairs[,1]=="p25c2"),2]
	"black"->taxa_color_pairs[which(taxa_color_pairs[,1]=="low_abundance"),2]
	"grey"->taxa_color_pairs[which(taxa_color_pairs[,1]=="remainder"),2]
	"green"->taxa_color_pairs[which(taxa_color_pairs[,1]=="FR1"),2]
	"orange"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Sphingomonas"),2]
	"blue"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Pseudomonas"),2]
	"magenta"->taxa_color_pairs[which(taxa_color_pairs[,1]=="S63D7"),2]
	"#c7f214"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Sphingobacteriaceae"),2]
	taxa_color_pairs[match(levels(histDL$organism), taxa_color_pairs[,1]),2]->DerekColors
	histDL$alphas=rep(1, nrow(histDL))
	histDL$alphas[match(histDL$organism, abundant_families, nomatch=0)>0]=1

	
histDL3=histDL

#histDL3$organism=factor(histDL3$organism, levels=c("Pseudomonadaceae", "Sphingomonadaceae", "FR1", "S63D7", "remainder", "p25c2", "DC3000")) 
histDL3$organism=factor(histDL3$organism, levels=c("remainder", "Pseudomonas", "Sphingomonas", "FR1", "S63D7", "p25c2", "DC3000")) 
taxa_color_pairs[match(levels(histDL3$organism), taxa_color_pairs[,1]),2]->DerekColors

#histDL3$organism=factor(histDL3$organism, levels=c("Sphingomonadaceae", "Pseudomonadaceae", "FR1", "S63D7", "remainder", "DC3000", "p25c2"))


#pdf(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/soil_protection_photos_202107/16S_hamPCR/Test_Sphingo_soil_", date, ".pdf", sep="", collapse=""), width = 10, height = 4, useDingbats=FALSE)
ggplot(histDL3,aes(x=combined_name, y=abundance, fill=organism)) + theme(legend.position="none") +    #took out alpha=alphas
 	 geom_bar(position="stack", stat="identity", width=1.0) + scale_fill_manual(values=DerekColors) +
 	  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
 	 xlab("") + ylab("") + guides(fill=FALSE) +
 	 theme_classic() +
 	 theme(axis.text.x = element_text(size=2, angle = 90, vjust = 0.5, hjust=1)) 
		scale_alpha_continuous(range = c(0.3, 1)) 
		
		#scale_y_continuous(name = "abundance", breaks = c(0, .5, 1), labels = c(0, 50, 100))
#dev.off()








	
#GET CFU INFO AND ADD IT

CFU=data.frame(
	organism=factor(rep(1, 48), levels=c(1)),
	sample_name=factor(49:96, levels=49:96),
	abundance=as.numeric(as.matrix(metadata[match(49:96, metadata[,"display_order"]),"PSE_CFU_perG"]))
)

#pdf(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/pseudo_vs_sphingo_lysates/Pseudo_Sphing_compHIST_CFU2_", date, ".pdf", sep="", collapse=""), width = 4, height = 4, useDingbats=FALSE)
ggplot(CFU,aes(x=sample_name, y=abundance, fill=organism)) + theme(legend.position="none") + 
 	 geom_bar(position="stack", stat="identity", width=1.0) + scale_fill_manual(values="black") +
 	  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
 	 xlab("") + ylab("") + guides(fill=FALSE) +
 	 theme_classic() +
 	 theme(axis.text.x = element_text(size=2, angle = 90, vjust = 0.5, hjust=1)) 
#dev.off()




###MAKE BOXPLOT DC3000 vs OTU5 vs MgCl on 

#make more informative name column

histDL$categories=apply(cbind(as.matrix(histDL$genotype), as.matrix(histDL$sphingo), as.matrix(histDL$pseudo)), 1, paste, collapse="_")


### THIS PLOT SHOWS PSEUDOMONAS ON THE DIFFERENT PLANT GENOTYPES
### SHOWS COI1 IS RESISTANT TO DC3000
histDL4=histDL
#histDL4=histDL4[c(which(histDL4$organism=="DC3000"), which(histDL4$organism=="p25c2")),]
histDL_DC=histDL4[intersect(which(histDL4$pseudo=="DC"), which(histDL4$organism=="DC3000")),]
histDL_p25=histDL4[intersect(which(histDL4$pseudo=="p25c2"), which(histDL4$organism=="p25c2")),]
histDL4=rbind(histDL_DC, histDL_p25)
histDL4$organism=factor(histDL4$organism, levels=unique(as.matrix(histDL4$organism)))

#pdf(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/pseudo_vs_sphingo_lysates/Pseudo_Sphing_compBOX_Athaliana", date, ".pdf", sep="", collapse=""), width = 10, height = 6, useDingbats=FALSE)

ggplot(histDL4,aes(x=genotype,y=abundance, fill=organism)) + 
  	 	geom_boxplot(aes(color=organism), outlier.size=0, outlier.shape=NA, width=.7, position=position_dodge(.7)) + 
   	 	geom_point(aes(fill=organism, alpha=0.4), size = 3, stroke=0.1, shape = 21, position = position_jitterdodge(jitter.width=0.4)) +
   	 	theme_classic() + scale_color_manual(values=rep(c("black", "black"), 3)) + scale_fill_manual(values=rep(c("lightblue", "pink"), 3)) + 
    	theme(axis.text=element_text(size=3, color="black"), axis.title=element_text(size=14,face="bold")) + 
    	theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, color="black")) + 
    	theme(axis.ticks = element_line(colour = 'black')) +
    	scale_x_discrete(name = "injection") #+ 
 		#scale_y_continuous(limits = c(0,100)) 

#dev.off()


#MAKE VERSION WHERE LESS BELIEVABLE ARE DIFFERENT COLOR

histDL4$exclude=rep(1, nrow(histDL4))
histDL4$exclude[which(match(histDL4$sample_name, exclude_inaccurate_load, nomatch=0)>0)]<-24
histDL4$exclude[which(match(histDL4$sample_name, keep_accurate_load, nomatch=0)>0)]<-21
#histDL4$exclude=factor(histDL4$exclude, levels=c("K", "E"))

pdf(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/pseudo_vs_sphingo_lysates/Pseudo_ColEdsCoi_Athaliana", date, ".pdf", sep="", collapse=""), width = 5, height = 6, useDingbats=FALSE)
ggplot(histDL4,aes(x=genotype,y=abundance, fill=organism)) + 
  	 	geom_boxplot(aes(color=organism), outlier.size=0, outlier.shape=NA, width=.7, position=position_dodge(.7)) + 
   	 	geom_point(aes(fill=organism), alpha=.8, size = 3, stroke=0.2, shape = histDL4$exclude, position = position_jitterdodge(jitter.width=0.4)) +
   	 	theme_classic() + scale_color_manual(values=rep(c("black", "black"), 3)) + scale_fill_manual(values=rep(c("cyan", "magenta"), 3)) + 
    	theme(axis.text=element_text(size=3, color="black"), axis.title=element_text(size=14,face="bold")) + 
    	theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, color="black")) + 
    	theme(axis.ticks = element_line(colour = 'black')) +
    	scale_x_discrete(name = "injection") #+ 
 		#scale_y_continuous(limits = c(0,100)) 
dev.off()


### THIS PLOT SHOWS SPHINGOMONAS LOAD ON THE DIFFERENT PLANT GENOTYPES
histDL4=histDL

#select only plants sprayed with MgCl
histDL4=histDL4[which(histDL4$pseudo=="MgCl"),]
#histDL4=histDL4[which(histDL4$genotype=="C"),]


histDL_FR1=histDL4[intersect(which(histDL4$sphingo=="F"), 
	c(which(histDL4$organism=="FR1"),which(histDL4$organism=="S63D7"),which(histDL4$organism=="Sphingomonadaceae"))),]
histDL_S63D7=histDL4[intersect(which(histDL4$sphingo=="S"), 
	c(which(histDL4$organism=="FR1"),which(histDL4$organism=="S63D7"),which(histDL4$organism=="Sphingomonadaceae"))),]
histDL4=rbind(histDL_FR1, histDL_S63D7)
histDL4$organism=factor(histDL4$organism, levels=unique(as.matrix(histDL4$organism)))

#pdf(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/pseudo_vs_sphingo_lysates/Pseudo_Sphing_compBOX_Athaliana", date, ".pdf", sep="", collapse=""), width = 10, height = 6, useDingbats=FALSE)

ggplot(histDL4,aes(x=genotype,y=sqrt(sqrt(abundance)), fill=organism)) + 
  	 	geom_boxplot(aes(color=organism), outlier.size=0, outlier.shape=NA, width=.7, position=position_dodge(.7)) + 
   	 	geom_point(aes(fill=organism, alpha=0.4), size = 3, stroke=0.1, shape = 21, position = position_jitterdodge(jitter.width=0.4)) +
   	 	theme_classic() + scale_color_manual(values=rep(c("black", "black", "black"), 3)) + scale_fill_manual(values=rep(c("#C576B0", "#F5EB28", "green"), 3)) + 
    	theme(axis.text=element_text(size=3, color="black"), axis.title=element_text(size=14,face="bold")) + 
    	theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, color="black")) + 
    	theme(axis.ticks = element_line(colour = 'black')) +
    	scale_x_discrete(name = "injection") #+ 
 		#scale_y_continuous(limits = c(0,100)) 

#dev.off()

#DC3000 = C576B0
#p25c2 = F5EB28

histDL2=histDL

#focus="DC3000"
#focus="p25c2"

focus="hi"


histDL2=histDL2[grep("DC|p25c2", histDL2$pseudo),]

if(focus=="DC3000"){ histDL2=histDL2[grep("DC|p25c2", histDL2$pseudo),] }
if(focus=="p25c2"){ histDL2=histDL2[grep("p25c2|DC3000", histDL2$pseudo),] }
if(focus=="MgCl"){ histDL2=histDL2[grep("none|MgCl", histDL2$pseudo),] }

unique(paste(histDL2$genotype, histDL2$sphingo, histDL2$pseudo, sep="_"))

if(focus=="DC3000"){ histDL2=histDL2[grep("DC3000|Sph131_Sph18|S63D7|FR1|remainder", histDL2$organism, invert=FALSE),] 
	histDL2$organism=factor(histDL2$organism, levels=c("FR1", "S63D7", "remainder", "DC3000")) 
	histDL2$sample_category<-paste(histDL2$genotype, histDL2$sphingo, histDL2$pseudo, sep="_")
	histDL2$sample_category=factor(histDL2$sample_category, levels=c(
	"WT_F_DC",  "WT_S_DC", "WT_B_DC", "WT_M_DC",
	"eds_F_DC",  "eds_S_DC", "eds_B_DC", "eds_M_DC",
	"coi_F_DC",  "coi_S_DC", "coi_B_DC", "coi_M_DC"))
}
if(focus=="p25c2"){ histDL2=histDL2[grep("p25c2|Sph131_Sph18|S63D7|FR1|remainder", histDL2$organism, invert=FALSE),] 
	histDL2$organism=factor(histDL2$organism, levels=c("FR1", "S63D7", "remainder", "p25c2"))
	histDL2$sample_category<-paste(histDL2$genotype, histDL2$sphingo, histDL2$pseudo, sep="_")
	histDL2$sample_category=factor(histDL2$sample_category, levels=c(
	"WT_F_p25c2",  "WT_S_p25c2", "WT_B_p25c2", "WT_M_p25c2",
	"eds_F_p25c2",  "eds_S_p25c2", "eds_B_p25c2", "eds_M_p25c2",
	"coi_F_p25c2",  "coi_S_p25c2", "coi_B_p25c2", "coi_M_p25c2")) 
}
histDL2=histDL2[grep("DC3000|p25c2|S63D7|FR1|Sphingomonadaceae|Pseudomonadaceae|remainder", histDL2$organism, invert=FALSE),] 
histDL2$organism=factor(histDL2$organism, levels=c("FR1", "S63D7", "Sphingomonadaceae","Pseudomonadaceae", "remainder", "DC3000", "p25c2")) 
histDL2$sample_category<-paste(histDL2$genotype, histDL2$sphingo, histDL2$pseudo, sep="_")
histDL2$sample_category=factor(histDL2$sample_category, levels=c(
"WT_F_DC",  "WT_S_DC", "WT_B_DC", "WT_M_DC",
"eds_F_DC",  "eds_S_DC", "eds_B_DC", "eds_M_DC",
"coi_F_DC",  "coi_S_DC", "coi_B_DC", "coi_M_DC",
"WT_F_p25c2",  "WT_S_p25c2", "WT_B_p25c2", "WT_M_p25c2",
"eds_F_p25c2",  "eds_S_p25c2", "eds_B_p25c2", "eds_M_p25c2",
"coi_F_p25c2",  "coi_S_p25c2", "coi_B_p25c2", "coi_M_p25c2")) 




#if(focus=="OTU5"){ histDL2=histDL2[grep("OTU5|Sph131_Sph18|S63D7|FR1|remainder", histDL2$organism, invert=FALSE),] }
#if(focus=="BOIL"){ histDL2=histDL2[grep("OTU5|DC3000|Sph131_Sph18|S63D7|FR1|remainder", histDL2$organism, invert=FALSE),] }
histDL2$organism=factor(histDL2$organism, levels=unique(as.character(histDL2$organism)))
	
#log10(histDL2$abundance+.00000000000001)->histDL2$abundance
#histDL2=histDL2[histDL2$abundance>-5,]

	as.matrix(unique(histDL2$organism))->taxa_color_pairs
	COLUMN_IN_COLORSLIST=3    # 3 is colors, 4 is black / white / green / purple
	#ALPHA_COLUMN_IN_COLORSLIST=7
	source("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Regalado2017/family_colorscheme.R")
	"green"->taxa_color_pairs[which(taxa_color_pairs[,1]=="HOST"),2]
	"gray"->taxa_color_pairs[which(taxa_color_pairs[,1]=="low_abundance"),2]
	"#FFBC00"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Caulobacteraceae"),2]
	"#f43df0"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Flavobacteriaceae"),2]
	"#c7f214"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Sphingobacteriaceae"),2]
	"lightblue"->taxa_color_pairs[which(taxa_color_pairs[,1]=="p25c2"),2]
	"blue"->taxa_color_pairs[which(taxa_color_pairs[,1]=="DC3000"),2]
	"purple"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Sph131_Sph18"),2]
	"orange"->taxa_color_pairs[which(taxa_color_pairs[,1]=="S63D7"),2]
	"green"->taxa_color_pairs[which(taxa_color_pairs[,1]=="FR1"),2]
	"yellow"->taxa_color_pairs[which(taxa_color_pairs[,1]=="remainder"),2]


	taxa_color_pairs[match(levels(histDL2$organism), taxa_color_pairs[,1]),2]->histDL2$colors
	
	histDL2$alphas=rep(0.3, nrow(histDL2))



date=format(Sys.Date(), format="%Y%m%d")
#pdf(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/soil_protection_photos_20201007/Soil_DC3000_growth_", date, ".pdf", sep="", collapse=""), width = 7, height = 3.5, useDingbats=FALSE)
ggplot(histDL2, aes(fill=organism, y=sqrt(abundance), x=sample_category)) + 
   	 	geom_boxplot(outlier.size=0, outlier.shape=NA, alpha=0.6, width=.7, lwd=.1, position=position_dodge(.7)) + 
   	 	geom_point(aes(fill=organism), size = 1, stroke=0.1, shape = 21, position = position_jitterdodge(jitter.width=0.1, dodge=.7)) +
   	 	theme_classic() + scale_color_manual(values=histDL2$colors[1:5]) + scale_fill_manual(values=histDL2$colors[1:8]) + 
    	theme(axis.text=element_text(size=5, color="black"), axis.title=element_text(size=14,face="bold")) + 
    	theme(axis.text.x = element_text(angle=90, hjust=0.99, vjust=0.9, color="black")) + 
    	theme(axis.ticks = element_line(colour = 'black')) 
 		#scale_y_continuous(name = "log10 abundance", breaks=c(-3:3), labels=10^c(-3:3), limits = c(-3,3)) 
#dev.off()



 	 	
 	



