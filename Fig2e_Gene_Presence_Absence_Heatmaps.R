
library("ape")
library("phangorn")
library("phytools")
library("geiger")
library("ape")
library("Biostrings")
library("ggplot2")

genus="Sphingomonas"
date=format(Sys.Date(), format="%Y%m%d")
source("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/scripts/R/functions/microbiome_custom_functions.R")
	
#Pseudomonas
if(genus=="Pseudomonas"){
	PA=read.table(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Pseudomonas/GenePresenceAbsence/genePresence.aln"), sep="\t")
	core=read.tree(file = "/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Pseudo_runs/Pseudo_BigRun_20190510/strain_tree.nwk") 
	gsub("_.*", "", core$tip.label)->core$tip.label

}
if(genus=="Sphingomonas"){

	#PA=read.table(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_buscoGT85/genePresence.aln"), sep="\t")
	#PAorder=read.table(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_buscoGT85/PA_cluster_order.txt"), sep="\t")
	#PAkey=as.matrix(read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_buscoGT85/cluster_to_gene_table.txt", sep="\t", quote = ""))
	#PA_interesting_genelist=as.matrix(read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/PA_interesting_genelist.txt", sep="\t", quote = ""))
	
	basefolder="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS/"
	core=read.tree(file = paste(basefolder, "strain_tree.nwk", sep="", collapse=""))
	PA=read.table(paste(basefolder, "genePresence.aln", sep="", collapse=""),sep="\t")
	
	if("burrito"=="Xburrito"){
		install.packages("ape")
		library("ape")
		PA=read.table("/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/PanX/pan-genome-analysis/data/Sph_busco85_wONT_lessREFSEQ_noIBVSS/geneCluster/genePresence.aln",sep="\t")
		core=read.tree("/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/PanX/pan-genome-analysis/data/Sph_busco85_wONT_lessREFSEQ_noIBVSS/vis/strain_tree.nwk")
	}
}

subset=FALSE

#source("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/scripts/R/functions/microbiome_custom_functions.R")

as.matrix(PA)->PA
split_FASTA(PA, add_integer_column=FALSE)->PA
PA[,1]->samplenames
do.call(rbind,strsplit(PA[,2], split=""))->PA
apply(PA,2, as.numeric)->PA


#drop those also sequenced via nanopore from the tree which will later be used to order
Drop=which(match(core$tip.label,c("S337H113","S132H113","S136H113","S380H113","S237H113","S230H113","S190H113","S213H113","S18H113","S127H113","S133H113","S216H113"), nomatch=0)>0)
drop.tip(core, as.numeric(Drop))->core

gsub(" <unknown description>", "", samplenames)->samplenames
gsub(">", "", samplenames)->samplenames
gsub("GCF_", "GCF", samplenames)->samplenames
gsub("_.*", "", samplenames)->samplenames
gsub("GCF", "GCF_", samplenames)->samplenames


1:ncol(PA)->PAorder
as.matrix(PAorder)->PAorder

#get rid of genes present in 70% of strains (core genome)or those only absent in X strains
PAorder=PAorder[which(colSums(PA)<=(nrow(PA)-3))]  
PA=PA[,which(colSums(PA)<=(nrow(PA)-3))]  

#only consider genes found in at least X strains
PAorder=PAorder[which(colSums(PA)>=3)]
PA=PA[,which(colSums(PA)>=3)]


#library(gplots)
#library(heatmap3)
t(PA)->PA

PA->PAsaved
#PAsaved->PA


subset=FALSE
if(subset==TRUE){
	#here process the interesting gene list.
	PA_interesting_genelist
	#c("acsF","bchZ","fliD","parD","virB3","bchC","ccdb","fliE","prsD","virB4","bchE","epsF","fliP","prsE","virB8","bchL","epsM","gspD","relE","bchX","fitB","gspF","vgrG","bchY","fliC","gspK","virB1")
	genes=PA_interesting_genelist[,2]
	gene_categories=PA_interesting_genelist[,1]
	geneClusterIDs=vector(length=0)
	geneClusterIDs_detailed_categories=vector(length=0)
	geneClusterIDs_basic_categories=vector(length=0)
	for(p in genes){
		geneClusterIDs=c(geneClusterIDs, PAkey[grep(p, PAkey[,2]),1])
		geneClusterIDs_detailed_categories=c(geneClusterIDs_detailed_categories, grep(p, PAkey[,2], value=TRUE))
		geneClusterIDs_basic_categories=c(geneClusterIDs_basic_categories, rep(gene_categories[which(p==genes)], length=length(grep(p, PAkey[,2]))))
	}
	#the ones that match
	matchingClusterIDs=which(match(geneClusterIDs,PAorder, nomatch=0)>0)
	geneClusterIDs_detailed_categories[matchingClusterIDs]->geneClusterIDs_detailed_categories
	geneClusterIDs_basic_categories[matchingClusterIDs]->geneClusterIDs_basic_categories
	geneClusterIDs[matchingClusterIDs]->geneClusterIDs
	PA[match(geneClusterIDs,PAorder, nomatch=0),]->PA
	
	t(PA)->PA  #now samples are rows
	geneClusterIDs_basic_categories->colnames(PA)

	#put PA in same sample order as the core genome tree
	PA[match(core$tip.label, samplenames),]->PA
	
	#assign unique color to each gene group
	gene_group_colors=vector(length=length(geneClusterIDs_basic_categories))
	DerekPalette<-c("#000000","#ffff00","#f922b9","#32fff4","#1F78B4","#FB9A99","#FDBF6F","#E31A1C","#6A3D9A","#CAB2D6","#FF7F00","#FFFF99","#A6CEE3","#000000","#ffff00","#f922b9","#32fff4","#1F78B4","#FB9A99","#FDBF6F","#E31A1C","#6A3D9A","#CAB2D6","#FF7F00","#FFFF99","#A6CEE3")
	c(DerekPalette, DerekPalette)->DerekPalette
	for(g in 1:length(unique(geneClusterIDs_basic_categories))){
		gene_group_colors[which(geneClusterIDs_basic_categories==unique(geneClusterIDs_basic_categories)[g])]<-DerekPalette[g]
	}
	
	#PUT BLANKS BETWEEN COLUMNS CATEGORIES
	blankwidth=5
	geneClusterIDs_basic_categories->gbc
	unique(gbc)->ucols
	for(co in 1:length(ucols)){
		print(co)
		max(which(gbc==ucols[co]))->startcol
		print(startcol)
		if(co!=length(ucols)){
			PA=cbind(
				PA[,1:startcol],
				matrix(data=-1, nrow=nrow(PA), ncol=blankwidth),
				PA[,(startcol+1):ncol(PA)]
			)
			gene_group_colors=c(gene_group_colors[1:startcol], rep("white", blankwidth), gene_group_colors[(startcol+1):length(gene_group_colors)])	
			gbc=c(gbc[1:startcol], rep("white", blankwidth), gbc[(startcol+1):length(gbc)])	
		}
	}
	

	
}


#c("acsF","bchZ","fliD","parD","virB3","bchC","ccdb","fliE","prsD","virB4","bchE","epsF","fliP","prsE","virB8","bchL","epsM","gspD","relE","bchX","fitB","gspF","vgrG","bchY","fliC","gspK","virB1")




#if(subset!=TRUE){
	if(nrow(PA)>ncol(PA)){t(PA)->PA}
	if(colSums(PA)[1]<colSums(PA)[ncol(PA)]){PA[,ncol(PA):1]->PA}
	if("clusterCols"=="clusterCols"){
		t(PA)->PA
		dist(PA)->PAdist
		hclust(PAdist)->hclustPAdist
		t(PA[hclustPAdist$order,])->PA
	}
	
	if("saved"=="saved"){
		library("ape")
		core=read.tree("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS/strain_tree.nwk")
		PA=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS/PA.txt",sep="\t")
		samplenames=rownames(PA)=PA[,1]
			gsub("Sph0349","Sph0349_S337H113", samplenames)->samplenames
			gsub("Sph0136","Sph0136_S132H113", samplenames)->samplenames
			gsub("Sph0140","Sph0140_S136H113", samplenames)->samplenames
			gsub("Sph0184","Sph0184_S380H113", samplenames)->samplenames
			gsub("Sph0247","Sph0247_S237H113", samplenames)->samplenames
			gsub("Sph0240","Sph0240_S230H113", samplenames)->samplenames
			gsub("Sph0219","Sph0219_S190H113", samplenames)->samplenames
			gsub("Sph0223","Sph0223_S213H113", samplenames)->samplenames
			gsub("Sph0018","Sph0018_S18H113", samplenames)->samplenames
			gsub("Sph0131","Sph0131_S127H113", samplenames)->samplenames
			gsub("Sph0137","Sph0137_S133H113", samplenames)->samplenames
			gsub("Sph0226","Sph0226_S216H113", samplenames)->samplenames
			samplenames->rownames(PA)	
		PA[,2:ncol(PA)]->PA
		as.matrix(PA)->PA
	}
	
	nanopore_control=FALSE
	#kill the HiSeq version of the nanopore genomes from the PA matrix
	if(nanopore_control==FALSE){   
		PA=PA[which(match(samplenames,c("S337H113","S132H113","S136H113","S380H113","S237H113","S230H113","S190H113","S213H113","S18H113","S127H113","S133H113","S216H113"), nomatch=0)==0),]
		samplenames=samplenames[which(match(samplenames,c("S337H113","S132H113","S136H113","S380H113","S237H113","S230H113","S190H113","S213H113","S18H113","S127H113","S133H113","S216H113"), nomatch=0)==0)]
		#gsub("_S.*", "", core$tip.label)->tmp.tip.label   #changed this in 2022 to the below
		gsub("_nanopore.*", "", core$tip.label)->tmp.tip.label 
		PA[match(tmp.tip.label, samplenames, nomatch=0),]->PA
	}
	

		
	
	

	present=col2rgb("#000000", alpha = FALSE)/255  #black
	absent=col2rgb("#FFD6A9", alpha = FALSE)/255 #light yellow   #FFFEF4"
	
	
	#PA[,which(colSums(PA)>=10)]->PA
		#PA[,which(colSums(PA)<=(nrow(PA)-10))]->PA
	#PA[,which(colSums(PA)<=(round(nrow(PA)*.7)))]->PA
	
	
	if(nanopore_control==TRUE){
			toMatch=c("S337H113","S132H113","S136H113",
										"S380H113","S237H113","S230H113",
										"S190H113","S213H113","S18H113",
										"S127H113","S133H113","S216H113",
										".*S337H113.*",".*S132H113.*",".*S136H113.*",
										".*S380H113.*",".*S237H113.*",".*S230H113.*",
										".*S190H113.*",".*S213H113.*",".*S18H113.*",
										".*S127H113.*",".*S133H113.*",".*S216H113.*")
			Drop=grep(paste(toMatch, collapse="|"), core$tip.label, invert=TRUE)
			drop.tip(core, as.numeric(Drop))->core
			gsub("_nanopore.*","", core$tip.label)->core$tip.label
			is_tip <- core$edge[,2] <= length(core$tip.label)
			ordered_tips <- core$edge[is_tip, 2]
			gsub("-", "_", as.matrix(core$tip.label[ordered_tips]))->actualorder
			correct_order=as.matrix(c("S337H113","Sph0349_S337H113","S132H113","Sph0136_S132H113","S136H113","Sph0140_S136H113","S380H113","Sph0184_S380H113","S237H113","Sph0247_S237H113","S230H113","Sph0240_S230H113","S190H113","Sph0219_S190H113","S213H113","Sph0223_S213H113","S18H113","Sph0018_S18H113","S127H113","Sph0131_S127H113","S133H113","Sph0137_S133H113","S216H113","Sph0226_S216H113"))
			if(all(actualorder==correct_order)!=TRUE){
				print("rotating_branch")
				rotate(core, node=40)->core
			}
			is_tip <- core$edge[,2] <= length(core$tip.label)
				ordered_tips <- core$edge[is_tip, 2]
				gsub("-", "_", as.matrix(core$tip.label[ordered_tips]))->actualorder
			if(all(actualorder==correct_order)==TRUE){
				print("correct")
			}
			plot(core)
			PA[match(actualorder, samplenames, nomatch=0),]->PA
	}
	
	
	rowSums(PA)
	
	S18H113 nanopore = 4093
	S18H113 HiSeq = 5143
	
	S213H113 nanopore = 3645
	S213H113 HiSeq = 3654
	
	S127H113 nanopore = 3821
	S127H113 HiSeq = 4017
	
	PA[grep("S380H113", PA),]->S380
	length(which((S380["Sph0184_S380H113",]-S380["S380H113",])==-1))  #shows non nanopore has 56 more genes
	which((S380["Sph0184_S380H113",]-S380["S380H113",])==-1)
	
	#convert to array to rasterize
	P=PA
	A=PA+1
	A[which(A==2)]=0
	#start from HERE to regenerate PA
	array(c((A*absent[1]+P*present[1]), 
		(A*absent[2]+P*present[2]), 
		(A*absent[3]+P*present[3])), c(nrow(P), ncol(P), 3))->PA2
		
	#samplenames->rownames(PA)
	#write.table(PA, "/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS/PA.txt",
	#		  quote =FALSE, sep = "\t", row.names =TRUE, col.names = FALSE)
	
	library(png)
	

	PA2[,ncol(PA2[,,1]):1,]->PA2	

	
	
	pdf(file=paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/AUGUSTUS/BLAST_pdfs/", "allPA_noIBVSS_", date, "_noncore_NANO.pdf", sep="", collapse=""), width = 5, height = 5, useDingbats=FALSE)
		plot.new() 
		rasterImage(PA2[nrow(PA2):1,,], 0,0,1,1, interpolate=TRUE) 
	dev.off()	
	
	
	#2022 make dist matrix to compare trees. #dist looks for distances between rows.
	#trying to support claim that presence absence correlates with overall genome distance.

	#now upload the MASH distance matrix of whole genome similarity
	#fullmatrix_Sph<-read.table(file="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/MASH/fullmatrix.txt", sep="\t")
	fullmatrix_Sph<-read.table(file="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/coreGenome/coreGenome_MASHmat.txt", sep="\t")
	as.matrix(fullmatrix_Sph)->fullmatrix_Sph
	fullmatrix_Sph[1,]->samplenames_Sph
	fullmatrix_Sph[2:nrow(fullmatrix_Sph),]->fullmatrix_Sph
	apply(fullmatrix_Sph, 2, as.numeric)->fullmatrix_Sph
	gsub("_pilon.fasta", "", samplenames_Sph)->samplenames_Sph
	gsub("_coreAln.fasta", "", samplenames_Sph)->samplenames_Sph
	gsub("_nanopore.*", "", samplenames_Sph)->samplenames_Sph
	gsub(".*_", "", samplenames_Sph)->samplenames_Sph
	
	samplenames_Sph->colnames(fullmatrix_Sph)
	samplenames_Sph->rownames(fullmatrix_Sph)
		
	#Next make the distance matrix of the PA data.  BLOCK THIS OUT IF YOU ALREADY DEFINED IT BECAUSE IT TAKES FOREVER
	
	library(vegan)
	testDIST=vegdist(PA, method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
	#testDIST=dist(PA, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
	
	
	#DEFINING M
	m <- data.frame(t(combn(rownames(PA),2)), as.numeric(testDIST))
	names(m) <- c("c1", "c2", "distance")
	m[,1]=gsub(".*_", "", m[,1])
	m[,2]=gsub(".*_", "", m[,2])
	m[,3]=as.numeric(as.matrix(m[,3]))
	m=m[grep("S", m[,1]),]
	m=m[grep("S", m[,2]),]
	
	library("reshape2")	
	#paste the first and second columns together to make a unique comparison string
	m$comparison_string=apply(m[,1:2], 1, paste, collapse="_", sep="_")
	
	
	#install.packages("reshape2")


	#now reduce the MASH matrix to inclu#SCP files from burrito
#scp dlundberg@burrito.eb.local:/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/process_scrape_metagenomes/*omonas_S*/mapping_table/*omonas_S*_full_mapping_table.txt /Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes

{   #TO CLEAN UP VIEWING AREA AND SHRING THIS STUFF...

library(ape)


date=format(Sys.Date(), format="%Y%m%d")
user="Derek"
genus="Pseudomonas"

#Pratchaya paths
if(user=="Pratchaya"){
	microbiome_custom_functions="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/scripts/R/functions/microbiome_custom_functions.R"
	full.table.spring="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/PseudoSpring_GENOMES_table.txt"
}



#Derek paths
if(user=="Derek"){
	microbiome_custom_functions="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/scripts/R/functions/microbiome_custom_functions.R"
	if(genus=="Pseudomonas"){
		#full.table.spring="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/PseudoSpring_GENOMES_table.txt"
		#full.table.summer="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/PseudoSummer_GENOMES_table.txt"
		full.table.spring.path="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/Pseudomonas_Spring_full_mapping_table.txt"
		full.table.summer.path="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/Pseudomonas_Summer_full_mapping_table.txt"
			basefolder="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Pseudo_runs/Pseudo_BigRun_20190510/"
			core=read.tree(file = paste(basefolder, "strain_tree.nwk", sep="", collapse=""))

		location_to_write_combined_SpringSummer_table="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/PseudoSpringSummer_COMBINED_table_20190812.txt"
	}	
	if(genus=="Sphingomonas"){
		#full.table.spring="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/SphingoSpring_GENOMES_table.txt"
		#full.table.summer="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/SphingoSummer_GENOMES_table.txt"
		full.table.spring.path="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/Sphingomonas_Spring_full_mapping_table.txt"
		full.table.summer.path="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/Sphingomonas_Summer_full_mapping_table.txt"
			basefolder="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS/"
			core=read.tree(file = paste(basefolder, "strain_tree.nwk", sep="", collapse=""))
		location_to_write_combined_SpringSummer_table="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/SphingoSpringSummer_COMBINED_table_20190812.txt"
	}	
	all_bulk_metagenomes_metadata="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/all_bulk_metagenomes_metadata.txt"
	location_to_write_finished_figures="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/"
	GCF_metadata=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/NCBI_GCF_metadata.txt", sep="\t")

}

blankcolor="black"


source(microbiome_custom_functions)
library(gplots)

#Pseudomonas
#LOAD THE MATRICES 
full.table.spring=read.table(full.table.spring.path) #, header = T, row.names = 1, sep = "\t")
full.table.summer=read.table(full.table.summer.path) #, header = T, row.names = 1, sep = "\t")

#READ METADATA FOR THE METAGENOMES
metagenomes_metadata=as.matrix(read.table(all_bulk_metagenomes_metadata, sep = "\t", comment.char = ""))
	gsub("NK_", "NeckarKing_", metagenomes_metadata)->metagenomes_metadata
#READ METADATA FOR THE MAPPING GENOMES
#genomes_metadata=read.table(all_mapping_genomes_metadata, sep = "\t", comment.char = "")

#READ IN 16S SEQUENCES FOR THE LOCAL MAPPING GENOMES
culture_fasta=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/KarasovAlmario16S/all_cultured_16S_20190410_reoriented_trimmed.fa")
as.matrix(culture_fasta)->culture_fasta
if(genus=="Pseudomonas"){
	culture_fasta[c(
						sort(c(grep("Pseudo", culture_fasta), (grep("Pseudo", culture_fasta)+1))),
				  		sort(c(grep("p.*16S", culture_fasta), (grep("p.*16S", culture_fasta)+1)))
				  	),]->Pseudo16S_culture
	split_FASTA(Pseudo16S_culture, add_integer_column=FALSE, leave_carrot=FALSE, reverse=FALSE)->Pseudo16S_table
	gsub("H133.Pseudo_16S", "H133" , Pseudo16S_table[,1])->Pseudo16S_table[,1]
	gsub("Pseudo", "p" , Pseudo16S_table[,1])->Pseudo16S_table[,1]
	gsub("_16S", "" , Pseudo16S_table[,1])->Pseudo16S_table[,1]
}
if(genus=="Sphingomonas"){
	culture_fasta[sort(c(grep("Sphingo", culture_fasta), (grep("Sphingo", culture_fasta)+1))),]->Sphingo16S_culture
	split_FASTA(Sphingo16S_culture, add_integer_column=FALSE, leave_carrot=FALSE, reverse=FALSE)->Sphingo16S_table
	gsub("H133.Sphingo_16S", "H133" , Sphingo16S_table[,1])->Sphingo16S_table[,1]
	gsub("Sphingo", "p" , Sphingo16S_table[,1])->Sphingo16S_table[,1]
	gsub("_16S", "" , Sphingo16S_table[,1])->Sphingo16S_table[,1]
}

#FORMAT COLUMN NAMES
as.character(colnames(full.table.spring))->rnames.spring
as.character(colnames(full.table.summer))->rnames.summer
if(genus=="Sphingomonas"){
	gsub("_export.txt", "H104", rnames.spring)->rnames.spring
	gsub("uniqueCounts.txt", "H104", rnames.spring)->rnames.spring	
}
if(genus=="Pseudomonas"){
	gsub("_export.txt", "H103", rnames.spring)->rnames.spring
	gsub("uniqueCounts.txt", "H103", rnames.spring)->rnames.spring
}
gsub("_export.txt", "H135", rnames.summer)->rnames.summer
gsub("uniqueCounts.txt", "H135", rnames.summer)->rnames.summer
colnames(full.table.spring)=rnames.spring
colnames(full.table.summer)=rnames.summer

#COMBINE SPRING AND SUMMER TABLES INTO ONE MATRIX
as.matrix(full.table.spring)->full.table.spring
as.matrix(full.table.summer)->full.table.summer
combine_taxa_tables(full.table.spring, full.table.summer)->combined


#RENAME COLNAMES (PLANT SAMPLES)
metagenomes_metadata[match(colnames(combined), metagenomes_metadata[,3], nomatch=0),18]


#THRESHOLDS
#FILTER SAMPLES WITHOUT ENOUGH READS
if(genus=="Pseudomonas"){
	combined[,which(colSums(combined)>=20000)]->combined
}

if(genus=="Sphingomonas"){
	#combined[,which(colSums(combined)>=20000)]->combined
	combined[,which(colSums(combined)>=15000)]->combined
}

#RENAME COLNAMES (PLANT SAMPLES)
metagenomes_metadata[match(colnames(combined), metagenomes_metadata[,3], nomatch=0),18]


#ALSO THROW THOSE SAMPLES THAT HAD LESS THAN 10 VISIBLE COLONIES PER METAGENOME
#MAY2022 I realized the existing code was matching both sphingomonas and pseudomonas CFU counts to the same sample ID.
#I created a new table with the same information in which each sample ID was unique.

CFU<-read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/Final_CFU_counts_SpringSummer.txt", sep="\t")
as.matrix(CFU[1,])->names(CFU)
gsub("NK", "NeckarKing", CFU[,1])->CFU[,1]
CFU[2:nrow(CFU),]->CFU
as.numeric(as.matrix(CFU$Pseudo2))->CFU$Pseudo2
as.numeric(as.matrix(CFU$LB))->CFU$LB
as.numeric(as.matrix(CFU$Sphingo))->CFU$Sphingo


CFU2022<-read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/Final2022_metagenome_CFU_counts_SpringSummer.txt", sep="\t")
as.matrix(CFU2022[1,])->names(CFU2022)
gsub("NK", "NeckarKing", CFU2022[,1])->CFU2022[,1]
CFU2022[2:nrow(CFU2022),]->CFU2022
as.matrix(CFU2022[,1])->CFU2022[,1]
as.numeric(as.matrix(CFU2022$Colonies2))->CFU2022$Colonies2

metagenomes_metadata[match(as.matrix(CFU2022[,1]), metagenomes_metadata[,18]),3]->CFU2022$alias

length(colnames(combined))

#remove NA
CFU2022[is.na(CFU2022$alias)==FALSE,]->CFU2022
CFU2022$alias[which(CFU2022$Colonies2>=10)]->keep
combined[,match(keep, colnames(combined), nomatch=0)]->combined2

combined2->combined


length(colnames(combined))




#RENAME COLNAMES (PLANT SAMPLES)
metagenomes_metadata[match(colnames(combined), metagenomes_metadata[,3], nomatch=0),18]


#DELETE #
#metagenomes_metadata[match(colnames(combined), metagenomes_metadata[,3]),18]
#H135s=grep("H135", colnames(combined), value=TRUE)
#H103s=grep("H103", colnames(combined), value=TRUE)
#metagenomes_metadata[match(H103s, metagenomes_metadata[,3]),18]
#metagenomes_metadata[match(H135s, metagenomes_metadata[,3]),18]
#metagenomes_metadata[,c(3, 18, 19)]

#make everything mapping a 3 column field, separated by underscore
#field 1 = local or RDP or Decoy
#field 2 = genome ID
#field 3 = contig ID

#FORMAT GENOME ROWNAMES
if(genus=="Pseudomonas"){
	gsub("_", "-", rownames(combined))->rownames(combined)  #change all underscores to hyphen
	gsub("^DL133-OTU5", "DL133.OTU5", rownames(combined))->rownames(combined) 
	gsub("^DL133-nonOTU5", "DL133.nonOTU5", rownames(combined))->rownames(combined) 
	gsub("^nonOTU5-TK165", "TK165.nonOTU5", rownames(combined))->rownames(combined)
	gsub("^OTU5-TK165", "TK165.OTU5", rownames(combined))->rownames(combined)
	gsub("(^[A-Za-z0-9.]+)-(.*)", "\\1_\\2", rownames(combined))->rownames(combined) #replace first hyphen with underscore 
	gsub("(.*)-([A-Za-z0-9.]+$)", "\\1_\\2", rownames(combined))->rownames(combined) #replace last hyphen with underscore 
	gsub("(.*)(PsyRun133)-(S[0-9]+)(.*)", "\\1H133\\3\\4", rownames(combined))->rownames(combined) #replace first hyphen with underscore 
	gsub("(H133)(S[0-9]+)", "\\2\\1", rownames(combined))->rownames(combined) #replace first hyphen with underscore 
}
if(genus=="Sphingomonas"){
	gsub("_", "-", rownames(combined))->rownames(combined)  #change all underscores to hyphen
	gsub("^DL133-OTU5", "DL133.OTU5", rownames(combined))->rownames(combined) 
	gsub("^DL133-nonOTU5", "DL133.nonOTU5", rownames(combined))->rownames(combined) 
	gsub("^nonOTU5-TK165", "TK165.nonOTU5", rownames(combined))->rownames(combined)
	gsub("^OTU5-TK165", "TK165.OTU5", rownames(combined))->rownames(combined)
	gsub("(^[A-Za-z0-9.]+)-(.*)", "\\1_\\2", rownames(combined))->rownames(combined) #replace first hyphen with underscore 
	gsub("(.*)-([A-Za-z0-9.]+$)", "\\1_\\2", rownames(combined))->rownames(combined) #replace last hyphen with underscore 
	gsub("(.*)(PsyRun133)-(S[0-9]+)(.*)", "\\1H133\\3\\4", rownames(combined))->rownames(combined) #replace first hyphen with underscore 
	gsub("(H133)(S[0-9]+)", "\\2\\1", rownames(combined))->rownames(combined) #replace first hyphen with underscore 
}




#WRITE COMBINED TABLE
write.table(combined, paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/PseudoSpringSummer_COMBINED_table", "_", date, ".txt", sep=""), sep = '\t')


as.matrix(metagenomes_metadata)->metagenomes_metadata
pal1 <- c(blankcolor, colorRampPalette(c("black","#3f324f","#6A1B9A","#E85285", "#FFECB3"))(n=300))


#N"#888888"


#simplify genome names for color bar labels
as.matrix(GCF_metadata)->GCF_metadata
GCF_metadata[match(rownames(combined), GCF_metadata[,2]),3]->GCF_genus
GCF_metadata[match(rownames(combined), GCF_metadata[,2]),4]->GCF_group
genomelabels2=genomelabels1=rownames(combined)
if(genus=="Sphingomonas"){
	genomelabels1[which(GCF_genus=="Sphingomonas")]<-"REFSEQ"
	genomelabels1[which(GCF_genus=="OTHER" | GCF_genus=="Pseudomonas")]<-"DECOY"
	genomelabels1[grep("S.*H.*", genomelabels1)]<-"LOCAL"
}
if(genus=="Pseudomonas"){
	#gsub("(p[0-9]+.[A-Z][0-9]+)", "TK.Pseudo", genomelabels1)->genomelabels1
	genomelabels1[which(GCF_genus=="Pseudomonas")]<-"REFSEQ"
	genomelabels1[which(GCF_genus=="OTHER" | GCF_genus=="Sphingomonas")]<-"DECOY"
	genomelabels1[grep("S.*H.*", genomelabels1)]<-"LOCAL"
}

#gsub(".*H133.*", "Local_Eyach", genomelabels1)->genomelabels1
#gsub(".*S[0-9]+.*", "Local_Roger", genomelabels1)->genomelabels1
#gsub("At_ROOT.*", "AtROOT", genomelabels1)->genomelabels1
#gsub("At_LSPHERE-Leaf.*", "AtLSPHERE", genomelabels1)->genomelabels1
#gsub("At_SOIL-Soil.*", "AtSOIL", genomelabels1)->genomelabels1
#gsub(".*NZ.*", "REFSEQ", genomelabels1)->genomelabels1
#gsub(".*NC.*", "REFSEQ", genomelabels1)->genomelabels1



#FIND THE OTU CLASSIFICATION FOR THE LOCAL GENOMES
	#Turn genome names into the actual V3V4 ASV sequences
	if(genus=="Pseudomonas"){
		Pseudo16S_table[match(rownames(combined), Pseudo16S_table[,1]),2]->associated_seqs_of_genomelabels
	}
	if(genus=="Sphingomonas"){
		Sphingo16S_table[match(rownames(combined), Sphingo16S_table[,1]),2]->associated_seqs_of_genomelabels
	}
	#Rename cultured seqs into ASVs
	as.matrix(table(associated_seqs_of_genomelabels))->unique_seqs
	unique_seqs[order(unique_seqs[,1], decreasing=TRUE),]->unique_seqs
	cbind(names(unique_seqs), unique_seqs)->unique_seqs
	rownames(unique_seqs)<-NULL
	paste("000", 1:nrow(unique_seqs), sep="")->ASVnums
	substr(ASVnums, nchar(ASVnums)-1, nchar(ASVnums))->ASVnums
	unique_seqs=cbind(unique_seqs, paste("ASV", ASVnums, sep=""))
	unique_seqs[match(associated_seqs_of_genomelabels, unique_seqs[,1]),3]->ASV_ID
	ASV_ID[is.na(ASV_ID)==TRUE]<-""
	#only mark ASV1 (the most abundant one)
	#for(asv in 2:50){
	#	ASV_ID[grep(paste("^ASV", asv, "$", sep="", collapse=""), ASV_ID)]<-""
	#}
	cbind(rownames(combined), genomelabels1, ASV_ID)->genome_name_replacement_table  
	#paste(genome_name_replacement_table[,3], genome_name_replacement_table[,2], sep="")->genome_name_replacement_table[,2]

#RENAME COLNAMES (PLANT SAMPLES)
metagenomes_metadata[match(colnames(combined), metagenomes_metadata[,3], nomatch=0),18]->colnames(combined)

#REMOVE ROGUE SAMPLES
if(genus=="Pseudomonas"){
	combined[,grep("SmallClover_S59.1_Pseudobulk", colnames(combined), invert=TRUE)]->combined
	combined[,grep("SpringKraut_S3.1_Pseudobulk", colnames(combined), invert=TRUE)]->combined
	combined[,grep("BigClover", colnames(combined), invert=TRUE)]->combined
}
#REMOVE ROGUE SAMPLES
if(genus=="Sphingomonas"){
	combined[,grep("BigClover", colnames(combined), invert=TRUE)]->combined
}

#NORMALIZE
#combined <- normalize100(combined)
#2022 NOT SURE THIS IS THE RIGHT PLACE TO NORMALIZE............


#sort columns by abundance
combined[,order(colSums(combined), decreasing=TRUE)]->combined

gsub("NK", "NeckarKing", colnames(combined))->colnames(combined)

#arrange plants in the following order
reorder_columns=c(
	grep("Athaliana", colnames(combined)),
	grep("Draba", colnames(combined)),
	grep("Cardamine", colnames(combined)),
	grep("Dandelion", colnames(combined)),
	grep("Thistle", colnames(combined)),	
	grep("Grass", colnames(combined)),
	grep("NeckarKing", colnames(combined)),
	grep("Clover", colnames(combined)),
	grep("Fuzzy", colnames(combined)),
	grep("Plantago", colnames(combined)),
	grep("SpringKraut", colnames(combined)),
	grep("Moss", colnames(combined)),
	grep("Soil", colnames(combined))
)	
combined[,reorder_columns]->combined




#Remove S168H133 from HEATMAP cause it's the only local genome not in the tree. 
if("removeS168H133"=="removeS168H133"){
	keep=grep("S168H133", genome_name_replacement_table[,1], invert=TRUE)
	combined[keep,]->combined
	genome_name_replacement_table[keep,]->genome_name_replacement_table
}

#Remove p8.A2 from HEATMAP cause it's probably a hybrid genome and the 16S sequence is messed up
if("remove_p8.A2"=="remove_p8.A2"){
	keep=grep("p8.A2", genome_name_replacement_table[,1], invert=TRUE)
	combined[keep,]->combined
	genome_name_replacement_table[keep,]->genome_name_replacement_table
}

#order rows first by total mapping abundance
combined=combined[order(rowSums(combined), decreasing=TRUE),]
#put genome name replacement table in same order
genome_name_replacement_table=genome_name_replacement_table[match(rownames(combined), genome_name_replacement_table[,1]),]
	
	
#Sort the local genomes by their order in the tree, and drop tips in tree not found in the condensed mapping
if(genus=="Sphingomonas"){
	match(gsub("_nanopore.*", "", gsub("Sph.*_S", "S", core$tip.label)), genome_name_replacement_table[,1], nomatch=0)->inTree
}
if(genus=="Pseudomonas"){
	match(gsub("_.*", "", core$tip.label), genome_name_replacement_table[,1], nomatch=0)->inTree
}	
	inTree[length(inTree):1]->inTree
	combined[inTree,]->combinedinTree
	#drop from tree those genomes not in the heatmap. Remove all nanopore so it's not duplicated.
	if(genus=="Sphingomonas"){
		keep_in_tree=match(rownames(combinedinTree), gsub(".*_nanopore.*", "", gsub("Sph.*_S", "S", core$tip.label)))
	}
	if(genus=="Pseudomonas"){
		keep_in_tree=match(rownames(combinedinTree), gsub("_.*", "", core$tip.label))
	}
		drop_from_tree=setdiff(1:length(core$tip.label), keep_in_tree)
		drop.tip(core, drop_from_tree)->core
	genome_name_replacement_table[inTree,]->genome_name_replacement_tableinTree
		rep("TREE", nrow(genome_name_replacement_tableinTree))->genome_name_replacement_tableinTree[,2]
setdiff(1:nrow(genome_name_replacement_table), inTree)->notinTree
	combined[notinTree,]->combinedNOTinTree
	genome_name_replacement_table[notinTree,]->genome_name_replacement_tableNOTinTree
rbind(combinedinTree, combinedNOTinTree)->combined
	rbind(genome_name_replacement_tableinTree, genome_name_replacement_tableNOTinTree)->genome_name_replacement_table 


#ignore heirarchical clustering and order heatmap rows by categories in the order below
genome_categories=genome_name_replacement_table[,2]
if(genus=="Pseudomonas"){	
	reorder_vector=c(
	sort(c(grep("ASV01TK.Pseudo", genome_categories),
		grep("ASV01Local_Eyach", genome_categories)
		), decreasing=FALSE),
	sort(c(grep("ASV[0123456789][023456789]TK.Pseudo", genome_categories),
		grep("ASV[123456789]1TK.Pseudo", genome_categories),
		grep("ASV[0123456789][023456789]Local_Eyach", genome_categories),
		grep("ASV[123456789]1Local_Eyach", genome_categories),
		grep("^Local_Eyach", genome_categories)
		), decreasing=FALSE),
	grep("TREE", genome_categories),
	grep("LOCAL", genome_categories),
	grep("REFSEQ", genome_categories),
	grep("DECOY", genome_categories),
	grep("AtLSPHERE", genome_categories),
	grep("AtROOT", genome_categories),
	grep("AtSOIL", genome_categories)
	)
}
if(genus=="Sphingomonas"){	
	reorder_vector=c(
			sort(c(grep("ASV01Local_Roger", genome_categories),
				grep("ASV01Local_Eyach", genome_categories)
				), decreasing=FALSE),
	sort(c(grep("ASV[0123456789][023456789]Local_Roger", genome_categories),
		grep("ASV[123456789]1Local_Roger", genome_categories),
		grep("ASV[0123456789][023456789]Local_Eyach", genome_categories),
		grep("ASV[123456789]1Local_Eyach", genome_categories),
		grep("^Local_Eyach", genome_categories)
		), decreasing=FALSE),
	grep("TREE", genome_categories),
	grep("LOCAL", genome_categories),
	grep("REFSEQ", genome_categories),
	grep("DECOY", genome_categories),
	grep("AtLSPHERE", genome_categories),
	grep("AtROOT", genome_categories),
	grep("AtSOIL", genome_categories)
	)
}
combined=combined[reorder_vector,]
genome_name_replacement_table=genome_name_replacement_table[reorder_vector,]


#Remove DECOY from HEATMAP?
if("removeDECOY"=="rRemoveDECOY"){
	keep=grep("DECOY", genome_name_replacement_table[,2], invert=TRUE)
	combined[keep,]->combined
	genome_name_replacement_table[keep,]->genome_name_replacement_table
}



#GENOME COLORS
if(genus=="Pseudomonas"){
	unique(genome_name_replacement_table[,2])
	#c("REFSEQ",      "AtLSPHERE",    "AtROOT",     "AtSOIL", "Local_Eyach", "TK.Pseudo")->  categories
	#c("#000000",   "#36B000",     "#CACACA",      "#996903",   "black",       "black")->names(categories)
	c("REFSEQ",       "DECOY",    "LOCAL", "TREE")->  categories
	c("#000000",     "#36B000",     "#CACACA", "#FFFFFF")->names(categories)
	#add colors to genome name replacement table
	names(categories)[match(genome_name_replacement_table[,2], categories)]->genome_name_replacement_table[,3]
	#make all TKpseudo black
	genome_name_replacement_table[grep("TK.Pseudo", genome_name_replacement_table[,2]),3]<-"#6715ED"
	#make all Local_Eyach black
	genome_name_replacement_table[grep("Local_Eyach", genome_name_replacement_table[,2]),3]<-"#6715ED"
	#make all "OTU5" lightgrey
	genome_name_replacement_table[grep("ASV01Local_Eyach", genome_name_replacement_table[,2]),3]<-"#15D0ED"
	genome_name_replacement_table[grep("ASV01TK.Pseudo", genome_name_replacement_table[,2]),3]<-"#15D0ED"

	#assign more colors based on ASV sequence
	unique_seqs=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/pseudo_16Sseqs_to_color_20190910.txt", sep="\t", comment.char="")
	culture_fasta=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/KarasovAlmario16S/all_cultured_16S_20190410_reoriented_trimmed.fa")
	as.matrix(culture_fasta)->culture_fasta
	culture_fasta[sort(c(grep(">", culture_fasta), grep(">", culture_fasta)+1))]->culture_fasta
	split_FASTA(culture_fasta)->culture_fasta
	gsub(">", "", culture_fasta)->culture_fasta
	gsub("GCF_", "GCF", culture_fasta)->culture_fasta
	gsub("_.*", "", culture_fasta)->culture_fasta
	gsub(".Sphingo", "", culture_fasta)->culture_fasta
	gsub(".Pseudo", "", culture_fasta)->culture_fasta
	gsub("GCF", "GCF_", culture_fasta)->culture_fasta
	as.matrix(unique_seqs[match(culture_fasta[,2], unique_seqs[,1]),4])->culture_fasta[,3]
	culture_fasta->genome_name_to_color_key
	genomes_locations_to_give_ASV_color=which(match(genome_name_replacement_table[,1], genome_name_to_color_key[,1], nomatch=0)!=0)
	
	colors_to_give_them=genome_name_to_color_key[match(genome_name_replacement_table[,1], genome_name_to_color_key[,1], nomatch=0),3]
	colors_to_give_them->genome_name_replacement_table[genomes_locations_to_give_ASV_color,3]
	#fix REFSEQ blacks
	#genome_name_replacement_table[grep("REFSEQ",genome_name_replacement_table[,2]),3]<-"#000000"
	#fix LOCAL greys
	genome_name_replacement_table[grep("LOCAL",genome_name_replacement_table[,2]),3]<-"#CACACA"
}
if(genus=="Sphingomonas"){
	unique(genome_name_replacement_table[,2])
	#c("REFSEQ",      "AtLSPHERE",    "AtROOT",     "AtSOIL", "Local_Eyach", "Local_Roger")->  categories
	#c("#000000",     "#36B000",     "#CACACA",     "#996903",   "black",       "black")->names(categories)
	c("REFSEQ",       "DECOY",    "LOCAL", "TREE")->  categories
	c("#000000",     "#36B000",     "#CACACA", "#FFFFFF")->names(categories)
	#add colors to genome name replacement table
	names(categories)[match(genome_name_replacement_table[,2], categories)]->genome_name_replacement_table[,3]
	#make all TKpseudo black
	genome_name_replacement_table[grep("TK.Pseudo", genome_name_replacement_table[,2]),3]<-"#6715ED"
	#make all Local_Eyach black
	genome_name_replacement_table[grep("Local_Eyach", genome_name_replacement_table[,2]),3]<-"#6715ED"
	#make all "OTU5" lightgrey
	genome_name_replacement_table[grep("ASV01Local_Eyach", genome_name_replacement_table[,2]),3]<-"#15D0ED"
	genome_name_replacement_table[grep("ASV01Local_Roger", genome_name_replacement_table[,2]),3]<-"#15D0ED"

	#assign more colors based on ASV sequence
	unique_seqs=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/sphingo_16Sseqs_to_color_20190909.txt", sep="\t", comment.char="")
	culture_fasta=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/KarasovAlmario16S/all_cultured_16S_20190410_reoriented_trimmed.fa")
	as.matrix(culture_fasta)->culture_fasta
	culture_fasta[sort(c(grep(">", culture_fasta), grep(">", culture_fasta)+1))]->culture_fasta
	split_FASTA(culture_fasta)->culture_fasta
	gsub(">", "", culture_fasta)->culture_fasta
	gsub("GCF_", "GCF", culture_fasta)->culture_fasta
	gsub("_.*", "", culture_fasta)->culture_fasta
	gsub(".Sphingo", "", culture_fasta)->culture_fasta
	gsub(".Pseudo", "", culture_fasta)->culture_fasta
	gsub("GCF", "GCF_", culture_fasta)->culture_fasta
	as.matrix(unique_seqs[match(culture_fasta[,2], unique_seqs[,1]),4])->culture_fasta[,3]
	culture_fasta->genome_name_to_color_key
	genomes_locations_to_give_ASV_color=which(match(genome_name_replacement_table[,1], genome_name_to_color_key[,1], nomatch=0)!=0)
	colors_to_give_them=genome_name_to_color_key[match(genome_name_replacement_table[,1], genome_name_to_color_key[,1], nomatch=0),3]
	colors_to_give_them->genome_name_replacement_table[genomes_locations_to_give_ASV_color,3]
	#fix REFSEQ blacks
	genome_name_replacement_table[grep("REFSEQ",genome_name_replacement_table[,2]),3]<-"#000000"
	
}




#match rownames to the appropriate colors
genome_name_replacement_table[match(rownames(combined), genome_name_replacement_table[,1]),3]->genome_colors



#extract just local strains for the main text
	#keep=c(grep("TK.Pseudo", genome_name_replacement_table[,2]), grep("Local_Eyach", genome_name_replacement_table[,2]))
	#genome_name_replacement_table[keep,]->genome_name_replacement_table
	#combined[keep,]->combined
	#genome_colors[keep]->genome_colors
#Mark the main Pseudomonas ASVs in a special way
if("mark_ASV1"=="ma22"){
	unique(genome_name_replacement_table[,2])
	sort(grep("TK.Pseudo", unique(genome_name_replacement_table[,2]), value=TRUE))->TK_ASVs
	#rep(c("lightgrey", "black"), 100)[1:length(TK_ASVs)]->names(TK_ASVs)
	sort(grep("Eyach", unique(genome_name_replacement_table[,2]), value=TRUE))->EY_ASVs
	#rep(c("lightgrey", "black"), 100)[1:length(EY_ASVs)]->names(EY_ASVs)
	sort(c(TK_ASVs, EY_ASVs))->allASV
	#assign alternating colors in these ASVs
		for(a in 101:199){
			if(a%%2==1){"black"->col}
			if(a%%2==0){"lightgrey"->col}
			substr(a, nchar(a)-1, nchar(a))->a
			names(allASV)[grep(paste("ASV",a,sep="", collapse=""), allASV)]<-col
		} 	
	names(allASV)[match(genome_name_replacement_table[,2], allASV)]->genome_name_replacement_table[,3]

	genome_name_replacement_table[match(rownames(combined), genome_name_replacement_table[,1]),3]->genome_colors
	#ignore heirarchical clustering and order heatmap rows by categories in the order above.
		reorder_vector=order(match(genome_name_replacement_table[,2],allASV))
		combined=combined[reorder_vector,]
		genome_colors[reorder_vector]->genome_colors
}

#SAMPLE SEASON COLORS
colnames(combined) -> rnametest
gsub(".*S[0-9]+.*", "Summer", rnametest)->rnametest
gsub(".*Og[0-9]+.*", "Spring", rnametest)->rnametest
c("Spring", "Summer")->categories
c("cadetblue3", "firebrick3")->names(categories)
cbind(colnames(combined), rnametest)->sample_name_replacement_table
#add colors to sample name replacement table
cbind(sample_name_replacement_table, names(categories)[match(sample_name_replacement_table[,2], categories)])->sample_name_replacement_table
	#ignore heirarchical clustering and order heatmap rows by categories in the order above.
reorder_vector=order(match(sample_name_replacement_table[,2],categories))
combined=combined[,reorder_vector]
#must also reorder the sample replacement table by the same factor? no...? but doesn't hurt
sample_name_replacement_table=sample_name_replacement_table[reorder_vector,]

sample_name_replacement_table[match(colnames(combined), sample_name_replacement_table[,1]),3]->season_colors



plantcols <- rep('black', length(colnames(combined)))

#plant_color_list<-matrix(c(
#	'Athaliana', "#AACF37",
#    'Moss', "#2078B4", 
#    'FuzzyRosemary', "#F59999",
#    'Soil', "#9F4222", 
#    'SmallClover', "#ACACAC", 
#    'Draba', "#3DB549", 
#    'SpringKraut', "#DA4699", 
#    'Cardamine', "#0C6D38", 
#    'Plantago', "#FF0000",
#    'Thistle', "#FF8B00", 
#    'BigClover', "#727272",
#    'Grass', "#5A00EE",
#    'Dandelion', "#fbe955",
#    "?_S4.2_Pseudobulk", "red",  
#    'Clover', "#ACACAC",
#    'NeckarKing', "#000000"
#), ncol=2, byrow=TRUE)

for(i in 1:nrow(plant_color_list)){
  plantcols[grepl(plant_color_list[i,1], colnames(combined))] <- plant_color_list[i,2]
}


#here 2022 I'm finding which decoys have the most reads
classifications=genome_name_replacement_table[match(names(rowSums(combined)), genome_name_replacement_table[,1]),2]
classifications=cbind(names(rowSums(combined)), rowSums(combined), classifications)
classifications=as.data.frame(classifications)
as.numeric(as.matrix(classifications[,2]))->classifications[,2]
answer=classifications[which(classifications[,3]=="DECOY"),]

#here 2022 I'm making a chart of overall mapping to decoy vs local vs refseq
#pdf(file=paste(location_to_write_finished_figures, genus, "_", "SpringSummer_maptoTREE_REF_DECOY_", genus, date, ".pdf", sep="", collapse=""), width = 5, height = 5, useDingbats=FALSE)
	tr=sum(classifications[which(classifications[,3]=="TREE"),2])/10^6
	re=sum(classifications[which(classifications[,3]=="REFSEQ"),2])/10^6
	d=sum(classifications[which(classifications[,3]=="DECOY"),2])/10^6
	plot(c(tr, re, d), type="h", lwd=10, lend=2, ylim=c(0, 400))
#dev.off()

d/(tr+re)   # pseudomonas is 2%  sphingoonas 30%
#





combined->ct

#2022 added normalization here. 
ct=normalize100(ct)

genome_colors->gc
genome_name_replacement_table[,2]->gc2

#PUT BLANKS BETWEEN ROWS
blankwidth=4
#unique(genome_colors)->ucols   #to split columns by color
unique(genome_name_replacement_table[,2])->ucols   #to split columns by category
for(co in 1:length(ucols)){
	print(co)
	#max(which(gc==ucols[co]))->startrow
	max(which(gc2==ucols[co]))->startrow
	print(startrow)
	if(co!=length(ucols)){
		ct=rbind(
			ct[1:startrow,],
			matrix(data=-1, ncol=ncol(ct), nrow=blankwidth),
			ct[(startrow+1):nrow(ct),]
		)
		gc=c(gc[1:startrow], rep("white", blankwidth), gc[(startrow+1):length(gc)])	
		gc2=c(gc2[1:startrow], rep("white", blankwidth), gc2[(startrow+1):length(gc2)])	
	}
}


#PUT BLANKS BETWEEN COLUMNS
gsub("_Og.*", "_Og", colnames(ct))->pc
gsub("_S.*", "_S", pc)->pc
season_colors->sc
blankwidth=1
unique(pc)->ucols
for(co in 1:length(ucols)){
	print(co)
	max(which(pc==ucols[co]))->startcol
	print(startcol)
	if(co!=length(ucols)){
		ct=cbind(
			ct[,1:startcol],
			matrix(data=-1, nrow=nrow(ct), ncol=blankwidth),
			ct[,(startcol+1):ncol(ct)]
		)
		sc=c(sc[1:startcol], rep("white", blankwidth), sc[(startcol+1):length(sc)])	
		pc=c(pc[1:startcol], rep("white", blankwidth), pc[(startcol+1):length(pc)])	
		plantcols=c(plantcols[1:startcol], rep("white", blankwidth), plantcols[(startcol+1):length(plantcols)])	
		
	}
}

#PUT WHITE SPACE BETWEEN SEASONS
blankwidth=20
unique(pc)->ucols
max(grep(ucols[max(grep("_Og", ucols))], pc))->upper
min(grep(ucols[min(grep("_S", ucols))], pc))->lower
mean(c(upper, lower), round=0)->startcol
		ct=cbind(
			ct[,1:startcol],
			matrix(data=-1, nrow=nrow(ct), ncol=blankwidth),
			ct[,(startcol+1):ncol(ct)]
		)
	plantcols=c(plantcols[1:startcol], rep("white", blankwidth), plantcols[(startcol+1):length(plantcols)])





date=format(Sys.Date(), format="%Y%m%d")
sqrt(sqrt(ct))->ct2   #will get warning.. .it's OK. 
-1->ct2[is.na(ct2)]      #if you don't use this, the spacers will be white


 quantile( ct2[ct2>0], prob=c(0, .2, .4, .6, .8, 1))
 
 ct2[ct2>0]


pal1 <- c(blankcolor, colorRampPalette(c("black","#3f324f","#6A1B9A","#E85285", "#FFECB3"))(n=300), "white")



pdf(file=paste(location_to_write_finished_figures, genus, "_", "SpringSummer_heatmap_", date, ".pdf", sep="", collapse=""), width = 5, height = 5, useDingbats=FALSE)
tmp.pdf<- heatmap.2(t(ct2), density.info = "none",
                    trace = "none", dendrogram="none", colRow = plantcols,
                    col = pal1, main = paste(genus, " reads per genome", sep="", collapse=""), 
                    margins = c(10,9), Colv=FALSE, Rowv=FALSE, cexRow = 0.41, 
                    keysize = 0.7,breaks = c(-1, seq(0, 1.77828, length=301), 100),
                    key.title = NA, key.xlab = NA, sepcolor=NA,
                    ColSideColors= gc, RowSideColors=plantcols,
                    labCol = F, sepwidth=c(0,0),
                    key.par=list(mar=c(2.5,0,0,30))   )
dev.off()



cbind(c(-1 , seq(0, 1.77828, length=301), 100), c(pal1, "end"))-> colorkey

pal1[c(2, 50, 100, 150, 200, 250, 302)]->pal2

as.numeric(colorkey[c(2, 50, 100, 150, 200, 250, 302),1])^4

pdf(file=paste(location_to_write_finished_figures, genus, "_", "SpringSummer_heatmap_", date, "_COLORKEY.pdf", sep="", collapse=""), width = 5, height = 5, useDingbats=FALSE)
image(1:length(pal2), 1, as.matrix(1:length(pal2)), 
   	   col=pal2, xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
dev.off()

pdf(file=paste(location_to_write_finished_figures, genus, "_", "SpringSummer_localTREE_", date, ".pdf", sep="", collapse=""), width = 5, height = 10, useDingbats=FALSE)
	plot(core, cex=0.2, show.tip.label = TRUE)
dev.off()


pdf(file=paste(location_to_write_finished_figures, genus, "_", "LEGEND_", date, ".pdf", sep="", collapse=""), width = 5, height = 5, useDingbats=FALSE)


heatmap.2(t(ct2), density.info = "none",
                    trace = "none", dendrogram="none", colRow = plantcols,
                    col = pal1, main = paste(genus, " reads per genome", sep="", collapse=""), 
                    margins = c(10,9), Colv=FALSE, Rowv=FALSE, cexRow = 0.41, 
                    breaks = c(-1, -0.1 , seq(0, 1, length=300)),
                    keysize = 4,
                    key.title = "hi", sepcolor=NA,
                    ColSideColors= gc, RowSideColors=plantcols,
                    labCol = F, sepwidth=c(0,0))
                dev.off()
                    ,
                    key.xtickfun=function() {
               	cex <- par("cex")*par("cex.axis")
               side <- 1
               line <- 0
               col <- par("col.axis")
               font <- par("font.axis")
               mtext("low", side=side, at=0, adj=0,
                     line=line, cex=cex, col=col, font=font)
               mtext("high", side=side, at=1, adj=1,
                     line=line, cex=cex, col=col, font=font)
               return(list(labels=FALSE, tick=FALSE))
          }
          )
dev.off()


mean(ct2[1:49,grep("lover", colnames(ct2))] ) / mean(ct2[59:125,grep("lover", colnames(ct2))])  
mean(ct2[1:49,grep("thaliana", colnames(ct2))] ) / mean(ct2[59:125,grep("thaliana", colnames(ct2))])  
mean(ct2[1:49,grep("Draba", colnames(ct2))] ) / mean(ct2[59:125,grep("Draba", colnames(ct2))])  
mean(ct2[1:49,grep("ardamine", colnames(ct2))] ) / mean(ct2[59:125,grep("ardamine", colnames(ct2))])  
mean(ct2[1:49,grep("andelion", colnames(ct2))] )/mean(ct2[59:125,grep("andelion", colnames(ct2))]) 




ct3=ct
ct3[ct3<0]=0


#Histograms for each genotype
gsub("_Og.*", "_Og", colnames(ct3))->pc3
gsub("_S.*", "_S", pc3)->pc3
unique(pc3)->ucols
ucols[which(ucols!="")]->ucols
grep("_Og", ucols, value=TRUE)->ucols
#pdf(file="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/PseudoSpringSummer_hist.pdf", width = 5, height = 5)
par(mfrow=c(10,1))
par(mar=c(.5,.5,.5,.5))
for(co in 1:length(ucols)){
	matrix_format(ct3[,which(pc3==ucols[co])])->each_plant
	plot(rowSums(each_plant), type="h", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")	
}
#dev.off()

#Histograms for spring vs. summer
gsub(".*Og.*", "Og", colnames(ct3))->pc3
gsub(".*S.*", "S", pc3)->pc3
ct3[1:length(grep("TREE", genome_name_replacement_table[,2])),]->ct3   #if only the tree
unique(pc3)->ucols
ucols[which(ucols!="")]->ucols
grep("Og|S", ucols, value=TRUE)->ucols
pdf(file=paste(location_to_write_finished_figures, genus, "_", "SpringSummer_sum_histogram_", date, ".pdf", sep="", collapse=""), width = 5, height = 5, useDingbats=FALSE)
par(mfrow=c(length(unique(ucols)),1))
par(mar=c(.5,.5,.5,.5))
for(co in 1:length(ucols)){   #first time through determinate max sum to set consistent X axis
	print(co)
	matrix_format(ct3[,which(pc3==ucols[co])])->each_plant
	plot(rowSums(each_plant), type="h", xaxt="n", yaxt="n", xlab="", ylab="", bty="n", lwd=2)	
}
dev.off()

 
 
 
 
ct3[,grep("_S", colnames(ct3))]->ct3summer
ct3[,grep("_Og", colnames(ct3))]->ct3spring


for(uc in unique(genome_colors)){
	print(round(100*sum(ct3spring[which(uc==genome_colors),])/sum(ct3spring),1))
}
for(uc in unique(genome_colors)){
	print(round(100*sum(ct3summer[which(uc==genome_colors),])/sum(ct3summer),1))
}




if(genus=="Pseudomonas"){histcol="#3281BC"}
if(genus=="Sphingomonas"){histcol="#DD2E27"}

pdf(file=paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/", genus,"_", date, "Spring_hist.pdf", sep="", collapse=""), width = 5, height = 2)
plot(rowSums(ct3spring),type="h", xaxt="n", yaxt="n", xlab="", ylab="", bty="n", col=histcol, lend="square", lwd=0.4)	
dev.off()
pdf(file=paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Bulk_metagenomes/", genus, "_", date,"Summer_hist.pdf", sep="", collapse=""), width = 5, height = 2)
plot(rowSums(ct3summer),type="h", xaxt="n", yaxt="n", xlab="", ylab="", bty="n", col=histcol, lend="square", lwd=0.4)				
dev.off()



c(51, 204)*2/3



}   #TO CLEAN UP VIEWING AREA AND SHRING THIS STUFF...


#de only comparisons of the same strains as in m
	#2022 m was defined quite some distance above, and is the euclidean distance in the presence absense matrix of accessory genes

	fullmatrix_Sph=fullmatrix_Sph[,match(unique(m[,1]), colnames(fullmatrix_Sph))]
	fullmatrix_Sph=fullmatrix_Sph[match(unique(m[,1]), rownames(fullmatrix_Sph)),]
	mash_m=melt(fullmatrix_Sph)
	#paste the first and second columns together to make a unique comparison string
	mash_m$comparison_string=apply(mash_m[,1:2], 1, paste, collapse="_", sep="_")
	mash_m=mash_m[which(mash_m[,3]!=1),]
	
	
	mash_m=mash_m[match(m$comparison_string, mash_m$comparison_string, nomatch=0),]
	m=m[match(mash_m$comparison_string, m$comparison_string, nomatch=0),]
	
	
	
	plot(as.numeric(as.matrix(mash_m[,3])), as.numeric(as.matrix(m[,3])))
	plot(sqrt(as.numeric(as.matrix(mash_m[,3]))), sqrt(as.numeric(as.matrix(m[,3]))))
	
	#make dataframe combining presence absence distance and mash distance
	PA_vs_mash=cbind(m, mash_m)
	PA_vs_mash=as.data.frame(PA_vs_mash)
	names(PA_vs_mash)=c("PA_c1", "PA_c2", "PA_d", "PA_string","MASH_c1", "MASH_c2", "MASH_d", "MASH_string")
	
	pdf(file=paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PNAS_revision/", "presenceAbsence_vs_MASH_", date, ".pdf", sep="", collapse=""), width = 5, height = 5, useDingbats=FALSE)
		plot(PA_vs_mash$MASH_d, PA_vs_mash$PA_d, cex=1, xlab="MASH",  ylab="PA", col=add.alpha("black", alpha=0.03), pch=16,
			xaxt="n", yaxt="n", xlim=c(0, .3), ylim=c(0, .8))	
		axis(side=1, at=c(0, .1, .2, .3), labels=c(0, .1, .2, .3))
		axis(side=2, at=c(0, .2, .4, .6, .8), labels=c(0, .2, .4, .6, .8))
		abline(lm(PA_vs_mash$PA_d~PA_vs_mash$MASH_d), lwd=2, col="red")
		cor(as.numeric(as.matrix(PA_vs_mash$PA_d)), as.numeric(as.matrix(PA_vs_mash$MASH_d)))^2
	dev.off()	
	


	
	mantel.rtest(as.numeric(as.matrix(mash_m[,3])),  as.numeric(as.matrix(m[,3])), nrepet = 999)
	
	
	as.dist(as.numeric(as.matrix(mash_m[,3])))
	
	my_nj <- ape::nj(testDIST)
	plot(my_nj,  cex=0.2)
	
	gsub("_nanopore.*", "", core$tip.label)->core$tip.label 
	Drop=setdiff(core$tip.label, my_nj$tip.label)
	Drop=match(Drop, core$tip.label)
	drop.tip(core, as.numeric(Drop))->core
	
	
	
	
	install.packages("TreeDist")
	library('TreeDist')


	setdiff(core$tip.label, my_nj$tip.label)
	distance <- TreeDistance(core, my_nj)
	TreeDistance(core, core)
	
	VisualizeMatching(ClusteringInfoDistance, core, my_nj, cex=.2)
	
	
	pal <- c(colorRampPalette(c("black","#27004E", "#340069","#6A1B9A", "#CB018E", "#EB6E99", "#FFBC7A", "#FFD6A9"))(n=8))
	pal1=pal[length(pal):1]
	library(gplots)
	
	#This will totally stall a laptop...too much processing for a big matrix
   heatmap.2(PA[nrow(PA):1,], density.info = "none",
                    trace = "none", dendrogram="none", main="hi",
                    col = pal1, Colv=FALSE, Rowv=FALSE, key=FALSE, 
                    key.title = NA, key.xlab = NA, sepcolor="#FFFFFF",
                    breaks = c(0, .30, .40, .50, .60, .70, .80, .90, 100),
                    labCol = F, labRow=F, sepwidth=c(0,.01), 
                    colsep=NA, rowsep=1:nrow(PA))      

	
#}





#PA->PAsaved2
#gene_group_colors->gene_group_colors2
#gbc->gbc2

PA<-PAsaved2
gene_group_colors<-gene_group_colors2
gbc<-gbc2

#condensed can use heatmap
if(subset==TRUE){
	P=PA
	A=PA+1
	A[which(A==2)]=0
	#purge rows where only a couple strains have it or are missing it
	#gene_group_colors=gene_group_colors[which(colSums(PA)<=nrow(PA)-5)]
	#	PA=PA[,which(colSums(PA)<=nrow(PA)-5)]
	#gene_group_colors=gene_group_colors[which(colSums(PA)>=5  |  colSums(PA)< -1)]
	#	PA=PA[,which(colSums(PA)>=5  |  colSums(PA)< -1)]



	#extract specific groups
	matrix_format(PA[,grep("xanthine synthesis",gbc)], type="V")->PA
	gene_group_colors[grep("xanthine synthesis", gbc)]->gene_group_colors
	cbind(PA, rep(-1, nrow(PA)))->PA
	c(gene_group_colors, "white")->gene_group_colors
	
	
	library(gplots)
	pal1 <- c("white", "#FBCE33","#3f324f")
	tmp.pdf<- heatmap.2(PA[nrow(PA):1,], density.info = "none",
                    trace = "none", dendrogram="none", 
                    col = pal1, Colv=FALSE, Rowv=FALSE,  
                    ColSideColors=gene_group_colors,
                    key.title = NA, key.xlab = NA, sepcolor=NA,
                    labCol = F, sepwidth=c(0,0),
                    )
                    
                    
                    
}


2B2461 = blue
FBCE33 = yellow





duplicates=c("S165H113","S67H113","S167H113","S169H113","S170H113","S171H113","S168H113","S74H113","S172H113","S372H113")      
setdiff(1:length(rownames(PA)), match(duplicates, rownames(PA), nomatch=0))->keep

PA2[keep,,]->PA2

plot.new() 
rasterImage(PA2[,1:1000,], 0,0,1,1)   
          
        

          
          
          
          
          
          
          
          
