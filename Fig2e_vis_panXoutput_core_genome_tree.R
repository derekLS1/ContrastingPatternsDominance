#visualize core genome tree in APE
#install.packages("ape")
#install.packages("phangorn")
#install.packages("phytools")
#install.packages("geiger")
#install.packages("ape")
#install.packages("Biostrings")
#install.packages("ggplot2")
#install.packages("gplots")

library("ape")
library("phangorn")
library("phytools")
library("geiger")
library("ape")
library("Biostrings")
library("ggplot2")
library("gplots")
#library("ggtree")
source("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/scripts/R/functions/microbiome_custom_functions.R")



#core=read.tree(file = "/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sphingo_ASV1ASV4/strain_tree.nwk")



#choose genus
genus="Sphingomonas"
date=format(Sys.Date(), format="%Y%m%d")
nanopore_control=TRUE #TRUE
otu99=FALSE #TRUE  #TRUE

if(genus=="Sphingomonas"){
	#basefolder="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_buscoGT85/"
	#basefolder="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ/"
	#basefolder="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_lessDivergent/"
	basefolder="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS/"
	#basefolder="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS_ASV1/"
	core=read.tree(file = paste(basefolder, "strain_tree.nwk", sep="", collapse=""))
}

length(core$tip.label)

if(genus=="Pseudomonas"){ 
	basefolder="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Pseudo_runs/Pseudo_BigRun_20190510/"
	#basefolder="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Pseudo_runs/Pseudo_ASV1/"
	core=read.tree(file = paste(basefolder, "strain_tree.nwk", sep="", collapse=""))
}

#LOAD METADATA
metadata=read.table(file="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Single_genomes/core_genome_tree_metadata.txt", comment.char="", quote = "")
as.matrix(metadata[1,])->names(metadata)
metadata[2:nrow(metadata),]->metadata


#COLOR GENOME NAME BY OTU IDENTITY
#LOAD CULTRE FASTA AND EXTRACT SPHINGO SEQS
culture_fasta=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/KarasovAlmario16S/all_cultured_16S_20190410_reoriented_trimmed.fa")
if(otu99==TRUE){
	culture_fasta=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/KarasovAlmario16S/cluster_cultured_to_99pcnt/all_cultured_16S_20190818_reoriented_trimmed_FL_995otus.fa") 
	#subset by those which have the SphASV1 sequence

}
	
	
as.matrix(culture_fasta)->culture_fasta

if(genus=="Sphingomonas"){ 
	#drop those also duplicated by nanopore sequencing
		if(nanopore_control==FALSE){
			Drop=which(match(core$tip.label,c("S337H113","S132H113","S136H113",
										"S380H113","S237H113","S230H113",
										"S190H113","S213H113","S18H113",
										"S127H113","S133H113","S216H113"), nomatch=0)>0)
			drop.tip(core, as.numeric(Drop))->core
			if(basefolder=="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS_ASV1/"){
				drop.tip(core, which(core$tip.label=="S150H133"))->core
			}
		}
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
		}

	
	print(length(core$tip.label))
		
	if(otu99!=TRUE){
		culture_fasta[c(
			sort(c(grep("Sphingo", culture_fasta), (grep("Sphingo", culture_fasta)+1))),
			sort(c(grep("GCF_", culture_fasta), (grep("GCF_", culture_fasta)+1)))),]->Sphingo16S_culture
			split_FASTA(Sphingo16S_culture, add_integer_column=FALSE, leave_carrot=FALSE, reverse=FALSE)->Sphingo16S_table
			gsub("H133.Sphingo_16S", "H133" , Sphingo16S_table[,1])->Sphingo16S_table[,1]
			gsub("H113.Sphingo_16S", "H113" , Sphingo16S_table[,1])->Sphingo16S_table[,1]
			gsub("Sphingo", "" , Sphingo16S_table[,1])->Sphingo16S_table[,1]
			gsub("_16S", "" , Sphingo16S_table[,1])->Sphingo16S_table[,1]
			gsub("GCF_", "GCF" , Sphingo16S_table[,1])->Sphingo16S_table[,1]
			gsub("_.*", "" , Sphingo16S_table[,1])->Sphingo16S_table[,1]
			gsub("GCF", "GCF_" , Sphingo16S_table[,1])->Sphingo16S_table[,1]
	}
	if(otu99==TRUE){
		culture_fasta[c(
			sort(c(grep("S.*H.*", culture_fasta), (grep("S.*H.*", culture_fasta)+1))),
			sort(c(grep("GCF_", culture_fasta), (grep("GCF_", culture_fasta)+1)))),]->Sphingo16S_culture
			split_FASTA(Sphingo16S_culture, add_integer_column=FALSE, leave_carrot=FALSE, reverse=FALSE)->Sphingo16S_table
				SphingoASV1="TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTTACCCGGGATGATAATGACAGTACCGGGAGAATAAGCTCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGAGCTAGCGTTATTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTAAGAGGTGAAAGCCCAGAGCTCAACTCTGGAATTGCCTTTTAGACTGCATCGCTTGAATCATGGAGAGGTCAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGGCTGACTGGACATGTATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG"
			#subset the full culture fasta of 16S by only those sharing hte V3V4 ASV1 sequence
			Sphingo16S_table[grep(SphingoASV1, Sphingo16S_table[,2]),]->Sphingo16S_table
	}
	#Turn tip labels into the actual V3V4 ASV sequences
	#make temp core$tip.label that renames the nanopore seqs
	gsub("Sph.*_S", "S", core$tip.label)->tmp.tip.label
	gsub("_nanopore.*", "", tmp.tip.label)->tmp.tip.label
	Sphingo16S_table[match(tmp.tip.label, Sphingo16S_table[,1]),2]->associated_seqs_of_rownames
}



	#20220525 for some reason this tree is not only ASV1. Therefore need to pull the correct strains out. 
	#drop those which are not among the 174 SphASV1
	if(genus=="Sphingomonas"){
	if(basefolder=="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS_ASV1/"){
	print("hi")
	keep=grep("TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTTACCCGGGATGATAATGACAGTACCGGGAGAATAAGCTCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGAGCTAGCGTTATTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTAAGAGGTGAAAGCCCAGAGCTCAACTCTGGAATTGCCTTTTAGACTGCATCGCTTGAATCATGGAGAGGTCAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGGCTGACTGGACATGTATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG", associated_seqs_of_rownames)
	throw=grep("TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTTACCCGGGATGATAATGACAGTACCGGGAGAATAAGCTCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGAGCTAGCGTTATTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTAAGAGGTGAAAGCCCAGAGCTCAACTCTGGAATTGCCTTTTAGACTGCATCGCTTGAATCATGGAGAGGTCAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGGCTGACTGGACATGTATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG", associated_seqs_of_rownames, invert=TRUE)
	drop.tip(core, throw)->core	
	
	throw=grep("GCF_", associated_seqs_of_rownames)
	keep=grep("GCF_", associated_seqs_of_rownames, invert=TRUE)
	drop.tip(core, throw)->core	
	
	gsub("Sph.*_S", "S", core$tip.label)->tmp.tip.label
	gsub("_nanopore.*", "", tmp.tip.label)->tmp.tip.label
	Sphingo16S_table[match(tmp.tip.label, Sphingo16S_table[,1]),2]->associated_seqs_of_rownames
	
	}}
	




if(genus=="Pseudomonas"){ 
	#drop two outlier genomes not similar to mine
	#Drop=which(match(gsub("_.*", "", core$tip.label),c("p9.H8","p9.H9"), nomatch=0)>0)
	#drop.tip(core, as.numeric(Drop))->core

	culture_fasta[
		c(sort(c(grep("Pseudo", culture_fasta), (grep("Pseudo", culture_fasta)+1))), 
		    sort(c(grep(">p", culture_fasta), (grep(">p", culture_fasta)+1)))),   ]->Pseudo16S_culture
		split_FASTA(Pseudo16S_culture, add_integer_column=FALSE, leave_carrot=FALSE, reverse=FALSE)->Pseudo16S_table
		gsub("H133.Pseudo_16S", "H133" , Pseudo16S_table[,1])->Pseudo16S_table[,1]
		gsub("H113.Pseudo_16S", "H113" , Pseudo16S_table[,1])->Pseudo16S_table[,1]
		gsub("Pseudo", "p" , Pseudo16S_table[,1])->Pseudo16S_table[,1]
		gsub("_16S", "" , Pseudo16S_table[,1])->Pseudo16S_table[,1]
	#Turn tip labels into the actual V3V4 ASV sequences
	gsub("_.*", "", core$tip.label)->core$tip.label
	Pseudo16S_table[match(core$tip.label, Pseudo16S_table[,1]),2]->associated_seqs_of_rownames
}


#Rename cultured seqs into ASVs
as.matrix(table(associated_seqs_of_rownames))->unique_seqs
unique_seqs[order(unique_seqs[,1], decreasing=TRUE),]->unique_seqs
cbind(names(unique_seqs), unique_seqs)->unique_seqs
rownames(unique_seqs)<-NULL
unique_seqs=cbind(unique_seqs, paste("ASV", 1:nrow(unique_seqs), sep=""))

#calculate how many are unique looking at 127 bp of the front and back (HiSeq)
first127=apply(do.call(rbind, lapply(strsplit(associated_seqs_of_rownames, split=""), head, n=127)), 1, paste, collapse="")
last127=apply(do.call(rbind, lapply(strsplit(associated_seqs_of_rownames, split=""), tail, n=127)), 1, paste, collapse="")
HiSeq_resolution=apply(cbind(first127, last127), 1, paste, collapse="")
SynCom_options=cbind(core$tip.label, associated_seqs_of_rownames, HiSeq_resolution)
write.table(SynCom_options, "/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/random_design_tables/SynCom_options.txt", col.names=FALSE, row.names=FALSE, sep="\t")

#calculate V4 trimmed uniques
associated_seqs_of_rownamesV4=gsub(".*GTG.CAGC.GCCGCGGTAA", "", associated_seqs_of_rownames)
first127V4=apply(do.call(rbind, lapply(strsplit(associated_seqs_of_rownamesV4, split=""), head, n=127)), 1, paste, collapse="")
last127V4=apply(do.call(rbind, lapply(strsplit(associated_seqs_of_rownamesV4, split=""), tail, n=127)), 1, paste, collapse="")
HiSeq_resolutionV4=apply(cbind(first127V4, last127V4), 1, paste, collapse="")
SynCom_optionsV4=cbind(core$tip.label, associated_seqs_of_rownamesV4, HiSeq_resolutionV4)
write.table(SynCom_optionsV4, "/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/random_design_tables/SynCom_optionsV4.txt", col.names=FALSE, row.names=FALSE, sep="\t")



#SynCom_optionsV4[3,3]
#SynCom_optionsV4[3,4]


library(RColorBrewer)
n <- nrow(unique_seqs)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
if(genus=="Sphingomonas"){
	if(otu99!=TRUE){
		col_vector[17:length(col_vector)]->col_vector
		"yellow"->col_vector[11]
		gsub("#A6CEE3", "#B71E57", col_vector)->col_vector
	}
	if(otu99==TRUE){
		#col_vector[25:length(col_vector)]->col_vector
		#"yellow"->col_vector[11]
		library(viridis)
		sample(c(inferno(40), viridis(40)), 80)->col_vector
		#col_vector=c("#E3E418FF","#E35932FF","#440154FF","#781C6DFF","#472D7BFF","#110A31FF","#9A2865FF","#27AD81FF","#CA404AFF","#BB3754FF","yellow",,"#31688EFF","#C7E020FF","#330A5FFF","#D84C3EFF","#56106EFF","#21908CFF","#8FD744FF","#20A486FF","#68166EFF","#F98C0AFF","#1F9A8AFF","#450A69FF","#F4DE52FF","#000004FF","#47C16EFF","#24868EFF","#2C728EFF","#FCFFA4FF","#365D8DFF","#481F70FF","#F1F17AFF","#75D054FF","#060418FF","#AB2F5EFF","#210C4AFF","#89226AFF","#443A83FF","#287C8EFF","#AADC32FF","#404688FF","#FCA108FF","#FCB519FF","#FDE725FF","#35B779FF","#ED6925FF","#3B528BFF","#471164FF","#F9C932FF","#5DC863FF")
		#col_vector=c("#460B5DFF","#ED6925FF","#B63458FF","#FCA209FF","#67CC5CFF","#FB9606FF","#297B8EFF","#A22B62FF","#472E7CFF","#F67E14FF","#404588FF","#375B8DFF","#F98A0BFF","#2E6F8EFF","#481467FF","#BADE28FF","#6E186EFF","#C93F4BFF","#24878EFF","#4F0D6CFF","#21A585FF","#482576FF","#3D4D8AFF","#433E85FF","#34618DFF","#3A538BFF","#26818EFF","#87D549FF","#781C6DFF","#2D0B5AFF","#83206BFF","#97D83FFF","#440A68FF","#59106EFF","#59C864FF","#F4E258FF","#C03A51FF","#FCAE12FF","#0F092CFF","#76D153FF","#390963FF","#DDE318FF","#64156EFF","#228D8DFF","#AD305DFF","#A9DB33FF","#EDE51BFF","#F6D543FF","#1FA088FF","#CBE11EFF","#35B779FF","#440154FF","#1F998AFF","#FDE725FF","#E15635FF","#2DB17EFF","#481D6FFF","#453581FF","#08051DFF","#F3F78CFF","#25AC82FF","#F1ED71FF","#982766FF","#31688EFF","#8D2369FF","#230C4BFF","#D24644FF","#FBBC21FF","#180C3CFF","#40BC72FF","#F9C830FF","#F1741CFF","#D94D3DFF","#2B748EFF","#000004FF","#02020FFF","#E75F2EFF","#20938CFF","#FCFFA4FF","#4CC26CFF")
		col_vector=c("#ED6925FF","#460B5DFF","#ffff20","#FCA209FF","#67CC5CFF","#FB9606FF","#297B8EFF","#A22B62FF","#472E7CFF","#F67E14FF","#404588FF","#375B8DFF","#F98A0BFF","#2E6F8EFF","#481467FF","#BADE28FF","#6E186EFF","#C93F4BFF","#24878EFF","#4F0D6CFF","#21A585FF","#482576FF","#3D4D8AFF","#433E85FF","#34618DFF","#3A538BFF","#26818EFF","#87D549FF","#781C6DFF","#2D0B5AFF","#83206BFF","#97D83FFF","#440A68FF","#59106EFF","#59C864FF","#F4E258FF","#C03A51FF","#FCAE12FF","#0F092CFF","#76D153FF","#390963FF","#DDE318FF","#64156EFF","#228D8DFF","#AD305DFF","#A9DB33FF","#EDE51BFF","#F6D543FF","#1FA088FF","#CBE11EFF","#35B779FF","#440154FF","#1F998AFF","#FDE725FF","#E15635FF","#2DB17EFF","#481D6FFF","#453581FF","#08051DFF","#F3F78CFF","#25AC82FF","#F1ED71FF","#982766FF","#31688EFF","#8D2369FF","#230C4BFF","#D24644FF","#FBBC21FF","#180C3CFF","#40BC72FF","#F9C830FF","#F1741CFF","#D94D3DFF","#2B748EFF","#000004FF","#02020FFF","#E75F2EFF","#20938CFF","#FCFFA4FF","#4CC26CFF")
		
		key99=c("AGAGTTTGATCATGGCTCAGAATGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGAAGCCTTCGGGTTTAGTGGCGCACGGGTGCGTAACGCGTGGGAATCTGCCCCTTGGTTCGGAATAACAGTTGGAAACGACTGCTAATACCGGATGATGACGTAAGTCCAAAGATTTATCGCCGAGGGATGAGCCCGCGTAGGATTAGGTAGTTGGTGTGGTAAAGGCGCACCAAGCCGACGATCCTTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTTACCCGGGATGATAATGACAGTACCGGGAGAATAAGCTCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGAGCTAGCGTTATTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTAAGAGGTGAAAGCCCAGAGCTCAACTCTGGAATTGCCTTTTAGACTGCATCGCTTGAATCATGGAGAGGTCAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGGCTGACTGGACATGTATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGATAACTAGCTGTCCGGGGACTTGGTCTTTGGGTGGCGCAGCTAACGCATTAAGTTATCCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCTGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGCAGAACCTTACCAGCGTTTGACATGGCAGGACGACTTCTGGAGACAGATTTCTTCCCTTCGGGGACCTGCACACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTCGACCTTAGTTGCCATCATTTAGTTGGGCACTTTAAGGTAACCGCCGGTGATAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTACGCGCTGGGCTACACACGTGCTACAATGGCGGTGACAGTGGGCAGCAAGCACGCGAGTGTGAGCTAATCTCCAAAAGCCGTCTCAGTTCGGATTGTTCTCTGCAACTCGAGAGCATGAAGGCGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCAGGCCTTGTACACACCGCCCGTCACACCATGGGAGTTGGATTCACCCGAAGGCGTTGCGCTAACTCAGCAATGAGAGGCAGGCGACCACGGTGGGTTTAGCGACTGGGGTGAAGTCGTAACAAGGT",    
	  "AGAGTTTGATCATGGCTCAGAATGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGAATTACCTTCGGGTAGTTAGTGGCGCACGGGTGCGTAACGCGTGGGAATCTGCCCCTTGGTTCGGAATAACAGTTGGAAACGACTGCTAATACCGGATGATGACGTAAGTCCAAAGATTTATCGCCGAGGGATGAGCCCGCGTAGGATTAGGTAGTTGGTGTGGTAAAGGCGCACCAAGCCGACGATCCTTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTTACCCGGGATGATAATGACAGTACCGGGAGAATAAGCTCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGAGCTAGCGTTATTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTAAGAGGTGAAAGCCCAGAGCTCAACTCTGGAATTGCCTTTTAGACTGCATCGCTTGAATCATGGAGAGGTCAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGGCTGACTGGACATGTATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGATAACTAGCTGTCCGGACACTTGGTGTTTGGGTGGCGCAGCTAACGCATTAAGTTATCCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCTGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGCAGAACCTTACCAGCGTTTGACATGGCAGGACGACTTCCAGAGATGGATTTCTTCCCTTCGGGGACCTGCACACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTCGCCTTTAGTTGCCATCATTTAGTTGGGCACTTTAAAGGAACCGCCGGTGATAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTACGCGCTGGGCTACACACGTGCTACAATGGCGGTGACAGTGGGCAGCAAGCACGCGAGTGTGAGCTAATCTCCAAAAGCCGTCTCAGTTCGGATTGTTCTCTGCAACTCGAGAGCATGAAGGCGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCAGGCCTTGTACACACCGCCCGTCACACCATGGGAGTTGGATTCACCCGAAGGCGTTGCGCTAACTCAGCAATGAGAGGCAGGCGACCACGGTGGGTTTAGCGACTGGGGTGAAGTCGTAACAAGGT",
	  "AGAGTTTGATCATGGCTCAGAATGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGAAGTCCTTCGGGGCTTAGTGGCGCACGGGTGCGTAACGCGTGGGAATCTGCCCCTTGGTTCGGAATAACAGTTGGAAACGACTGCTAATACCGGATGATGACGTAAGTCCAAAGATTTATCGCCGAGGGATGAGCCCGCGTAGGATTAGGTAGTTGGTGTGGTAAAGGCGCACCAAGCCGACGATCCTTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTTACCCGGGATGATAATGACAGTACCGGGAGAATAAGCTCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGAGCTAGCGTTATTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTAAGAGGTGAAAGCCCAGAGCTCAACTCTGGAATTGCCTTTTAGACTGCATCGCTTGAATCATGGAGAGGTCAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGGCTGACTGGACATGTATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGATAACTAGCTGTCCGGGGACTTGGTCTTTGGGTGGCGCAGCTAACGCATTAAGTTATCCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCTGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGCAGAACCTTACCAGCGTTTGACATGGCAGGACGACTTCCGGAGACGGATTTCTTCCCTTCGGGGACCTGCACACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTCGCCTTTAGTTGCCATCATTTAGTTGGGCACTTTAAAGGAACCGCCGGTGATAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTACGCGCTGGGCTACACACGTGCTACAATGGCGGTGACAGTGGGCAGCAAGCACGCGAGTGTGAGCTAATCTCCAAAAGCCGTCTCAGTTCGGATTGTTCTCTGCAACTCGAGAGCATGAAGGCGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCAGGCCTTGTACACACCGCCCGTCACACCATGGGAGTTGGATTCACCCGAAGGCGTTGCGCTAACTCGTAAGAGAGGCAGGCGACCACGGTGGGTTTAGCGACTGGGGTGAAGTCGTAACAAGGT",    
	  "AGAGTTTGATCATGGCTCAGAATGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGAGACCTTCGGGTCTAGTGGCGCACGGGTGCGTAACGCGTGGGAATCTGCCCCTTGGTTCGGAATAACAGTTGGAAACGACTGCTAATACCGGATGATGACGTAAGTCCAAAGATTTATCGCCGAGGGATGAGCCCGCGTAGGATTAGGTAGTTGGTGTGGTAAAGGCGCACCAAGCCGACGATCCTTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTTACCCGGGATGATAATGACAGTACCGGGAGAATAAGCTCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGAGCTAGCGTTATTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTAAGAGGTGAAAGCCCAGAGCTCAACTCTGGAATTGCCTTTTAGACTGCATCGCTTGAATCATGGAGAGGTCAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGGCTGACTGGACATGTATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGATAACTAGCTGTCCGGACACTTGGTGTTTGGGTGGCGCAGCTAACGCATTAAGTTATCCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCTGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGCAGAACCTTACCAGCGTTTGACATGGCAGGACGACTTCCAGAGATGGATTTCTTCCCTTCGGGGACCTGCACACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTCGCCTTTAGTTGCCATCATTCAGTTGGGCACTTTAAAGGAACCGCCGGTGATAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTACGCGCTGGGCTACACACGTGCTACAATGGCGGTGACAGTGGGCAGCAAACTCGCAAGGGTGAGCTAATCTCCAAAAGCCGTCTCAGTTCGGATTGTTCTCTGCAACTCGAGAGCATGAAGGCGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCAGGCCTTGTACACACCGCCCGTCACACCATGGGAGTTGGATTCACCCGAAGGCGTTGCGCTAACTCGCAAGAGAGGCAGGCGACCACGGTGGGTTTAGCGACTGGGGTGAAGTCGTAACAAGGT",      
	  "AGAGTTTGATCATGGCTCAGAATGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGAAGGCTTCGGCCTTAGTGGCGCACGGGTGCGTAACGCGTGGGAATCTGCCCTCAGGTTCGGAATAACAGTTGGAAACGACTGCTAATACCGGATGATGACGTAAGTCCAAAGATTTATCGCCTGAGGATGAGCCCGCGTAGGATTAGCTAGTTGGTGTGGTAAAGGCGCACCAAGGCGACGATCCTTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTTACCCGGGATGATAATGACAGTACCGGGAGAATAAGCTCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGAGCTAGCGTTATTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTAAGAGGTGAAAGCCCAGAGCTCAACTCTGGAATTGCCTTTTAGACTGCATCGCTTGAATCATGGAGAGGTCAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGGCTGACTGGACATGTATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGATAACTAGCTGTCCGGACACTTGGTGTTTGGGTGGCGCAGCTAACGCATTAAGTTATCCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCTGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGCAGAACCTTACCAGCGTTTGACATGGCAGGACGACTTCCAGAGATGGATTTCTTCCCTTCGGGGACCTGCACACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTCGACCTTAGTTGCCATCATTTAGTTGGGCACTTTAAGGTAACCGCCGGTGATAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTACGCGCTGGGCTACACACGTGCTACAATGGCGGTGACAGTGGGCAGCAAACTCGCAAGGGTGAGCTAATCTCCAAAAGCCGTCTCAGTTCGGATTGTTCTCTGCAACTCGAGAGCATGAAGGCGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCAGGCCTTGTACACACCGCCCGTCACACCATGGGAGTTGGATTCACCCGAAGGCGTTGCGCTAACTCGCAAGAGAGGCAGGCGACCACGGTGGGTTTAGCGACTGGGGTGAAGTCGTAACAAGGT",      
	  "AGAGTTTGATCATGGCTCAGAATGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGAAGGCTTCGGCCTTAGTGGCGCACGGGTGCGTAACGCGTGGGAATCTGCCCCTTGGTTCGGAATAACAGTTGGAAACGACTGCTAATACCGGATGATGACGTAAGTCCAAAGATTTATCGCCGAGGGATGAGCCCGCGTAGGATTAGCTAGTTGGTGGGGTAAAGGCCTACCAAGGCGACGATCCTTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTTACCCGGGATGATAATGACAGTACCGGGAGAATAAGCTCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGAGCTAGCGTTATTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTAAGAGGTGAAAGCCCAGAGCTCAACTCTGGAATTGCCTTTTAGACTGCATCGCTTGAATCATGGAGAGGTCAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGGCTGACTGGACATGTATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGATAACTAGCTGTCCGGGCACTTGGTGCCTGGGTGGCGCAGCTAACGCATTAAGTTATCCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCTGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGCAGAACCTTACCAGCGTTTGACATGGCAGGACGACTTCCAGAGATGGATTTCTTCCCTTCGGGGACCTGCACACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTCGCCTTTAGTTGCCATCATTTAGTTGGGCACTTTAAAGGAACCGCCGGTGATAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTACGCGCTGGGCTACACACGTGCTACAATGGCGGTGACAGTGGGCAGCAAGCACGCGAGTGTGAGCTAATCTCCAAAAGCCGTCTCAGTTCGGATTGTTCTCTGCAACTCGAGAGCATGAAGGCGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCAGGCCTTGTACACACCGCCCGTCACACCATGGGAGTTGGATTCACCCGAAGGCGTTGCGCTAACTCAGCAATGAGAGGCAGGCGACCACGGTGGGTTTAGCGACTGGGGTGAAGTCGTAACAAGGT",
		"#ED6925FF", "#460B5DFF", "#ffff20",   "#FCA209FF", "#67CC5CFF", "#FB9606FF")
		matrix(key99, ncol=2)->key99
		
		key99[match(unique_seqs[,1], key99[,1]),2]->col_vector
	}
}
if(genus=="Pseudomonas"){
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
   #maxColors[17:80]->col_vector
   #col_vector[11:length(col_vector)]<-"#FFFFFF"
   col_vector[17:length(col_vector)]->col_vector
	"yellow"->col_vector[11]
	gsub("#A6CEE3", "#B71E57", col_vector)->col_vector
      col_vector[setdiff(1:length(col_vector), c(4, 15:43))]->col_vector
       "#77C043"->col_vector[1]
       gsub("#B2DF8A", "#455636", col_vector)->col_vector
       gsub("#B3B3B3", "#C027CC", col_vector)->col_vector
       gsub("#D9D9D9", "#0034A8", col_vector)->col_vector
       c(col_vector, c("#297B8EFF","#A22B62FF","#472E7CFF","#A9DB33FF","#1F998AFF","#FDE725FF"))->col_vector
	image(1:length(col_vector), 1, as.matrix(1:length(col_vector)), 
      col=col_vector, xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
}

#image(1:length(col_vector), 1, as.matrix(1:length(col_vector)), 
 #     col=col_vector, xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
#Dpal<-colorRampPalette(colors=c("#41C4E2", "#6050A1"))
#Dpal2<-colorRampPalette(colors=c("#C13AA6", "#6050A1"))
#Dpal3<-colorRampPalette(colors=c("#C1A03A", "#3AC19A"))
#col_vector=c("#41C4E2", Dpal2(12), Dpal3(50))
cbind(unique_seqs, col_vector[1:n])->unique_seqs

#write this as the master color for each ASV sequence
if(genus=="Sphingomonas"){ 
	write.table(unique_seqs, file=paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/sphingo_16Sseqs_to_color_", date, ".txt", sep="", collapse=""), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
}
if(genus=="Pseudomonas"){ 
	write.table(unique_seqs, file=paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/pseudo_16Sseqs_to_color_", date, ".txt", sep="", collapse=""), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
}



#Match the rownames to an ASV number and rename the tips
unique_seqs[match(associated_seqs_of_rownames, unique_seqs[,1]),3]->ASVs_rownames


#ASVs_rownames->rownames(reads)



#Match the tip labels to an ASV color and rename the tips
unique_seqs[match(associated_seqs_of_rownames, unique_seqs[,1]),4]->colors_rownames
unique_seqs[match(associated_seqs_of_rownames, unique_seqs[,1]),3]->ASVnames
#anything that didn't match can be colored black
colors_rownames[is.na(colors_rownames)]="white"




#augment refseq genomes with Sphingomonas species names
if(genus=="Sphingomonas"){
	gsub("GCF_", "GCF", core$tip.label)->core$tip.label
	gsub("_S", "-S", core$tip.label)->core$tip.label   #so nanopore also not lost
	gsub("_.*", "", core$tip.label)->core$tip.label
	gsub("GCF", "GCF_", core$tip.label)->core$tip.label
	GCF_to_species<-read.table(file="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/REFSEQ_metadata/REFSEQ_GCF_to_species.txt", sep="\t")
	GCF_to_species[match(core$tip.label, GCF_to_species[,2], nomatch=0),5]->species
	GCF_to_species[match(core$tip.label, GCF_to_species[,2], nomatch=0),2]->OG
	paste(OG, species, sep="_")->extended_names
	extended_names->core$tip.label[match(core$tip.label, GCF_to_species[,2], nomatch=0)>0]
}

#make all font colored black
colors_edgelabels=rep("black", length(colors_rownames))

if(nanopore_control==TRUE){
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
}


#all(actualorder==correct_order)
	coredata=data.frame(tiplabels = core$tip.label, ASVnames=ASVnames, ASVcolors=colors_rownames)
	factor(coredata$ASVcolors, levels=unique(as.matrix(coredata$ASVcolors)))->coredata$ASVcolors
	
if(genus=="Sphingomonas"){
	#pdf(file=paste(basefolder, "strain_tree_colors.pdf", sep="", collapse=""), width = 10, height = 20, useDingbats=FALSE)
	##plot(core, cex=.2, tip.color=colors_edgelabels, align.tip.label = FALSE, show.tip.label=FALSE)
	#plot(core, cex=.2, tip.color=colors_rownames, align.tip.label = FALSE, show.tip.label=FALSE)
	#	tiplabels(pch=15, adj=c(.5, .5), col=colors_rownames, cex=0.4)
	#dev.off()
	#tree dataframe
	
	#search for sent strains in tree Vogel 27 Feb 2020
	#sent_strains=c("S8H113","S151H113","S214H113","S338H113","S132H113","S187H113","S328H113","S139H133","S337H113","S129H113","S135H113","S221H113","S252H113","S188H113","S123H113","S223H113","S133H113","S246H133","S227H113","S211H113","S191H113","S334H113","S122H113","S220H113","S230H133","S136H113","S259H113","S226H113","S98H133","S225H113","S339H113","S367H113","S219H113","S218H113","S127H113","S340H113","S217H113","S212H113","S231H113","S130H113","S222H113","S128H113","S213H113","S341H113","S224H113","S11H113","S139H113","S190H113","S330H113","S216H113","S215H113","S12H113","S228H113","S339H113","S18H113","S114H113","S230H113","S101H133","S111H133","S189H113","S380H113","S330H113","S126H113","S237H113")
	#sent_strains<-grep("S8H113|S151H113|S214H113|S338H113|S132H113|S187H113|S328H113|S139H133|S337H113|S129H113|S135H113|S221H113|S252H113|S188H113|S123H113|S223H113|S133H113|S246H133|S227H113|S211H113|S191H113|S334H113|S122H113|S220H113|S230H133|S136H113|S259H113|S226H113|S98H133|S225H113|S339H113|S367H113|S219H113|S218H113|S127H113|S340H113|S217H113|S212H113|S231H113|S130H113|S222H113|S128H113|S213H113|S341H113|S224H113|S11H113|S139H113|S190H113|S330H113|S216H113|S215H113|S12H113|S228H113|S339H113|S18H113|S114H113|S230H113|S101H133|S111H133|S189H113|S380H113|S330H113|S126H113|S237H113",core$tip.label)
	#sent_record=rep("", length=length(core$tip.label))
	#sent_record[sent_strains]<-"sent"
	#paste(core$tip.label, sent_record)->core$tip.label
		
	pdf(file=paste(basefolder, "strain_tree_colors_", date, ".pdf", sep="", collapse=""), width = 10, height = 20, useDingbats=FALSE)
	
	alphas=rep(1, length(coredata$ASVcolors))
	#alphas[grep("-S", core$tip.label)]<-1 #nanopores darker	
	similar_to_pseudomonas=c("S282H133","S128H113","S227H113","S220H113","S363H113","S301H113","S283H113","S249H113","S300H113","S226H113","S377H113")
	alphas[match(similar_to_pseudomonas, core$tip.label)]<-1  #similar to pseudomonas darker
	
	plot(core, cex=0.2)
	#p=ggtree(core) %<+% coredata 
	#p + geom_tippoint(color=coredata$ASVcolors, shape="square", size=0.7) +  xlim(0, 2) + ylim(-10, 500) +
	#	geom_tiplab(color=coredata$ASVcolors, size=0, align=TRUE, linetype='solid', linesize=1.2, alpha=alphas) + 
	#	geom_tiplab(color="black", size=.6, align=TRUE, linesize=NA) + 
	#	geom_treescale(x=0, y=-10, width=1, color='black') 
	
	dev.off()
	
	write.table(core$tip.label[length(core$tip.label):1], file=paste(basefolder, "strain_tree_order_", date, ".txt", sep="", collapse=""), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")

}
if(genus=="Pseudomonas"){
	pdf(file=paste(basefolder, "strain_tree_colors_", date, ".pdf", sep="", collapse=""), width = 10, height = .56*20, useDingbats=FALSE)
	
	alphas=rep(0.5, length(coredata$ASVcolors))
	alphas[grep("-S", core$tip.label)]<-1 #nanopores darker
	
	#plot(core, cex=.2, tip.color=colors_edgelabels, align.tip.label = FALSE, show.tip.label=FALSE)
	plot(core, cex=0.2)
	dev.off()
	
	write.table(core$tip.label[length(core$tip.label):1], file=paste(basefolder, "strain_tree_order_", date, ".txt", sep="", collapse=""), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
}		
		
	
library(gplots)

#1 Sept 2019
#pare down heatmap so that Athaliana genomes are reduced
if("pare_down"=="pare_deeeown"){
for(r in 1:3){
	drop_list=vector(length=0)
	as.matrix(metadata$host_species[match(core$tip.label,metadata$SeqID)])->actual_plant_tiplabels
	actual_plant_tiplabels[grep("Leaf", core$tip.label)]<-"Athaliana"
	actual_plant_tiplabels[is.na(actual_plant_tiplabels)]<-"unknown"
	for(i in 2:(length(actual_plant_tiplabels)-1)){
		if(actual_plant_tiplabels[i-1]=="Athaliana" & actual_plant_tiplabels[i]=="Athaliana" & actual_plant_tiplabels[i+1]=="Athaliana"){
			drop_list=c(drop_list, i)
		}
		if(actual_plant_tiplabels[i-1]=="unknown" & actual_plant_tiplabels[i]=="unknown" & actual_plant_tiplabels[i+1]=="unknown"){
			drop_list=c(drop_list, i)
		}
		if(actual_plant_tiplabels[i-1]=="Athaliana" & actual_plant_tiplabels[i]=="unknown" & actual_plant_tiplabels[i+1]=="Athaliana"){
			drop_list=c(drop_list, i)
		}
	}
	drop.tip(core, drop_list)->core
	coredata[setdiff(1:nrow(coredata), drop_list),]->coredata
   }
   actual_plant_tiplabels[actual_plant_tiplabels=="unknown"]<-NA
}
if("onlyH133"=="onlyH133"){
	if(genus=="Sphingomonas"){ drop_list=grep("H133", core$tip.label, invert=TRUE) }
	if(genus=="Pseudomonas"){ drop_list=grep("H133", core$tip.label, invert=TRUE) }
	drop.tip(core, drop_list)->core
	coredata[setdiff(1:nrow(coredata), drop_list),]->coredata
}

	
collection_color_list<-matrix(c(
	'spring', "#00FF23",
    'summer', "#D09400", 
    "TK", "#41C4E2", 
    "Roger", "#41C4E2", 
    'At-LSPHERE', "#39AF49",
    'Eyach', "#41C4E2",
    'REFSEQ', "#000000"
), ncol=2, byrow=TRUE)

as.matrix(coredata$ASVcolors)->colors_rownames #for some reason i changed to colors rownames... anyway...
#Put all heatmap colors together and assign and integer value to them
all_heatmap_colors=c(unique(plant_color_list[,2]), 
						unique(colors_rownames),
						unique(collection_color_list[,2]))
heatmap_color_key=cbind(all_heatmap_colors, 1:length(all_heatmap_colors))
colors=unique(colors_rownames)
black=which(heatmap_color_key[,1]=="#000000")[1]

#find color for the plant species of origin
as.matrix(metadata$host_species[match(gsub("Sph.*-", "", core$tip.label),metadata$SeqID)])->actual_plant_tiplabels
plant_color_list[match(actual_plant_tiplabels, plant_color_list[,1]),2]->actual_plant_colorlabels
as.numeric(heatmap_color_key[match(actual_plant_colorlabels, heatmap_color_key[,1]),2])->actual_plant_color_code


#find color for the OTU
as.numeric(heatmap_color_key[match(colors_rownames, heatmap_color_key[,1]),2])->actual_OTU_color_code
actual_OTU_color_code[which(actual_OTU_color_code==16)]<-NA   #turn black to NA

#find color for the collection
if(genus=="Sphingomonas"){
	gsub("S.*H133.*", "Eyach", core$tip.label)->collection_color
	gsub("S.*H113.*", "Eyach", collection_color)->collection_color
	gsub(".*Leaf.*", "REFSEQ", collection_color)->collection_color
	gsub(".*Root.*", "REFSEQ", collection_color)->collection_color
	gsub("GCF.*", "REFSEQ", collection_color)->collection_color
	gsub(".*GCF.*", "REFSEQ", collection_color)->collection_color
}
if(genus=="Pseudomonas"){
	gsub(".*H133.*", "Eyach", core$tip.label)->collection_color
	gsub(".*H113.*", "Eyach", collection_color)->collection_color
	gsub("^p.*", "Karasov", collection_color)->collection_color
	gsub(".*GCF.*", "REFSEQ", collection_color)->collection_color
}
#find color for the genome type
collection_color_list[match(collection_color, collection_color_list[,1]),2]->actual_collection_colorlabels
as.numeric(heatmap_color_key[match(actual_collection_colorlabels, heatmap_color_key[,1]),2])->actual_collection_color_code


as.matrix(metadata$plantID[match(gsub("Sph.*-", "", core$tip.label),metadata$SeqID)])->collection_color2
collection_color2[is.na(collection_color2)]<-""
gsub(".*Og.*", "spring", collection_color2)->collection_color2
gsub(".*S.*", "summer", collection_color2)->collection_color2
gsub(".*N.*", "Roger", collection_color2)->collection_color2
gsub(".*O.*", "Roger", collection_color2)->collection_color2
gsub("^102$", "Roger", collection_color2)->collection_color2
gsub("EY", "TK", collection_color2)->collection_color2
#find color for the collection time
collection_color_list[match(collection_color2, collection_color_list[,1]),2]->actual_time_colorlabels
as.numeric(heatmap_color_key[match(actual_time_colorlabels, heatmap_color_key[,1]),2])->actual_time_color_code


#add white and gray to heatmap color key as 25, 26 respectively
rbind(heatmap_color_key, c("#FFFFFF", 25), c("#A0A0A0", 26))->heatmap_color_key
heatmap_color_key[1:26,]->heatmap_color_key

#@@@@@@@
as.matrix(metadata$plantID[match(gsub("Sph.*-", "", core$tip.label),metadata$SeqID)])->collection_plant
collection_plant[is.na(collection_plant)]<-""
unique(collection_plant)->uniq_plants
cbind(uniq_plants, col_vector[1:length(uniq_plants)])->collection_plant_list
collection_plant_list[match(collection_plant, collection_plant_list[,1]),2]->actual_plant_colorlabels

matrix(nrow=length(actual_plant_colorlabels), ncol=length(sort(table(actual_plant_colorlabels))))->sample_ID_mat
sort(table(actual_plant_colorlabels), decreasing=TRUE)->all_plantIDs
for(cl in 1:ncol(sample_ID_mat)){
	if(cl%%2==0){sample_ID_mat[,cl]<-26}  #even coloumns grey
	if(cl%%2==1){sample_ID_mat[,cl]<-25}  #odd coloumns white
	sample_ID_mat[which(actual_plant_colorlabels==names(all_plantIDs[cl])), cl]<-black
}



image(1:length(actual_plant_colorlabels), 1, as.matrix(1:length(actual_plant_colorlabels)), 
     col=actual_plant_colorlabels, xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")


image(1:length(heatmap_color_key[,1]), 1, as.matrix(1:length(heatmap_color_key[,1])), 
     col=heatmap_color_key[,1], xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")

      


#change OTU information from colors into a different one each column
matrix(nrow=length(actual_OTU_color_code), ncol=length(sort(table(actual_OTU_color_code))))->OTU_ID_mat
sort(table(actual_OTU_color_code), decreasing=TRUE)->all_OTUs
for(cl in 1:ncol(OTU_ID_mat)){
	OTU_ID_mat[which(actual_OTU_color_code==names(all_OTUs[cl])), cl]<-as.numeric(names(all_OTUs[cl]))
}
	 
	 

#change plant species information so different species each column
matrix(nrow=length(actual_plant_color_code), ncol=length(sort(table(actual_plant_color_code))))->plant_ID_mat
sort(table(actual_plant_color_code), decreasing=TRUE)->all_plantIDs
actual_plant_tiplabels[match(names(all_plantIDs), actual_plant_color_code)]


for(cl in 1:ncol(plant_ID_mat)){
	place_to_put_color=( match(actual_plant_color_code, names(all_plantIDs[cl]))>0 )
	color_code_to_use=actual_plant_color_code[match(names(all_plantIDs[cl]), actual_plant_color_code)]
	plant_ID_mat[which(actual_plant_color_code==names(all_plantIDs[cl])), cl]<-color_code_to_use #as.numeric(names(all_plantIDs[cl]))
}


	
#change genome type information so different types each column
#distinguish ASV1 from the others

matrix(nrow=length(actual_collection_color_code), ncol=length(sort(table(actual_collection_color_code))))->collection_ID_mat
sort(table(actual_collection_color_code), decreasing=TRUE)->all_collectionIDs
for(cl in 1:ncol(collection_ID_mat)){
	collection_ID_mat[which(actual_collection_color_code==names(all_collectionIDs[cl])), cl]<-black #as.numeric(names(all_collectionIDs[cl]))
}

#change time information from colors into a different one each column
matrix(nrow=length(actual_time_color_code), ncol=length(sort(table(actual_time_color_code))))->time_ID_mat
sort(table(actual_time_color_code), decreasing=TRUE)->all_times
for(cl in 1:ncol(time_ID_mat)){
	time_ID_mat[which(actual_time_color_code==names(all_times[cl])), cl]<-black #as.numeric(names(all_times[cl]))
}
	 

REFSEQ_genomes=rep(NA, length=length(actual_OTU_color_code))
	black->REFSEQ_genomes[which(collection_color=="REFSEQ")]
	
NANOPORE_genomes=rep(NA, length=length(actual_OTU_color_code))
	black->NANOPORE_genomes[grep("S337H113|S132H113|S136H113|S380H113|S237H113|S230H113|S190H113|S213H113|S18H113|S127H113|S133H113|S216H113", core$tip.label)]
	
	
actual_OTU_color_code

	
#bind all the vectors into a heatmap
space_between_items=4
cbind(  actual_OTU_color_code,
		matrix(data=NA, nrow=length(actual_OTU_color_code), ncol=space_between_items),
		actual_plant_color_code,
		plant_ID_mat,
		REFSEQ_genomes,
		matrix(data=NA, nrow=length(actual_OTU_color_code), ncol=space_between_items),
		actual_collection_color_code,
		collection_ID_mat,
		matrix(data=NA, nrow=length(actual_OTU_color_code), ncol=space_between_items),
		actual_time_color_code,
		time_ID_mat,
		NANOPORE_genomes,
		matrix(data=NA, nrow=length(actual_OTU_color_code), ncol=space_between_items),
		sample_ID_mat
		)->mat
mat[nrow(mat):1,]->mat


if(genus=="Sphingomonas"){ pdf(file=paste(basefolder, "strain_tree_color_labels_", date, ".pdf", sep="", collapse=""), width = 5, height = 5, useDingbats=FALSE) }
if(genus=="Pseudomonas"){ pdf(file=paste(basefolder, "strain_tree_color_labels_", date, ".pdf", sep="", collapse=""), width = 5, height = 5, useDingbats=FALSE) }
heatmap.2(mat, breaks=c(0:nrow(heatmap_color_key)), col=heatmap_color_key[,1], 
	trace="none", Colv=FALSE, Rowv=FALSE, key=FALSE, dendrogram="none", 
	labRow = "",  labCol ="",
	sepwidth=c(0.01,0.01),
           sepcolor="black",
         colsep=0:ncol(mat))
dev.off()

if(genus=="Sphingomonas"){ pdf(file=paste(basefolder, "strain_tree_colors_mini_", date, ".pdf", sep="", collapse=""), width = 10, height = 20, useDingbats=FALSE) }
if(genus=="Pseudomonas"){ pdf(file=paste(basefolder, "strain_tree_colors_mini_", date, ".pdf", sep="", collapse=""), width = 10, height = 20, useDingbats=FALSE) }
	plot(core, cex=0.2)	
dev.off()



	 
	
	
