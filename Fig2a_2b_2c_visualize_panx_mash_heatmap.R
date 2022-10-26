install.packages("RColorBrewer")
install.packages("pheatmap")
install.packages("viridis")
install.packages("gplots")

library(RColorBrewer)
library(pheatmap)
library(viridis)
library(gplots)


source("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/scripts/R/functions/microbiome_custom_functions.R")


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
if(matrixType=="WGS"){ fullmatrix_Sph<-read.table(file="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/MASH/fullmatrix.txt", sep="\t") }
if(matrixType=="core"){ 
	fullmatrix_Sph<-read.table(file="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/coreGenome/coreGenome_MASHmat.txt", sep="\t")
	if(fake_core_pseudomonas==TRUE){ fullmatrix_Sph<-read.table(file="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/MASH/fullmatrix.txt", sep="\t") }
	WGS_MASH_order<-as.matrix(read.table(file="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/coreGenome/order_of_Sphing_in_WGS_heatmap.txt", sep="\t")	)
}

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

#plant heatmap  load a metadata table about which plant each genome came from
plantspecies=read.table(file="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/MASH/single_genomes_plantspecies.txt", sep="\t")

#augment refseq genomes with Sphingomonas species names
gsub("GCF_", "GCF", colnames(reads))->colnames(reads)
gsub("_.*", "", colnames(reads))->colnames(reads)
gsub("GCF", "GCF_", colnames(reads))->colnames(reads)
gsub("GCF_", "GCF", rownames(reads))->rownames(reads)
gsub("_.*", "", rownames(reads))->rownames(reads)
gsub("GCF", "GCF_", rownames(reads))->rownames(reads)
GCF_to_species<-read.table(file="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/REFSEQ_metadata/REFSEQ_GCF_to_species.txt", sep="\t")
GCF_to_species[match(rownames(reads), GCF_to_species[,2], nomatch=0),5]->species
GCF_to_species[match(rownames(reads), GCF_to_species[,2], nomatch=0),2]->OG
paste(OG, species, sep="_")->extended_names
extended_names->rownames(reads)[match(rownames(reads), GCF_to_species[,2], nomatch=0)>0]
GCF_to_species[match(colnames(reads), GCF_to_species[,2], nomatch=0),5]->species
GCF_to_species[match(colnames(reads), GCF_to_species[,2], nomatch=0),2]->OG
paste(OG, species, sep="_")->extended_names
extended_names->colnames(reads)[match(colnames(reads), GCF_to_species[,2], nomatch=0)>0]

#clean rownames
gsub(".fasta", "", rownames(reads))->rownames(reads)
colnames(reads)=rownames(reads)


#SELECT ONLY ATHALIANA STRAINS
if(matrixType=="WGS"){ 
	Athaliana_set1=which(as.matrix(plantspecies[match(rownames(reads), plantspecies[,1]),14])=="Athaliana")
	Athaliana_set2=grep("^p", rownames(reads))
	reads[c(Athaliana_set1, Athaliana_set2), c(Athaliana_set1, Athaliana_set2)]->reads
}


#appears to be trash, next few lines)
length(which(as.matrix(plantspecies[match(rownames(reads), plantspecies[,1]),4])=="Sphingomonas"))
Athaliana_set1=which(as.matrix(plantspecies[match(rownames(reads), plantspecies[,1]),14])=="Athaliana")
Athaliana_set2=which(as.matrix(plantspecies[match(rownames(reads), plantspecies[,1]),4])=="Sphingomonas")
length(intersect(Athaliana_set1, Athaliana_set2))

length(which(as.matrix(plantspecies[match(rownames(reads), plantspecies[,1]),4])=="Pseudomonas"))
Athaliana_set1=which(as.matrix(plantspecies[match(rownames(reads), plantspecies[,1]),14])!="Athaliana")


#  awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" }END { printf "%s", n }' OTUrepSeqs97.fa > OTUrepSeqs97.fixed.fa
#length(grep("Sphingomonas", plantspecies[match(rownames(reads2), plantspecies[,1]),4]))

Athaliana_set2=which(as.matrix(plantspecies[match(rownames(reads), plantspecies[,1]),4])=="Pseudomonas")
length(intersect(Athaliana_set1, Athaliana_set2))




column_scaling=ncol(reads)/nrow(reads)

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
mat_breaks <- quantile_breaks(reads, n = 11)


#reduce columns from refseq ID set until all columns are different enough
remove_similar_GCF=FALSE
if(remove_similar_GCF==TRUE){
GCF_threshold=0.2
for(g in rownames(reads)){
	if(length(reads[g,which(reads[g,]<=GCF_threshold)]>1)){
		reads[,c(which(reads[g,]>GCF_threshold), sample(which(reads[g,]==0),1))]->reads   #take all nonmatching columsn and one randomly chosen matching column
		print(dim(reads))
		#, invert=TRUE, value=TRUE),]->reads
	}
}
for(g in rownames(reads)){
	if(length(reads[g,which(reads[g,]<=GCF_threshold)]>1)){
		reads[,c(which(reads[g,]>GCF_threshold), sample(which(reads[g,]==0),1))]->reads   #take all nonmatching columsn and one randomly chosen matching column
		print(dim(reads))
		#, invert=TRUE, value=TRUE),]->reads
	}
}
}


#reads[reads>.05]<-1
#reads[reads<=.05]<-0



#COLOR GENOME NAME BY OTU IDENTITY
#LOAD CULTRE FASTA AND EXTRACT SPHINGO SEQS
culture_fasta=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/KarasovAlmario16S/all_cultured_16S_20190410_reoriented_trimmed.fa")
as.matrix(culture_fasta)->culture_fasta
culture_fasta[sort(c(grep("Sphingo", culture_fasta), (grep("Sphingo", culture_fasta)+1))),]->Sphingo16S_culture
	split_FASTA(Sphingo16S_culture, add_integer_column=FALSE, leave_carrot=FALSE, reverse=FALSE)->Sphingo16S_table
	gsub("H133.Sphingo_16S", "H133" , Sphingo16S_table[,1])->Sphingo16S_table[,1]
	gsub("H113.Sphingo_16S", "H113" , Sphingo16S_table[,1])->Sphingo16S_table[,1]
	gsub("Sphingo", "" , Sphingo16S_table[,1])->Sphingo16S_table[,1]
	gsub("_16S", "" , Sphingo16S_table[,1])->Sphingo16S_table[,1]
	if(otu99==TRUE){
		culture_fasta=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/KarasovAlmario16S/cluster_cultured_to_99pcnt/all_cultured_16S_20190818_reoriented_trimmed_FL_995otus.fa")
		as.matrix(culture_fasta)->culture_fasta
		split_FASTA(culture_fasta, add_integer_column=FALSE, leave_carrot=FALSE, reverse=FALSE)->Sphingo16S_table
	}	

	#gyrB info
		gyrB_fasta=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/Fabrice/all_gyrB_seqs.txt")
		as.matrix(gyrB_fasta)->gyrB_fasta
		split_FASTA(gyrB_fasta, add_integer_column=FALSE, leave_carrot=FALSE, reverse=FALSE)->Sphingo_gyrB_table
		Sphingo_gyrB_table->Pseudo_gyrB_table



		culture_fasta[c(sort(c(grep("Pseudo", culture_fasta), (grep("Pseudo", culture_fasta)+1))),
				sort(c(grep(">p[0-9]", culture_fasta), (grep(">p[0-9]", culture_fasta)+1)))),]->Pseudo16S_culture
								
		split_FASTA(Pseudo16S_culture, add_integer_column=FALSE, leave_carrot=FALSE, reverse=FALSE)->Pseudo16S_table
		gsub("H133.Pseudo_16S", "H133" , Pseudo16S_table[,1])->Pseudo16S_table[,1]
		gsub("H113.Pseudo_16S", "H113" , Pseudo16S_table[,1])->Pseudo16S_table[,1]
		gsub("Pseudo", "p" , Pseudo16S_table[,1])->Pseudo16S_table[,1]
		gsub("_16S", "" , Pseudo16S_table[,1])->Pseudo16S_table[,1]



rbind(Sphingo16S_table, Pseudo16S_table)->full16S_table
#Turn tip labels into the actual V3V4 ASV sequences
full16S_table[match(rownames(reads), full16S_table[,1]),2]->associated_seqs_of_rownames
#Rename cultured seqs into ASVs
as.matrix(table(associated_seqs_of_rownames))->unique_seqs
unique_seqs[order(unique_seqs[,1], decreasing=TRUE),]->unique_seqs
cbind(names(unique_seqs), unique_seqs)->unique_seqs
rownames(unique_seqs)<-NULL
unique_seqs=cbind(unique_seqs, paste("ASV", 1:nrow(unique_seqs), sep=""))


library(RColorBrewer)
n <- nrow(unique_seqs)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
cbind(unique_seqs, col_vector[1:n])->unique_seqs

#Match the rownames to an ASV number and rename the tips
unique_seqs[match(associated_seqs_of_rownames, unique_seqs[,1]),3]->ASVs_rownames
#ASVs_rownames->rownames(reads)
#Match the tip labels to an ASV color and rename the tips
unique_seqs[match(associated_seqs_of_rownames, unique_seqs[,1]),4]->colors_rownames
#anything that didn't match can be colored black
colors_rownames[is.na(colors_rownames)]="black"

SphingoASV1="TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTTACCCGGGATGATAATGACAGTACCGGGAGAATAAGCTCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGAGCTAGCGTTATTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTAAGAGGTGAAAGCCCAGAGCTCAACTCTGGAATTGCCTTTTAGACTGCATCGCTTGAATCATGGAGAGGTCAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGGCTGACTGGACATGTATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG"
if(otu99==TRUE){
	SphingoASV1="AGAGTTTGATCATGGCTCAGAATGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGAATTACCTTCGGGTAGTTAGTGGCGCACGGGTGCGTAACGCGTGGGAATCTGCCCCTTGGTTCGGAATAACAGTTGGAAACGACTGCTAATACCGGATGATGACGTAAGTCCAAAGATTTATCGCCGAGGGATGAGCCCGCGTAGGATTAGGTAGTTGGTGTGGTAAAGGCGCACCAAGCCGACGATCCTTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTTACCCGGGATGATAATGACAGTACCGGGAGAATAAGCTCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGAGCTAGCGTTATTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTAAGAGGTGAAAGCCCAGAGCTCAACTCTGGAATTGCCTTTTAGACTGCATCGCTTGAATCATGGAGAGGTCAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGGCTGACTGGACATGTATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGATAACTAGCTGTCCGGACACTTGGTGTTTGGGTGGCGCAGCTAACGCATTAAGTTATCCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCTGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGCAGAACCTTACCAGCGTTTGACATGGCAGGACGACTTCCAGAGATGGATTTCTTCCCTTCGGGGACCTGCACACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTCGCCTTTAGTTGCCATCATTTAGTTGGGCACTTTAAAGGAACCGCCGGTGATAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTACGCGCTGGGCTACACACGTGCTACAATGGCGGTGACAGTGGGCAGCAAGCACGCGAGTGTGAGCTAATCTCCAAAAGCCGTCTCAGTTCGGATTGTTCTCTGCAACTCGAGAGCATGAAGGCGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCAGGCCTTGTACACACCGCCCGTCACACCATGGGAGTTGGATTCACCCGAAGGCGTTGCGCTAACTCAGCAATGAGAGGCAGGCGACCACGGTGGGTTTAGCGACTGGGGTGAAGTCGTAACAAGGT"
}
PseudoASV1="TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAGCACTTTAAGTTGGGAGGAAGGGCAGTAACCTAATACGTTATTGTTTTGACGTTACCGACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGAATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCCAAAACTGGCAAGCTAGAGTAGGGCAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGGCTCATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACAGG"


pal1 <- colorRampPalette(c("black","#3f324f","#6A1B9A","#E85285", "#FFECB3"))(n=300)
pal1 <- inferno(100)
#pal1=pal1[length(pal1):1] #invert
colnames(reads)=rownames(reads)


#TO SELECT ONLY ASV1 sphingo AND ASV1 pseudo
ASV1sph=grep((SphingoASV1), associated_seqs_of_rownames)
ASV1pse=grep((PseudoASV1), associated_seqs_of_rownames)
#there are more sphingomonas than there are diverse OTU5 strains... so subsample sphingomonas to match. Shouldn't
	#need to do this... cause talia has so many OTU5...!  FIX
if(matrixType=="WGS"){
	low_quality_genomes=c("S165H113","S67H113","S167H113","S169H113","S170H113","S171H113","S168H113","S74H113","S172H113","S372H113")
		length(ASV1sph)
	ASV1sph=ASV1sph[match(rownames(reads)[ASV1sph], low_quality_genomes,  nomatch=0)==0]
		length(ASV1sph)
	#sample(ASV1sph, length(ASV1pse))->ASV1sph
	#for reproducibility of clustering for figure making, keep a single subsample effort
	#ASV1sph=c("126","80","188","83","222","191","176","178","229","322","366","298","237","116","148","272","367","77","177","186","201","161","252","183","223","140","130","11","139","172","163","248","174","185","194","64","207","203","292","142","324","369","236","357","187","334","338","71","158","8","197","273","170","361","347","284","141","267","85","166","5","213","196","79","243","156","25","152","251","348","208","214","235","225","61","233","12","362","204","364","82","133","327","165","28","119","231","226","180","263","68","153","179","199","280","9","265","60","209")
 	#as.numeric(ASV1sph)->ASV1sph 	
}

rep("black", length(colors_rownames))->colors_rownames

if(matrixType=="core"){ 
	reads2=reads
	keep=match(WGS_MASH_order, rownames(reads2), nomatch=0)
	reads2[keep, keep]->reads2
	dendrogram="none"
	Colv=Rowv=FALSE
}
if(matrixType=="WGS"){ 
	dendrogram="both"
	Colv=Rowv=TRUE
	ASV1sphpse=c(ASV1sph, ASV1pse)
	associated_seqs_of_rownames[ASV1sphpse]->associated_seqs_of_rownames
	reads[ASV1sphpse,ASV1sphpse]->reads2
	colors_rownames[ASV1sphpse]->colors_rownames
	ASVs_rownames[ASV1sphpse]->ASVs_rownames
}
	
100*(1-reads2)->reads2

dim(reads2)

#
	#
		#
			#  Calculate colors for OTUs based on 99% matches and gyrB matches if doing core genome 
		#   #to generate the OTU map for the side of the heatmap, need to do the heatmap first, which reorders "reads2". Then can do this.
	#		
#
if(matrixType=="core" & otu99==FALSE){  #needs to be false, because if true, means only ASV1 considering full length OTUs is selected and others not kept
	#use unique seqs as calculated on all reads, not the condensed ones
	#Drop=which(match(core$tip.label,c("S337H113","S132H113","S136H113","S380H113","S237H113","S230H113","S190H113","S213H113","S18H113","S127H113","S133H113","S216H113"), nomatch=0)>0)
	#drop.tip(core, as.numeric(Drop))->core
#SPHINGOMONAS
	#FL 16S
	culture_fasta=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/KarasovAlmario16S/cluster_cultured_to_99pcnt/all_cultured_16S_20190818_reoriented_trimmed_FL_995otus.fa")
	as.matrix(culture_fasta)->culture_fasta
		culture_fasta[c(
			sort(c(grep("S.*H.*", culture_fasta), (grep("S.*H.*", culture_fasta)+1))),
			sort(c(grep("GCF_", culture_fasta), (grep("GCF_", culture_fasta)+1)))),]->Sphingo16S_culture
			split_FASTA(Sphingo16S_culture, add_integer_column=FALSE, leave_carrot=FALSE, reverse=FALSE)->Sphingo16S_table
			Sphingo16S_table[grep(SphingoASV1, Sphingo16S_table[,2]),]->Sphingo16S_table  #select SPHINGOMONAS again by seq
	Sphingo16S_table[match(rownames(reads), Sphingo16S_table[,1]),2]->associated_seqs_of_rownames
	Sphingo16S_table[match(rownames(reads2), Sphingo16S_table[,1]),2]->associated_seqs_of_rownames_reads2  #this is here because these will be subselected for the color map.
	
	#gyrB... put also the seqs in the same order as the rownames
	Sphingo_gyrB_table[match(rownames(reads2), Sphingo_gyrB_table[,1]),2]->associated_gyrB_of_rownames_reads


#Rename cultured seqs into ASVs
	as.matrix(table(associated_seqs_of_rownames_reads2))->unique_seqs   #*
	unique_seqs[order(unique_seqs[,1], decreasing=TRUE),]->unique_seqs
	cbind(names(unique_seqs), unique_seqs)->unique_seqs
	rownames(unique_seqs)<-NULL
	unique_seqs=cbind(unique_seqs, paste("ASV", 1:nrow(unique_seqs), sep=""))
	#same for gyrB = rename cultured seqs into ASVs
		as.matrix(table(associated_gyrB_of_rownames_reads))->unique_gyrB_seqs   #*
		unique_gyrB_seqs[order(unique_gyrB_seqs[,1], decreasing=TRUE),]->unique_gyrB_seqs
		cbind(names(unique_gyrB_seqs), unique_gyrB_seqs)->unique_gyrB_seqs
		rownames(unique_gyrB_seqs)<-NULL
		unique_gyrB_seqs=cbind(unique_gyrB_seqs, paste("ASV", 1:nrow(unique_gyrB_seqs), sep=""))
		
		
library(RColorBrewer)
	n <- nrow(unique_seqs)
	qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
	col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
	col_vector[25:length(col_vector)]->col_vector
	"yellow"->col_vector[11]
	col_vector=c("#ED6925FF","#460B5DFF","#ffff20","#FCA209FF","#67CC5CFF","#FB9606FF","#297B8EFF","#A22B62FF","#472E7CFF","#F67E14FF","#404588FF","#375B8DFF","#F98A0BFF","#2E6F8EFF","#481467FF","#BADE28FF","#6E186EFF","#C93F4BFF","#24878EFF","#4F0D6CFF","#21A585FF","#482576FF","#3D4D8AFF","#433E85FF","#34618DFF","#3A538BFF","#26818EFF","#87D549FF","#781C6DFF","#2D0B5AFF","#83206BFF","#97D83FFF","#440A68FF","#59106EFF","#59C864FF","#F4E258FF","#C03A51FF","#FCAE12FF","#0F092CFF","#76D153FF","#390963FF","#DDE318FF","#64156EFF","#228D8DFF","#AD305DFF","#A9DB33FF","#EDE51BFF","#F6D543FF","#1FA088FF","#CBE11EFF","#35B779FF","#440154FF","#1F998AFF","#FDE725FF","#E15635FF","#2DB17EFF","#481D6FFF","#453581FF","#08051DFF","#F3F78CFF","#25AC82FF","#F1ED71FF","#982766FF","#31688EFF","#8D2369FF","#230C4BFF","#D24644FF","#FBBC21FF","#180C3CFF","#40BC72FF","#F9C830FF","#F1741CFF","#D94D3DFF","#2B748EFF","#000004FF","#02020FFF","#E75F2EFF","#20938CFF","#FCFFA4FF","#4CC26CFF")

#OTU99 colors
	cbind(unique_seqs, col_vector[1:n])->unique_seqs
			#image(1:length(col_vector), 1, as.matrix(1:length(col_vector)), 
   	   #col=col_vector, xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")

#Match the rownames to an ASV number and rename the tips
	unique_seqs[match(associated_seqs_of_rownames, unique_seqs[,1]),3]->ASVs_rownames
#Match the tip labels to an ASV color and rename the tips. Here we use the reads2 variant
	unique_seqs[match(associated_seqs_of_rownames_reads2, unique_seqs[,1]),4]->colors_rownames
	unique_seqs[match(associated_seqs_of_rownames_reads2, unique_seqs[,1]),3]->ASVnames
	#anything that didn't match can be colored black
	colors_rownames[is.na(colors_rownames)]="white"
	#image(1:nrow(reads2), 1, as.matrix(1:nrow(reads2)), 
	#      col=rep("#A7CDE0", nrow(reads2), xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n"))

	pdf(file=paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/coreGenome/otu99_colors_", date, ".pdf", sep="", collapse=""), width = 3, height=2, useDingbats=FALSE)
		image(1:length(colors_rownames), 1, as.matrix(1:length(colors_rownames)), 
   	   col=colors_rownames, xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
	dev.off()
	
	#gyrB colors
		cbind(unique_seqs, col_vector[1:n])->unique_seqs
			#image(1:length(col_vector), 1, as.matrix(1:length(col_vector)), 
   	   #col=col_vector, xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")

	#Match the rownames to an ASV number and rename the tips
	unique_seqs[match(associated_seqs_of_rownames, unique_seqs[,1]),3]->ASVs_rownames
	#Match the tip labels to an ASV color and rename the tips. Here we use the reads2 variant
	unique_seqs[match(associated_seqs_of_rownames_reads2, unique_seqs[,1]),4]->colors_rownames
	unique_seqs[match(associated_seqs_of_rownames_reads2, unique_seqs[,1]),3]->ASVnames
	#anything that didn't match can be colored black
	colors_rownames[is.na(colors_rownames)]="white"
	#image(1:nrow(reads2), 1, as.matrix(1:nrow(reads2)), 
	#      col=rep("#A7CDE0", nrow(reads2), xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n"))

	pdf(file=paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/coreGenome/otu99_colors_", date, ".pdf", sep="", collapse=""), width = 3, height=2, useDingbats=FALSE)
		image(1:length(colors_rownames), 1, as.matrix(1:length(colors_rownames)), 
   	   col=colors_rownames, xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
	dev.off()
	
	
#PSEUDOMONAS
	culture_fasta=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/KarasovAlmario16S/cluster_cultured_to_99pcnt/all_cultured_16S_20190818_reoriented_trimmed_FL_995otus.fa")
	as.matrix(culture_fasta)->culture_fasta
	headers=which(match(gsub(">", "", culture_fasta), rownames(reads2), nomatch=0)>0)
	headers_and_seqs=sort(c(headers, headers+1))
	culture_fasta[headers_and_seqs,]->Pseudo16S_culture
	split_FASTA(Pseudo16S_culture, add_integer_column=FALSE, leave_carrot=FALSE, reverse=FALSE)->Pseudo16S_table
	Pseudo16S_table[grep(PseudoASV1, Pseudo16S_table[,2]),]->Pseudo16S_table  #select pseudomonas again by seq
	Pseudo16S_table[match(rownames(reads), Pseudo16S_table[,1]),2]->associated_seqs_of_rownames
	Pseudo16S_table[match(rownames(reads2), Pseudo16S_table[,1]),2]->associated_seqs_of_rownames_reads2  #this is here because these will be subselected for the color map.
#Rename cultured seqs into ASVs
	as.matrix(table(associated_seqs_of_rownames))->unique_seqs
	unique_seqs[order(unique_seqs[,1], decreasing=TRUE),]->unique_seqs
	cbind(names(unique_seqs), unique_seqs)->unique_seqs
	rownames(unique_seqs)<-NULL
	unique_seqs=cbind(unique_seqs, paste("ASV", 1:nrow(unique_seqs), sep=""))
library(RColorBrewer)
	n <- nrow(unique_seqs)
	qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
	col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
	col_vector[25:length(col_vector)]->col_vector
	"yellow"->col_vector[11]
	col_vector=c("#460B5DFF","#ED6925FF","#B63458FF","#FCA209FF","#67CC5CFF","#FB9606FF","#297B8EFF","#A22B62FF","#472E7CFF","#F67E14FF","#404588FF","#375B8DFF","#F98A0BFF","#2E6F8EFF","#481467FF","#BADE28FF","#6E186EFF","#C93F4BFF","#24878EFF","#4F0D6CFF","#21A585FF","#482576FF","#3D4D8AFF","#433E85FF","#34618DFF","#3A538BFF","#26818EFF","#87D549FF","#781C6DFF","#2D0B5AFF","#83206BFF","#97D83FFF","#440A68FF","#59106EFF","#59C864FF","#F4E258FF","#C03A51FF","#FCAE12FF","#0F092CFF","#76D153FF","#390963FF","#DDE318FF","#64156EFF","#228D8DFF","#AD305DFF","#A9DB33FF","#EDE51BFF","#F6D543FF","#1FA088FF","#CBE11EFF","#35B779FF","#440154FF","#1F998AFF","#FDE725FF","#E15635FF","#2DB17EFF","#481D6FFF","#453581FF","#08051DFF","#F3F78CFF","#25AC82FF","#F1ED71FF","#982766FF","#31688EFF","#8D2369FF","#230C4BFF","#D24644FF","#FBBC21FF","#180C3CFF","#40BC72FF","#F9C830FF","#F1741CFF","#D94D3DFF","#2B748EFF","#000004FF","#02020FFF","#E75F2EFF","#20938CFF","#FCFFA4FF","#4CC26CFF")
	cbind(unique_seqs, col_vector[30:(29+n)])->unique_seqs
#Match the rownames to an ASV number and rename the tips
	unique_seqs[match(associated_seqs_of_rownames, unique_seqs[,1]),3]->ASVs_rownames
#Match the tip labels to an ASV color and rename the tips. Here we use the reads2 variant
	unique_seqs[match(associated_seqs_of_rownames_reads2, unique_seqs[,1]),4]->colors_rownames
	unique_seqs[match(associated_seqs_of_rownames_reads2, unique_seqs[,1]),3]->ASVnames
	#anything that didn't match can be colored black
	colors_rownames[is.na(colors_rownames)]="white"
	#image(1:nrow(reads2), 1, as.matrix(1:nrow(reads2)), 
	#      col=rep("#A7CDE0", nrow(reads2), xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n"))

	pdf(file=paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/coreGenome/otu99_colors_PSEUDO_", date, ".pdf", sep="", collapse=""), width = 3, height=2, useDingbats=FALSE)
		image(1:length(colors_rownames), 1, as.matrix(1:length(colors_rownames)), 
   	   col=colors_rownames, xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
	dev.off()	

}

#
	#
		#
			#   END COLOR SECTION
		#
	#
#



if(matrixType=="WGS"){ pdf(file=paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/MASH/compare_ASV1_Pseudo_to_ASV1_Sphingo_", date, ".pdf", sep="", collapse=""), width = 10*column_scaling, height = 10, useDingbats=FALSE) }
if(matrixType=="core"){
	if(otu99!=TRUE){
		pdf(file=paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/coreGenome/compare_ASV1_Sphingo_", date, ".pdf", sep="", collapse=""), width = 10*column_scaling, height = 10, useDingbats=FALSE)
		write.table(reads2, "/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/coreGenome/reads2_ASV1_V3V4.txt")
	}
	if(otu99==TRUE){
		write.table(reads2, "/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/coreGenome/reads2_ASV1FullLength99.txt")
	}
}

colors_rownames[which(colors_rownames=="white")]="black"
ab<-heatmap.2(reads2, density.info = "none",
                    trace = "none", dendrogram=dendrogram, colRow = colors_rownames,
                    colCol = colors_rownames, Colv=Colv, Rowv=Rowv, 
                    col = pal1, main = "", 
                    margins = c(9,9), cexRow = 0.1, cexCol = 0.1, 
                    keysize = 0.7,breaks = seq(50, 100, by=0.5),
                    key.title = NA, key.xlab = NA, 
                    sepwidth=c(0,0),
                    key.par=list(mar=c(1.5,0,0,5))   )
dev.off()
if(matrixType=="WGS"){
	write.table(rownames(reads2)[ab$rowInd], paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/coreGenome/order_of_Sphing_in_WGS_heatmap", date, ".txt", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE)
	#reads2=reads2[rownames(reads2)[ab$rowInd], rownames(reads2)[ab$rowInd]]
}




#
	#
		#
			#
#make the distance histogram to go underneath. Can only be done if matrixType="WGS"
reads2[grep(SphingoASV1, associated_seqs_of_rownames), grep(SphingoASV1, associated_seqs_of_rownames)]->ASV1sphREADS
	
	culture_fasta=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/KarasovAlmario16S/cluster_cultured_to_99pcnt/all_cultured_16S_20190818_reoriented_trimmed_FL_995otus.fa")
	as.matrix(culture_fasta)->culture_fasta
	culture_fasta[c(
			sort(c(grep("S.*H.*", culture_fasta), (grep("S.*H.*", culture_fasta)+1))),
			sort(c(grep("GCF_", culture_fasta), (grep("GCF_", culture_fasta)+1)))),]->Sphingo16S_culture
			split_FASTA(Sphingo16S_culture, add_integer_column=FALSE, leave_carrot=FALSE, reverse=FALSE)->Sphingo16S_table
			Sphingo16S_table[grep(SphingoASV1, Sphingo16S_table[,2]),]->Sphingo16S_table  #select SPHINGOMONAS again by seq
			Sphingo16S_table[match(rownames(reads2), Sphingo16S_table[,1]),2]->associated_seqs_of_rownames_reads2  
			SphingoASV1_FL="AGAGTTTGATCATGGCTCAGAATGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGAATTACCTTCGGGTAGTTAGTGGCGCACGGGTGCGTAACGCGTGGGAATCTGCCCCTTGGTTCGGAATAACAGTTGGAAACGACTGCTAATACCGGATGATGACGTAAGTCCAAAGATTTATCGCCGAGGGATGAGCCCGCGTAGGATTAGGTAGTTGGTGTGGTAAAGGCGCACCAAGCCGACGATCCTTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTTACCCGGGATGATAATGACAGTACCGGGAGAATAAGCTCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGAGCTAGCGTTATTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTAAGAGGTGAAAGCCCAGAGCTCAACTCTGGAATTGCCTTTTAGACTGCATCGCTTGAATCATGGAGAGGTCAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGGCTGACTGGACATGTATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGATAACTAGCTGTCCGGACACTTGGTGTTTGGGTGGCGCAGCTAACGCATTAAGTTATCCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCTGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGCAGAACCTTACCAGCGTTTGACATGGCAGGACGACTTCCAGAGATGGATTTCTTCCCTTCGGGGACCTGCACACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTCGCCTTTAGTTGCCATCATTTAGTTGGGCACTTTAAAGGAACCGCCGGTGATAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTACGCGCTGGGCTACACACGTGCTACAATGGCGGTGACAGTGGGCAGCAAGCACGCGAGTGTGAGCTAATCTCCAAAAGCCGTCTCAGTTCGGATTGTTCTCTGCAACTCGAGAGCATGAAGGCGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCAGGCCTTGTACACACCGCCCGTCACACCATGGGAGTTGGATTCACCCGAAGGCGTTGCGCTAACTCAGCAATGAGAGGCAGGCGACCACGGTGGGTTTAGCGACTGGGGTGAAGTCGTAACAAGGT"
reads2[grep(SphingoASV1_FL, associated_seqs_of_rownames_reads2), grep(SphingoASV1_FL, associated_seqs_of_rownames_reads2)]->ASV1sphREADS_FL
reads2[grep(PseudoASV1, associated_seqs_of_rownames), grep(PseudoASV1, associated_seqs_of_rownames)]->ASV1pseREADS

#get sphingo reads based on core genomes and based on nanopore genomes. Don't need to grep because it's already selected as ASV1
core_ASV1sphREADS=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/coreGenome/reads2_ASV1_V3V4.txt")
  core_nanopore=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_only_nanopore/Sanimat.txt")
		core_ASV1_fulllength=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/coreGenome/reads2_ASV1FullLength99.txt")
as.matrix(core_ASV1sphREADS)->core_ASV1sphREADS
  as.matrix(core_nanopore)->core_nanopore
  	as.matrix(core_ASV1_fulllength)->core_ASV1_fulllength

nrow(ASV1pseREADS)->Nps
nrow(ASV1sphREADS)->NspALL
nrow(core_ASV1sphREADS)->NspASV1_v3v4
nrow(core_nanopore)->NspASV1_nanopore
nrow(core_ASV1_fulllength)->NspASV1_FL


ASV1pseREADS[as.numeric(ASV1pseREADS)<100]->ASV1pseREADScount
	ASV1sphREADS[as.numeric(ASV1sphREADS)<100]->ASV1sphREADScount
		ASV1sphREADS_FL[as.numeric(ASV1sphREADS_FL)<100]->ASV1sphREADScount_FL

			
			core_ASV1sphREADS[as.numeric(core_ASV1sphREADS)<100]->core_ASV1sphREADScount
				core_nanopore[as.numeric(core_nanopore)<100]->core_nanopore
					core_ASV1_fulllength[as.numeric(core_ASV1_fulllength)<100]->core_ASV1_fulllength
		


mean(ASV1sphREADScount)
mean(ASV1pseREADScount)
mean(core_ASV1sphREADScount)
mean(core_nanopore)
mean(core_ASV1_fulllength)

library(ggplot2)
as.data.frame(rbind(
				cbind(rep("ASV1sph", length(ASV1sphREADScount)), ASV1sphREADScount),
				cbind(rep("ASV1pse", length(ASV1pseREADScount)), ASV1pseREADScount),
				cbind(rep("ASV1sFL", length(ASV1sphREADScount_FL)), ASV1sphREADScount_FL))
				)->histWGS
				as.numeric(as.character(histWGS[,2]))->histWGS[,2]
				histWGS[,1]=factor(histWGS[,1], levels=c("ASV1pse", "ASV1sph", "ASV1sFL"))
as.data.frame(rbind(
				cbind(rep("ASV1sCG", length(core_ASV1sphREADScount)), core_ASV1sphREADScount),
				cbind(rep("CGnanop", length(core_nanopore)), core_nanopore),
				cbind(rep("CGFL16S", length(core_ASV1_fulllength)), core_ASV1_fulllength))
				)->histCG
				as.numeric(as.character(histCG[,2]))->histCG[,2]
				histCG[,1]=factor(histCG[,1], levels=c("ASV1sCG", "CGnanop", "CGFL16S"))

			
c("group", "similarity")->names(histCG)

#de2e26	#000000	#000000
DerekColors=c("#77C043","#B71E57","#461E5C", "#F9A21A", "#461E5C", "#000000")
#WG MASH
pdf(file=paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/MASH/compare_ASV1_Pseudo_to_ASV1_Sphingo_DENSITY_WGS_", date, ".pdf", sep="", collapse=""), width = 5, height = 2.5, useDingbats=FALSE)
ggplot(histWGS, aes(x=histWGS[,2], fill=histWGS[,1], color=histWGS[,1])) +  
	geom_histogram(binwidth=.1) + 
	geom_density(alpha=.2) + theme_classic() + 
	theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15)) + 
	theme(axis.line = element_line(colour = "black", size = .15), axis.ticks = element_line(colour = "black", size = .15))+
	scale_color_manual(values=DerekColors) + scale_fill_manual(values=DerekColors) + 
	#scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), limits = c(0,1.2)) + 
	#scale_x_continuous(limits = c(80,100)) + 
	theme(axis.text.x = element_text(size=15, color="black"),axis.text.y = element_text(size=15, color="black")) + 
	theme(axis.line = element_line(color = "black", size = .15), axis.ticks = element_line(color = "black", size = .15))
dev.off()












DerekColors=c("#B71E57","#000000", "#461E5C")

pdf(file=paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/MASH/compare_ASV1_Pseudo_to_ASV1_Sphingo_DENSITY_CORE_", date, ".pdf", sep="", collapse=""), width = 5, height = 2.5, useDingbats=FALSE)
ggplot(histCG, aes(x=histCG[,2], fill=histCG[,1], color=histCG[,1])) +  
	geom_density(alpha=.2) + theme_classic() + 
	theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15)) + 
	theme(axis.line = element_line(colour = "black", size = .15), axis.ticks = element_line(colour = "black", size = .15))+
	scale_color_manual(values=DerekColors) + scale_fill_manual(values=DerekColors) + 
	scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), limits = c(0,1)) + 
	scale_x_continuous(limits = c(80,100)) + 
	theme(axis.text.x = element_text(size=15, color="black"),axis.text.y = element_text(size=15, color="black")) + 
	theme(axis.line = element_line(color = "black", size = .15), axis.ticks = element_line(color = "black", size = .15))
dev.off()
