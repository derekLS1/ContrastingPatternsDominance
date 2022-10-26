#visualize core genome tree in APE

scp dlundberg@burrito:/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/PanX/pan-genome-analysis/data/phylogeny/RAxML_bipartitionsBranchLabels.test /Users/dlundberg/Documents/abt6/Pratchaya/PanX/RAxML_bipartitionsBranchLabels.test     
scp dlundberg@burrito:/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/PanX/pan-genome-analysis/data/phylogeny/RAxML_bipartitions.test /Users/dlundberg/Documents/abt6/Pratchaya/PanX/RAxML_bipartitions.test     


library(ape)
install.packages("phangorn")
install.packages("phytools")
install.packages("geiger")
source("/Users/dlundberg/Documents/abt6/scripts/R/functions/microbiome_custom_functions.R")


core=read.tree(file = "/Users/dlundberg/Documents/abt6/Pratchaya/PanX/RAxML_bipartitionsBranchLabels.test")
	reroot(core, node.number=289)->core
	
meta=as.matrix(read.table("/Users/dlundberg/Documents/abt6/Pratchaya/PanX/ASV_tree_metadata.txt"))
plantmeta=as.matrix(read.table("/Users/dlundberg/Documents/abt6/Pratchaya/PanX/H133_plant_metadata.txt"))

culture_fasta=read.table("/Users/dlundberg/Documents/abt6/Pratchaya/KarasovAlmario16S/all_cultured_16S_20190410_reoriented_trimmed.fa")
as.matrix(culture_fasta)->culture_fasta
culture_fasta[sort(c(grep("Pseudo", culture_fasta), (grep("Pseudo", culture_fasta)+1))),]->Pseudo16S_culture
	split_FASTA(Pseudo16S_culture, add_integer_column=FALSE, leave_carrot=FALSE, reverse=FALSE)->Pseudo16S_table
	gsub("H133.Pseudo_16S", "H133" , Pseudo16S_table[,1])->Pseudo16S_table[,1]
	gsub("Pseudo", "p" , Pseudo16S_table[,1])->Pseudo16S_table[,1]
	gsub("_16S", "" , Pseudo16S_table[,1])->Pseudo16S_table[,1]
	
	
	#as.matrix(table(Pseudo16S_table[,2]))->Pseudo16S_table	
	

core$tip.label=gsub("_.*", "", core$tip.label)


#Turn tip labels into the actual V3V4 ASV sequences
Pseudo16S_table[match(core$tip.label, Pseudo16S_table[,1]),2]->associated_seqs_of_tiplabels

#Rename cultured seqs into ASVs
as.matrix(table(associated_seqs_of_tiplabels))->unique_seqs
unique_seqs[order(unique_seqs[,1], decreasing=TRUE),]->unique_seqs
cbind(names(unique_seqs), unique_seqs)->unique_seqs
rownames(unique_seqs)<-NULL
unique_seqs=cbind(unique_seqs, paste("ASV", 1:nrow(unique_seqs), sep=""))
library(RColorBrewer)
n <- nrow(unique_seqs)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
cbind(unique_seqs, col_vector[1:n])->unique_seqs


#Match the tip labels to an ASV number and rename the tips
#unique_seqs[match(core$tip.label, unique_seqs[,1]),3]->core$tip.label

#Match the tip labels to an ASV color and rename the tips
unique_seqs[match(associated_seqs_of_tiplabels, unique_seqs[,1]),4]->color


#classify seqs as "Pseudomonas_ASVs"




#color=rep("black", length(core$tip.label))
#"red"->color[grep("H133", core$tip.label)]
#"red"->color[grep("Zotu7_", core$tip.label)]
#"blue"->color[grep("Zotu51_", core$tip.label)]
#"green"->color[grep("Zotu299_", core$tip.label)]

plantmeta[match(core$tip.label, plantmeta[,1]),5]->matches
for(m in 1:length(matches)){
	if(is.na(matches[m]==TRUE)){
		core$tip.label[m]->matches[m]
	}
}		
matches->core$tip.label
core$tip.label=gsub("Athaliana.*", "************", core$tip.label)

#pdf("/Users/dlundberg/Documents/abt6/Pratchaya/PanX/TREE.pdf", width = 3, height = 10, useDingbats=FALSE)
	
	

	plot(core, lwd=1, tip.color=color,cex=.2)
	#nodelabels(cex=.2)
#dev.off()



for s in $(ls | grep "16S_of_")
do 
echo $s
cat $s >> all_16S_TK165.fa
done

