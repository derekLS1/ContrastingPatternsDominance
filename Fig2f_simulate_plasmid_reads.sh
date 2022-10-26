
destination=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/PLASMIDS/SIMULATE_READS

fsize=122
#NANOPORE GENOMES
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/nanopore
for n in Sph*pilon.nanopore*.fasta
do 
name=$(sed "s/_pilon.nanopore.*.fasta//g" <(echo $n))
echo $name
grep -v '^>' "$n" | tr -d '\n' | fold -w "$fsize" | nl -n rz -s "@" | sed "s/^/>${name}_/g" | sed -e 's/@/\n/g' > "$destination"/"$name"."$fsize".fasta
done
cd $destination

#H133 GENOMES
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/genome_assemblies_HiSeq0133/Sphingomonas
for n in S*H*
do 
name=$n
echo $name
grep -v '^>' ./"$name"/"$name"_pilon.fasta | tr -d '\n' | fold -w "$fsize" | nl -n rz -s "@" | sed "s/^/>${name}_/g"  | sed -e 's/@/\n/g' > "$destination"/"$name"."$fsize".fasta
done
cd $destination


#H113 GENOMES
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/genome_assemblies_HiSeq0113
for n in S*H*
do 
name=$n
echo $name
grep -v '^>' ./"$name"/"$name"_pilon.fasta | tr -d '\n' | fold -w "$fsize" | nl -n rz -s "@" | sed "s/^/>${name}_/g"  | sed -e 's/@/\n/g' > "$destination"/"$name"."$fsize".fasta
done
cd $destination

#REFSEQ
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/REFSEQ_assemblies
for n in GCF_*fna
do 
name=$(sed "s/GCF_/GCF/g" <(echo $n))
name=$(sed "s/_.*//g" <(echo $name))
name=$(sed "s/GCF/GCF_/g" <(echo $name))
echo $name
grep -v '^>' "$n" | tr -d '\n' | fold -w "$fsize" | nl -n rz -s "@" | sed "s/^/>${name}_/g"  | sed -e 's/@/\n/g' > "$destination"/"$name"."$fsize".fasta
done
cd $destination

#PLASMIDS
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/nanopore/EACH_PLASMID/
for n in Sph*pilon.fna
do
    name=$(sed "s/___.*//g" <(echo $n))
    echo $name
    grep -v '^>' "$n" | tr -d '\n' | fold -w "$fsize" | nl -n rz -s "@" | sed "s/^/>${name}_/g"  | sed -e 's/@/\n/g' > "$destination"/"$name"."$fsize".fasta
done
cd $destination


#################################################################################################### 
#################################################################################################### 
#################################################################################################### 
# MINIMAP
#################################################################################################### 
#################################################################################################### 
#################################################################################################### 

plasmids=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/nanopore/Sph_nanopore_plasmids.fna
destination=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/PLASMIDS/SIMULATE_READS
mkdir -p minimap_results
cd $destination
for n in *"$fsize".fasta
do 
echo $n
name=$(sed "s/.fasta//g" <(echo $n))
echo $name
nice -10 minimap2 $plasmids $n > "$destination"/minimap_results/"$name".tmp.paf
cat "$destination"/minimap_results/"$name".tmp.paf | cut -f1-12 > "$destination"/minimap_results/"$name".paf
rm "$destination"/minimap_results/"$name".tmp.paf
sed -i "s/^/${name}_/g" "$destination"/minimap_results/"$name".paf
done

#remove files with less than 1 line
cd "$destination"/minimap_results
find . -type f -exec awk -v x=1 'NR==x{exit 1}' {} \; -exec rm -f {} \;

###OPEN R


for i in *\.100\.paf
do
echo $i
cut -f1,6,11,12 $i >> all_100bp_condensed.paf
done


#########
#open R

#allPAF<-list.files(path = "/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/PLASMIDS/SIMULATE_READS/minimap_results", pattern = "\\.100.paf")


#condensed_mappings=matrix(ncol=12, nrow=0)
#for(p in allPAF){
#	print(p)
#	paf=read.table(paste("/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/PLASMIDS/SIMULATE_READS/minimap_results/", p, sep="", collapse=""), sep="\t", quote="", comment.char = "")
#	print(dim(paf))
#	rbind(condensed_mappings, paf)->condensed_mappings
#}
#condensed_mappings[, c(1, 6, 11, 12)]->condensed_mappings


condensed_mappings=read.table("/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/PLASMIDS/SIMULATE_READS/minimap_results/all_100bp_condensed.paf", sep="\t", quote="", comment.char = "")
names(condensed_mappings)=c("query", "target", "bases", "quality")

condensed_mappingsS=condensed_mappings

condensed_mappings=condensed_mappings[which(condensed_mappings$bases>=50),]
condensed_mappings=condensed_mappings[which(condensed_mappings$quality>=20),]




write.table(condensed_mappings, "/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/PLASMIDS/SIMULATE_READS/minimap_results/simulate_plasmid_minimap_100.txt", 
	quote =FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
            

#open R mac

scp dlundberg@burrito.eb.local:/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/PLASMIDS/SIMULATE_READS/minimap_results/simulate_plasmid_minimap_100.txt /Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/AUGUSTUS/


condensed_mappings=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/AUGUSTUS/simulate_plasmid_minimap_100.txt", sep="\t", quote="")


condensed_mappings=condensed_mappings[which(condensed_mappings$bases>=75),]
condensed_mappings=condensed_mappings[which(condensed_mappings$quality>=20),]


#rename seqs
gsub("GCF_", "GCF", condensed_mappings[,1])->condensed_mappings[,1]
gsub("_S", "-S", condensed_mappings[,1])->condensed_mappings[,1]
gsub("_.*", "", condensed_mappings[,1])->condensed_mappings[,1]
gsub("GCF", "GCF_", condensed_mappings[,1])->condensed_mappings[,1]
gsub("100.*", "100", condensed_mappings[,1])->condensed_mappings[,1]

gsub("utg00000", "", condensed_mappings[,2])->condensed_mappings[,2]
gsub(":1.*pilon", "", condensed_mappings[,2])->condensed_mappings[,2]


#make table
#NOW GET THE MINIMAP PLASMID RESULTS
condensed_mappings->plasmid_minimap
unique(plasmid_minimap[,2])->plasmids
unique(plasmid_minimap[,1])->genomes
matrix(data=0, nrow=length(genomes), ncol=length(plasmids))->plasmidPA

colnames(plasmidPA)=plasmids
rownames(plasmidPA)=genomes
for(r in 1:nrow(plasmid_minimap)){
	destROW=which(plasmid_minimap[r,1]==genomes)
	destCOL=which(plasmid_minimap[r,2]==plasmids)
	plasmidPA[destROW, destCOL] + plasmid_minimap[r,3] -> plasmidPA[destROW, destCOL]
}
		
plasmid_lengths=matrix(c(
 "Sph0137_S133H113_3c", "263190",
 "Sph0018_S18H113_2c","558423", 
 "Sph0140_S136H113_2c","58686", 
 "Sph0137_S133H113_2c","44079", 
 "Sph0240_S230H113_3c","89037", 
 "Sph0131_S127H113_3c","48992", 
 "Sph0018_S18H113_4c","39688",  
 "Sph0349_S337H113_2c","142166",
 "Sph0219_S190H113_2c","53300", 
 "Sph0240_S230H113_2c","61217", 
 "Sph0018_S18H113_3c","62857",  
 "Sph0226_S216H113_2c","57877", 
 "Sph0223_S213H113_2c","62051", 
 "Sph0137_S133H113_4c","60801", 
 "Sph0131_S127H113_2c","237536",
 "Sph0136_S132H113_2c","66142"), ncol=2, byrow=TRUE)
	
as.numeric(plasmid_lengths[match(colnames(plasmidPA), plasmid_lengths[,1]),2])->matched_plasmid_length
for(c in 1:ncol(plasmidPA)){
	100*plasmidPA[,c] / matched_plasmid_length[c] -> plasmidPA[,c] 
}
	
length(which(plasmidPA>0))
length(which(plasmidPA>80))
	
	

	
