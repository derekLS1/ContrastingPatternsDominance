

destination=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/PLASMIDS
plasmids=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/nanopore/Sph_nanopore_plasmids.fna

#cycle minimap through the nanopore genomes only
#cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/nanopore
#mkdir -p minimap_plasmids


#NANOPORE GENOMES
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/nanopore
for n in Sph*pilon.nanopore*.fasta
do 
name=$(sed "s/_pilon.nanopore.*.fasta//g" <(echo $n))
echo $name
nice -10 minimap2 $plasmids $n > "$destination"/"$name".plasmids.minimap2.paf
cat "$destination"/"$name".plasmids.minimap2.paf | cut -f1-12 > "$destination"/"$name".plasmids.minimap.paf
rm "$destination"/"$name".plasmids.minimap2.paf
sed -i "s/^/${name}_/g" "$destination"/"$name".plasmids.minimap.paf
done
cd $destination
#H133 GENOMES
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/genome_assemblies_HiSeq0133/Sphingomonas
for n in S*H*
do 
name=$n
nice -10 minimap2 $plasmids ./"$name"/"$name"_pilon.fasta  > "$destination"/"$name".plasmids.minimap2.paf
cat "$destination"/"$name".plasmids.minimap2.paf | cut -f1-12 > "$destination"/"$name".plasmids.minimap.paf
rm "$destination"/"$name".plasmids.minimap2.paf
sed -i "s/^/${name}_/g" "$destination"/"$name".plasmids.minimap.paf
done
cd $destination
#H113 GENOMES
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/genome_assemblies_HiSeq0113
for n in S*H*
do 
nice -10 minimap2 $plasmids ./"$name"/"$name"_pilon.fasta  > "$destination"/"$name".plasmids.minimap2.paf
cat "$destination"/"$name".plasmids.minimap2.paf | cut -f1-12 > "$destination"/"$name".plasmids.minimap.paf
rm "$destination"/"$name".plasmids.minimap2.paf
sed -i "s/^/${name}_/g" "$destination"/"$name".plasmids.minimap.paf
done
cd $destination
#REFSEQ
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/REFSEQ_assemblies
for n in GCF_*fna
do 
name=$(sed "s/GCF_/GCF/g" <(echo $n))
name=$(sed "s/_.*//g" <(echo $name))
name=$(sed "s/GCF/GCF_/g" <(echo $name))
nice -10 minimap2 $plasmids $n > "$destination"/"$name".plasmids.minimap2.paf
cat "$destination"/"$name".plasmids.minimap2.paf | cut -f1-12 > "$destination"/"$name".plasmids.minimap.paf
rm "$destination"/"$name".plasmids.minimap2.paf
sed -i "s/^/${name}_/g" "$destination"/"$name".plasmids.minimap.paf
done
cd $destination


#Remove files with 0 lines
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/PLASMIDS
for f in *paf
do
lines=$(wc -l $f | cut -d" " -f1) 
echo $lines
if $lines==0
then
   echo "0 lines"
else
   echo "more lines"
fi
done


#########
#open R

allPAF<-list.files(path = "/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/PLASMIDS", pattern = "paf")


condensed_mappings=matrix(ncol=12, nrow=0)
for(p in allPAF){
print(p)
paf=read.table(paste("/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/PLASMIDS/", p, sep="", collapse=""), sep="\t", quote="", comment.char = "")
rbind(condensed_mappings, paf)->condensed_mappings
}

condensed_mappings[, c(1, 6, 11, 12)]->condensed_mappings
names(condensed_mappings)=c("query", "target", "bases", "quality")

condensed_mappings=condensed_mappings[which(condensed_mappings$bases>30000),]


#rename seqs
gsub("GCF_", "GCF", condensed_mappings[,1])->condensed_mappings[,1]
gsub("_S", "-S", condensed_mappings[,1])->condensed_mappings[,1]
gsub("_.*", "", condensed_mappings[,1])->condensed_mappings[,1]
gsub("GCF", "GCF_", condensed_mappings[,1])->condensed_mappings[,1]

gsub("utg00000", "", condensed_mappings[,2])->condensed_mappings[,2]
gsub(":1.*pilon", "", condensed_mappings[,2])->condensed_mappings[,2]


