#copy all the individual FASTA collections to burrito
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/GENES

rm  matches.txt
#rename all the headers so that the file name (gene name) precedes everything else 
for f in ./*/*fasta
do 
    echo $f
    gene=$(sed "s/.fasta//g" <(echo $f))
    gene=$(sed "s/.*\///g" <(echo $gene))
    echo $gene
    cat $f > /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/GENES/REFORMATTED/"$gene"_reformatted.fasta    
    cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/GENES/REFORMATTED
    
    awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" }END { printf "%s", n }' "$gene"_reformatted.fasta  > "$gene"_reformatted.fna
	rm "$gene"_reformatted.fasta 
	
    grep "|ref|" "$gene"_reformatted.fna    -A 1 >> "$gene"_reformatted.fasta 
    grep "|gbDSL|" "$gene"_reformatted.fna -A 1 >> "$gene"_reformatted.fasta 
    grep "|embDSL|" "$gene"_reformatted.fna -A 1 >> "$gene"_reformatted.fasta 
    grep "partial" -v "$gene"_reformatted.fna > "$gene"_reformatted.fasta 
    rm "$gene"_reformatted.fna

    usearch11 -fastx_uniques "$gene"_reformatted.fasta -fastaout "$gene"_uniques.fa -sizeout -relabel Uniq
    usearch11 -sortbylength "$gene"_uniques.fa -fastaout "$gene"_sorted.fasta -minseqlength 20
	usearch11 -cluster_smallmem "$gene"_sorted.fasta  -minsize 1 -id 0.95 -centroids "$gene"_otus95.fa
	usearch11  -sortbysize "$gene"_otus95.fa -fastaout "$gene"_otus95_sorted.fa -minsize 1
    rm "$gene"_reformatted.fasta
    rm "$gene"_uniques.fa
    rm "$gene"_sorted.fasta
    rm "$gene"_otus95.fa
    
    
    sed -i "s/^>/>"$gene"_____/g" "$gene"_otus95_sorted.fa   
    awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" }END { printf "%s", n }' "$gene"_otus95_sorted.fa   > "$gene"_otus95_sorted.fna
	head "$gene"_otus95_sorted.fna -n 40 > "$gene"_otus95_sorted.fa 
	rm "$gene"_otus95_sorted.fna
	echo "$gene" $(wc -l "$gene"_otus95_sorted.fa ) >> matches.txt  
	cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/GENES
done


#concatenate all the files into one fasta file
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/GENES/REFORMATTED
rm sphingo_genes_refseq.fna
for f in *otus95_sorted.fa 
do 
echo $f
cat $f >> sphingo_genes_refseq.fna
rm $f
done



#make a BLAST database
makeblastdb=/ebio/abt6_projects7/bacterial_strain_analysis/code/ncbi-blast-2.9.0+/bin/makeblastdb
$makeblastdb -in sphingo_genes_refseq.fna -dbtype prot
$makeblastdb -in gyrB.fna -dbtype nucl
cp sphingo_genes_refseq* /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/BLAST_db
cp gyrB* /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/BLAST_db


#Now can run the following, which calls a shell script below to run AUGUSTUS and BLAST 
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################

#be sure to comment out augustus, if not running it again

#run on all the H133 genomes:
#1=genome basename
#2=directory where the fasta file can be found
assembly_dir=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/genome_assemblies_HiSeq0133/Sphingomonas
blast_db=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/BLAST_db/sphingo_genes_refseq.fna
#blast_db=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/BLAST_db/gyrB.fna
cd $assembly_dir
for g in *133
do
echo $g
#qsub -N $g -o $(pwd)/out$g.txt -e $(pwd)/err$g.txt -cwd /ebio/abt6_projects7/bacterial_strain_analysis/code/augustus_and_blast_effectors.sh "$g" "$assembly_dir"/"$g" "$blast_db" "$g"_pilon.fasta
qsub -N $g -o $(pwd)/out$g.txt -e $(pwd)/err$g.txt -cwd /ebio/abt6_projects7/bacterial_strain_analysis/code/augustus_and_blast_effectors_fabrice.sh "$g" "$assembly_dir"/"$g" "$blast_db" "$g"_pilon.fasta
done
#run on all the H113 genomes:
#1=genome basename
#2=directory where the fasta file can be found
assembly_dir=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/genome_assemblies_HiSeq0113
blast_db=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/BLAST_db/sphingo_genes_refseq.fna
#blast_db=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/BLAST_db/gyrB.fna
cd $assembly_dir
for g in *113
do
echo $g
#qsub -N $g -o $(pwd)/out$g.txt -e $(pwd)/err$g.txt -cwd /ebio/abt6_projects7/bacterial_strain_analysis/code/augustus_and_blast_effectors.sh "$g" "$assembly_dir"/"$g" "$blast_db" "$g"_pilon.fasta
qsub -N $g -o $(pwd)/out$g.txt -e $(pwd)/err$g.txt -cwd /ebio/abt6_projects7/bacterial_strain_analysis/code/augustus_and_blast_effectors_fabrice.sh "$g" "$assembly_dir"/"$g" "$blast_db" "$g"_pilon.fasta
done
#run on all the nanopore genomes:
#1=genome basename
#2=directory where the fasta file can be found
assembly_dir=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/nanopore
blast_db=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/BLAST_db/sphingo_genes_refseq.fna
cd $assembly_dir
for g in Sph*pilon.nanopore*.fasta
do
echo $g
basename=$(sed "s/_pilon.nanopore.*fasta//g" <(echo $g))
echo $basename
#basename=$1 #assembly_dir=$2 blast_database=$3 #genome_fasta=$4
qsub -N $g -o $(pwd)/out$g.txt -e $(pwd)/err$g.txt -cwd /ebio/abt6_projects7/bacterial_strain_analysis/code/augustus_and_blast_effectors.sh "$basename" "$assembly_dir" "$blast_db" "$g"
done
#run on all the REFSEQ genomes:
#1=genome basename
#2=directory where the fasta file can be found
assembly_dir=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/REFSEQ_assemblies
blast_db=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/BLAST_db/sphingo_genes_refseq.fna
cd $assembly_dir
for genome in GCF_*fna
do
echo $genome
basename=$(sed "s/GCF_/GCF/g" <(echo $genome))
basename=$(sed "s/_.*//g" <(echo $basename))
basename=$(sed "s/GCF/GCF_/g" <(echo $basename))
echo $basename
#basename=$1 #assembly_dir=$2 blast_database=$3 #genome_fasta=$4
qsub -N $genome -o $(pwd)/out$genome.txt -e $(pwd)/err$g.txt -cwd /ebio/abt6_projects7/bacterial_strain_analysis/code/augustus_and_blast_effectors.sh "$basename" "$assembly_dir" "$blast_db" "$genome"
done
#run on all the nanopore PLASMIDS separately
assembly_dir=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/nanopore/EACH_PLASMID
blast_db=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/BLAST_db/sphingo_genes_refseq.fna
cd $assembly_dir
for g in Sph*pilon.fna
do
basename=$(sed "s/___.*fna//g" <(echo $g))
echo $basename
#basename=$1 #assembly_dir=$2 blast_database=$3 #genome_fasta=$4
qsub -N "$g" -o $(pwd)/out$g.txt -e $(pwd)/err$g.txt -cwd /ebio/abt6_projects7/bacterial_strain_analysis/code/augustus_and_blast_effectors.sh "$basename" "$assembly_dir" "$blast_db" "$g"
done
#DO SOME H133 GENOMES ON PSEUDOMONAS EFFECTORS
#assembly_dir=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/genome_assemblies_HiSeq0133/Sphingomonas
#blast_db=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/BLAST_db/Hop-effector-db_reference.fasta
#cd $assembly_dir
#for g in S139H133 S217H133 S114H133 S134H133 S243H133 S241H133 
#do
#echo $g
#qsub -N $g -o $(pwd)/out$g.txt -e $(pwd)/err$g.txt -cwd /ebio/abt6_projects7/bacterial_strain_analysis/code/augustus_and_blast_effectors.sh "$g" "$assembly_dir"/"$g" "$blast_db" "$g"_pilon.fasta
#done
#run on all the H113 genomes:
#
#makeblastdb=/ebio/abt6_projects7/bacterial_strain_analysis/code/ncbi-blast-2.9.0+/bin/makeblastdb
#$makeblastdb -in /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/BLAST_db/Hop-effector-db_reference.fasta -dbtype prot



#after complete, add filename to first column for all
#put master BLAST files in the AUGUSTUS directory
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
destination_dir=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/Sphingomonas_BLAST_results


assembly_dir=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/genome_assemblies_HiSeq0133/Sphingomonas
cd $assembly_dir
for g in *133
do
echo $g
cd "$assembly_dir"/"$g"
#nawk '{print FILENAME"\t"$0}' "$g"_BLAST_results.txt > "$g"_BLASTtemp.txt
nawk '{print FILENAME"\t"$0}' "$g"_BLAST_results_fabrice.txt > "$g"_BLASTtemp.txt
#cat "$g"_BLASTtemp.txt >> "$assembly_dir"/all_BLAST_H133.txt
#rm "$g"_BLASTtemp.txt
#mv "$g"_BLASTtemp.txt "$destination_dir"/"$g"_BLAST_H133.txt
mv "$g"_BLASTtemp.txt "$destination_dir"/"$g"_BLAST_H133_fabrice.txt
done
cd "$assembly_dir"
rm err*txt
rm out*txt
#cp all_BLAST_H133.txt /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/Sphingomonas_BLAST_results
#
assembly_dir=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/genome_assemblies_HiSeq0113
cd $assembly_dir
for g in *113
do
echo $g
cd "$assembly_dir"/"$g"
nawk '{print FILENAME"\t"$0}' "$g"_BLAST_results.txt > "$g"_BLASTtemp.txt
#cat "$g"_BLASTtemp.txt >> $assembly_dir/all_BLAST_H113.txt
#rm "$g"_BLASTtemp.txt
#mv "$g"_BLASTtemp.txt "$destination_dir"/"$g"_BLAST_H113.txt
mv *fabrice.txt "$destination_dir"/"$g"_BLAST_H113_fabrice.txt
done
cd "$assembly_dir"
rm err*txt
rm out*txt
#cp all_BLAST_H113.txt /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/Sphingomonas_BLAST_results
#
assembly_dir=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/nanopore
cd $assembly_dir
for g in Sph*pilon.nanopore*.fasta
do
basename=$(sed "s/_pilon.nanopore.*fasta//g" <(echo $g))
echo $basename
nawk '{print FILENAME"\t"$0}' "$basename"_BLAST_results.txt > "$basename"_BLASTtemp.txt
mv "$basename"_BLASTtemp.txt "$destination_dir"/"$basename"_BLAST_nanopore.txt
#cat "$basename"_BLASTtemp.txt >> $assembly_dir/all_BLAST_nanopore.txt
#rm "$basename"_BLASTtemp.txt
#cp all_BLAST_nanopore.txt /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/Sphingomonas_BLAST_results
done
#
assembly_dir=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/REFSEQ_assemblies
cd $assembly_dir
for genome in GCF_*fna
do
basename=$(sed "s/GCF_/GCF/g" <(echo $genome))
basename=$(sed "s/_.*//g" <(echo $basename))
basename=$(sed "s/GCF/GCF_/g" <(echo $basename))
echo $basename
nawk '{print FILENAME"\t"$0}' "$basename"_BLAST_results.txt > "$basename"_BLASTtemp.txt
mv "$basename"_BLASTtemp.txt "$destination_dir"/"$basename"_BLAST_REFSEQ.txt
#cat "$basename"_BLASTtemp.txt >> $assembly_dir/all_BLAST_REFSEQ.txt
#rm "$basename"_BLASTtemp.txt
#cp all_BLAST_REFSEQ.txt /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/Sphingomonas_BLAST_results
done
#PLASMIDS
assembly_dir=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/nanopore/EACH_PLASMID
cd $assembly_dir
for g in Sph*pilon.fna
do
basename=$(sed "s/___.*fna//g" <(echo $g))
echo $basename
nawk '{print FILENAME"\t"$0}' "$basename"_BLAST_results.txt > "$basename"_BLASTtemp.txt
mv "$basename"_BLASTtemp.txt "$destination_dir"/"$basename"_BLAST_PLASMIDS.txt
#cat "$basename"_BLASTtemp.txt >> $assembly_dir/all_BLAST_nanopore.txt
#rm "$basename"_BLASTtemp.txt
#cp all_BLAST_nanopore.txt /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/Sphingomonas_BLAST_results
done
#BLAST S380H113 genomes, not the aa file 
assembly_dir=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/bacterial_genome_assemblies/nanopore/EACH_PLASMID
cd $assembly_dir
for g in Sph*pilon.fna
do
basename=$(sed "s/___.*fna//g" <(echo $g))
echo $basename
nawk '{print FILENAME"\t"$0}' "$basename"_BLAST_results.txt > "$basename"_BLASTtemp.txt
mv "$basename"_BLASTtemp.txt "$destination_dir"/"$basename"_BLAST_PLASMIDS.txt
#cat "$basename"_BLASTtemp.txt >> $assembly_dir/all_BLAST_nanopore.txt
#rm "$basename"_BLASTtemp.txt
#cp all_BLAST_nanopore.txt /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/Sphingomonas_BLAST_results
done



#filter and pick best hits for the blast results
working_dir=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/Sphingomonas_BLAST_results
cd $working_dir
#now run Rscript to pick the best hit in each
for filename in *_BLAST*txt
do 
echo $filename
qsub -N "$filename" -o $(pwd)/out.txt -e $(pwd)/err.txt -V -cwd /ebio/abt6_projects7/bacterial_strain_analysis/code/filter_BLAST_results.R.sh "$filename" "$working_dir"
done


#combine all the blast results into one file
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/Sphingomonas_BLAST_results/BestHits
rm allBLASTcondensed.txt
for f in *condensed*
do
cat $f >> allBLASTcondensed.txt
done



#send to MAC
scp dlundberg@burrito.eb.local:/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/Sphingomonas_BLAST_results/BestHits/allBLASTcondensed.txt /Users/dlundberg/Documents/abt6/Pratchaya/AUGUSTUS


for f in *113_fabrice*
do
cat $f >> fabricecombined2.txt
done


#open R on mac
source("/Users/dlundberg/Documents/abt6/scripts/R/functions/microbiome_custom_functions.R")
allBLAST<-read.table(file="/Users/dlundberg/Documents/abt6/Pratchaya/AUGUSTUS/allBLASTcondensed.txt", quote="", sep="\t")
names(allBLAST)=c("SAMPLE", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")


#make PA table
allBLAST[,c(1, 3, 4)]->PAtable
#rep(1, nrow(PAtable))->PAtable[,3]
as.numeric(as.matrix(PAtable[,3]))->PAtable[,3]
unique(PAtable[,1])->unique_samples
as.matrix(unique(PAtable[,2]))->unique_blast_hits
matrix(nrow=nrow(unique_blast_hits), ncol=length(unique_samples))->PAtable2
for(s in unique_samples){
	which(s==unique_samples)->current_column
	which(s==PAtable[,1])->current_rows_in_PAtable
	#subset just the sample S as a new table to avoid confusion
	as.matrix(PAtable[current_rows_in_PAtable,])->PAtable_current_sample
	#PAtable_current_sample[match(unique_blast_hits, PAtable_current_sample[,2]),3]->PAtable_current_sample
	as.numeric(PAtable_current_sample[,3])->PAtable2[match(PAtable_current_sample[,2], unique_blast_hits),current_column]
}
PAtable2[is.na(PAtable2)]<-0
rownames(PAtable2)=unique_blast_hits
colnames(PAtable2)=unique_samples

#condense the PA table such that each gene, regardless of which strain it came from, is a row.
do.call(rbind, strsplit(as.matrix(unique_blast_hits), split="_____"))[,1]->genelist2
matrix(ncol=ncol(PAtable2), nrow=length(unique(genelist2)))->PAtable2condensed
for(g in unique(genelist2)){
	current_row=which(g==unique(genelist2))
	apply(matrix_format(PAtable2[which(genelist2==g),]), 2, max)->PAtable2condensed[current_row,]
}
unique(genelist2)->rownames(PAtable2condensed)
colnames(PAtable2)->colnames(PAtable2condensed)
	
	
	

#put the gene lists in order
library("ape")
core=read.tree(file = "/Users/dlundberg/Documents/abt6/Pratchaya/PanX/Sphingo_runs/Sph_buscoGT85/strain_tree.nwk") 
basefolder="/Users/dlundberg/Documents/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS/"
	core=read.tree(file = paste(basefolder, "strain_tree.nwk", sep="", collapse=""))
		#drop those also duplicated by nanopore sequencing
	
	compare_ILMN_to_nanopore=FALSE #TRUE #TRUE# FALSE #TRUE #FALSE  #TRUE #TRUE #TRUE #FALSE
	if(compare_ILMN_to_nanopore==TRUE){
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
	if(compare_ILMN_to_nanopore!=TRUE){
		Drop=which(match(core$tip.label,c("S337H113","S132H113","S136H113",
										"S380H113","S237H113","S230H113",
										"S190H113","S213H113","S18H113",
										"S127H113","S133H113","S216H113"), nomatch=0)>0)
		drop.tip(core, as.numeric(Drop))->core
	}
	
gsub("_BLAST_results.txt" , "", colnames(PAtable2condensed))->colnames(PAtable2condensed)

t(PAtable2condensed)->PAtable2condensed

rownames(PAtable2condensed)->samplenames
#gsub("GCF", "GCF_", samplenames)->samplenames

if(compare_ILMN_to_nanopore!=TRUE){
	gsub("_S.*", "", samplenames)->samplenames
}

#"Sph0018_S18H113","Sph0131_S127H113","Sph0136_S132H113","Sph0137_S133H113","Sph0140_S136H113","Sph0184_S380H113","Sph0219_S190H113","Sph0223_S213H113","Sph0226_S216H113","Sph0240_S230H113","Sph0247_S237H113","Sph0349_S337H113")



photosynthesis=c( "bchE", "acsF", "bchC", "bchL", "bchX", "bchY", "bchZ")
toxinantitoxin=c("brnT","hicA","pemK","vapC","yoeB","ccdb","hipA","ratA","yafQ","fitB","parD","relE","yhaV")
#toxinantitoxin=c("brnT", "ccdb", "fitB", "hicA", "hipA" ,"parD" ,"pemK" ,"ratA", "relE" ,"vapC" ,"yafQ" ,"yhaV", "yoeB")
T4SS=c("virD4", "virB1", "virB2",  "virB3", "virB4", "virB5", "virB6","virB7","virB8","virB9","virB10", "virB11")
T3SS=c("hrcC", "hrcJ", "hrcN", "hrcR", "hrcS", "hrcT","hrcU", "hrcV")#, "escR", "hrcC_Ps", "hrcT_Ps", "hrcV_Ps", "hrpH_Ps", "hrpJ_Ps", "hrpL_Ps", "hrpZ_Ps")
#flagella=c(		fliF  	fliL    "fliP"	fliQ   fliR   "flhB"   "flhA"
flagella=c("flgE", "flgF", "flgG", "flhA", "flhB", "fliC", "fliD","fliE", "fliP")
T1SS=c("atrB", "atrD")
T2SS=c("epsF", "epsM", "gspC","gspF","gspI","gspL","gspD","gspG","gspJ","gspM","gspE","gspH","gspK")
T6SS=c("hcp","vgrG","tssA","tssB","tssC","tssE","tssF","tssG","tssH", "tssK","tssL","tssM")
Gellan_biosynthesis=c("gelG","gelS","gelR","gelQ","gelI","gelK","gelL","gelJ","gelF","gelD","gelC","gelE","gelM","gelN","gelB","gelA")
Duitan_biosynthesis=c("dpsS","dpsG","dpsR","dpsQ","dpsI","dpsK","dpsL","dpsJ","dpsF","dpsD","dpsC","dpsE","dpsM","dpsN","dpsB")
Sphingan_biosynthesis=c("spsG","spsS","spsR","spsQ","spsI","spsK","spsL","spsJ","spsF","spsD","spsC","spsE","spsM","spsB")
Carotenoid_biosynthesis=c("LOG","crtB","crtG","crtI","crtW","crtX","crtY","crtZ")
Phosphate_metabolism=c("oprO","phoQ","phoU","pstA","pstB","pstC","pstS")
H2S_production=c("cysA","cysP","cysT","cysW")
Phosphate_metabolism=c("cysA","cysP","cysQ","cysT","cysW")
IAA_production=c("trpA","trpB","trpD","trpF", "iaaM_PS", "iaaL_Ps", "iaah_Ps")
Nitrogen_fixation=c("nifH", "nifD", "nifK", "nifE", "nifN", "nifB", "nifA","nifJ","nifQ","nifS","nifT","nifU","nifV","nifW","nifX","nifZ")
	#https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-162
	
#virB1-virB11

library(gplots)

pal1 <- c(colorRampPalette(c("black","#27004E", "#340069","#6A1B9A", "#CB018E", "#EB6E99", "#FFBC7A", "#FFD6A9"))(n=8))
#pal1 <- c(colorRampPalette(c("black","#27004E", "#340069","#6A1B9A", "#F28149", "#F49C60", "#FFBC7A", "#FFD6A9"))(n=8))

pal1[length(pal1):1]->pal1
pdf(file=paste("/Users/dlundberg/Documents/abt6/Pratchaya/AUGUSTUS/BLAST_pdfs/", "ColorKey", "_noIBVSS.pdf", sep="", collapse=""), width = 3, height = 3, useDingbats=FALSE)
image(1:length(pal1), 1, as.matrix(1:length(pal1)), 
      col=pal1, xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")  
dev.off()		

#FFFEF4,  "#FFECB3"
	

focus=c("flagella")
genelist=vector(length=0)
for(f in focus){
	c(genelist, get(f))->genelist
}



if(compare_ILMN_to_nanopore==TRUE){
    PAtable2condensed[match(core$tip.label, samplenames),]->PA
}
if(compare_ILMN_to_nanopore!=TRUE){
    PAtable2condensed[match(gsub("_S.*", "", core$tip.label), samplenames),]->PA
}

	PA[,match(genelist, colnames(PA), nomatch=NA)]->PA
	0->PA[is.na(PA)]
	
	#PA[,order(colSums(PA))]->PA
	
	#find ones that are full in most columns
	tmp=PA
	tmp[PA>0]<-1
	PA[which(rowSums(tmp)>=7),]->tmp2
	tmp2[grep("Sph", rownames(tmp2)),]
	
	
focus2="plasmid"
if(focus2=="plasmid"){
	condensed_mappings=read.table("/Users/dlundberg/Documents/abt6/Pratchaya/AUGUSTUS/simulate_plasmid_minimap_2000.txt", sep="\t", quote="")	
	condensed_mappings=condensed_mappings[which(as.numeric(condensed_mappings$bases)>=200),]
	condensed_mappings=condensed_mappings[which(as.numeric(condensed_mappings$quality)>=20),]
	names(condensed_mappings)->futurecolnames
	as.matrix(condensed_mappings)->condensed_mappings
	
	#rename seqs
	gsub("\\.100.*", "", condensed_mappings[,1])->condensed_mappings[,1]
	gsub("_GCF_.*", "", condensed_mappings[,1])->condensed_mappings[,1]
	#gsub("_S", "-S", condensed_mappings[,1])->condensed_mappings[,1]
	#gsub("_.*", "", condensed_mappings[,1])->condensed_mappings[,1]
	gsub("\\.2000", "", condensed_mappings[,1])->condensed_mappings[,1]
	gsub("H113.*", "H113", condensed_mappings[,1])->condensed_mappings[,1]
	gsub("H133.*", "H133", condensed_mappings[,1])->condensed_mappings[,1]
	gsub("utg0*", "", condensed_mappings[,1])->condensed_mappings[,1]
	gsub("utg0*", "", condensed_mappings[,2])->condensed_mappings[,2]
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
		as.numeric(plasmidPA[destROW, destCOL]) + as.numeric(plasmid_minimap[r,3]) -> plasmidPA[destROW, destCOL]
	}
			
		
			
	plasmid_lengths=matrix(c(
	"Sph0226_S216H113_2c","57877",  
	  	"Sph0137_S133H113_2c","44079",
	  	"Sph0137_S133H113_4c","60801",
	 	"Sph0137_S133H113_3c", "263190",
	"Sph0131_S127H113_2c","237536",
		"Sph0131_S127H113_3c","48992", 
	"Sph0018_S18H113_4c","39688", 
	"Sph0018_S18H113_3c","62857",  
	"Sph0018_S18H113_2c","558423",  
	 	"Sph0223_S213H113_2c","62051",
	"Sph0219_S190H113_2c","53300", 
		"Sph0240_S230H113_2c","61217", 
		"Sph0240_S230H113_3c","89037", 
	"Sph0140_S136H113_2c","58686", 
		"Sph0136_S132H113_2c","66142",
	"Sph0349_S337H113_2c","142166"
	 ), ncol=2, byrow=TRUE)
	

	
	#put plasmids in above order
	plasmidPA[,match(plasmid_lengths[,1], colnames(plasmidPA))]->plasmidPA
	
	for(c in 1:ncol(plasmidPA)){
		100*plasmidPA[,c] / as.numeric(plasmid_lengths[c,2]) -> plasmidPA[,c] 
	}
		
	#plasmidPA[match(core$tip.label, rownames(plasmidPA)),]->PA	
	plasmidPA[match(gsub("_nanopore.*", "", core$tip.label), rownames(plasmidPA)),]->PA	
	
	PA[is.na(PA)]<-0
	
	
	#plasmidPA[grep("c", rownames(plasmidPA)),]
	
}

	
date=format(Sys.Date(), format="%Y%m%d")

if(compare_ILMN_to_nanopore==FALSE){
	sepcolor=NA
	sepwidth=c(0,0)
	colsep=NA
	rowsep=NA
	version="FULL"
	breaks = c(0, 30, 40, 50, 60, 70, 80, 90, 100)
}
if(compare_ILMN_to_nanopore==TRUE){
	sepcolor="#FFFFFF"
	sepwidth=c(.01,.01)
	colsep=1:ncol(PA)
	rowsep=1:nrow(PA)
	version="nano"
	breaks = c(0, 30, 40, 50, 60, 70, 80, 90, 100)
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
}

if(focus2=="plasmid"){
	#PA[PA>40]=100
	focus="plasmid"
	breaks=c(0, 10, 15, 20, 25, 30, 35, 40, 100)
	breaks = c(0, 30, 40, 50, 60, 70, 80, 90, 100)
}

pdf(file=paste("/Users/dlundberg/Documents/abt6/Pratchaya/AUGUSTUS/BLAST_pdfs/", focus, "_", date, "_", version, "_noIBVSS.pdf", sep="", collapse=""), width = ncol(PA)*0.3, height = 5, useDingbats=FALSE)
    tmp.pdf<- heatmap.2(PA[nrow(PA):1,], density.info = "none",
                    trace = "none", dendrogram="none", main=paste(focus, "_", length(genelist),sep="", collapse=""),
                    col = pal1, Colv=FALSE, Rowv=FALSE, key=FALSE, 
                    key.title = NA, key.xlab = NA, sepcolor=sepcolor,
                    breaks = breaks,
                    labCol = F, labRow=F, sepwidth=sepwidth, 
                    colsep=colsep, rowsep=rowsep)      
dev.off()		
   tmp.pdf<- heatmap.2(PA[nrow(PA):1,], density.info = "none",
                    trace = "none", dendrogram="none", main=paste(focus, "_", length(genelist),sep="", collapse=""),
                    col = pal1, Colv=FALSE, Rowv=FALSE, key=FALSE, 
                    key.title = NA, key.xlab = NA, sepcolor=NA,
                    breaks = c(0, 30, 40, 50, 60, 70, 80, 90, 100),
                    labCol = F, labRow=F, sepwidth=c(0,0))
  

#20210113 calculate how many of the 410 sphingo genomes have plasmid in them
Sph410=c("S100H113","S101H113","S102H113","S103H113","S105H113","S108H113","S109H113","S110H113","S111H113","S116H113","S117H113","S118H113","S119H113","S11H113","S120H113","S121H113","S122H113","S123H113","S124H113","S125H113","S126H113","S127H113","S128H113","S129H113","S130H113","S131H113","S132H113","S133H113","S134H113","S135H113","S136H113","S137H113","S138H113","S139H113","S140H113","S141H113","S142H113","S143H113","S144H113","S145H113","S145H133","S146H113","S146H133","S147H113","S148H113","S148H133","S149H133","S150H113","S150H133","S151H113","S151H133","S152H113","S152H133","S153H113","S153H133","S154H113","S154H133","S155H113","S155H133","S156H133","S157H113","S157H133","S158H113","S158H133","S159H133","S160H113","S160H133","S161H113","S161H133","S162H133","S163H113","S163H133","S164H113","S164H133","S165H133","S166H113","S166H133","S174H113","S175H113","S176H113","S177H113","S178H113","S179H113","S180H113","S181H113","S182H113","S183H113","S184H113","S186H113","S187H113","S188H113","S189H113","S18H113","S190H113","S193H113","S193H133","S194H113","S195H113","S196H113","S197H113","S197H133","S198H113","S199H113","S199H133","S200H113","S201H113","S202H113","S202H133","S203H113","S204H113","S204H133","S205H113","S206H113","S206H133","S207H113","S207H133","S208H113","S209H113","S210H113","S212H113","S212H133","S213H113","S214H133","S215H113","S216H113","S217H113","S218H113","S218H133","S219H113","S219H133","S21H113","S220H113","S221H113","S222H113","S222H133","S223H113","S224H113","S225H113","S225H133","S226H113","S226H133","S227H113","S228H113","S228H133","S229H113","S230H113","S231H113","S231H133","S232H113","S232H133","S233H113","S234H113","S235H113","S235H133","S236H113","S236H133","S237H113","S238H113","S239H113","S239H133","S240H113","S241H113","S241H133","S242H113","S243H113","S243H133","S244H113","S245H113","S245H133","S246H113","S246H133","S247H113","S247H133","S248H113","S249H113","S249H133","S250H113","S250H133","S251H113","S251H133","S252H113","S252H133","S253H113","S253H133","S254H113","S254H133","S255H113","S256H113","S257H113","S257H133","S258H113","S258H133","S259H113","S259H133","S25H113","S260H133","S261H113","S262H113","S262H133","S263H113","S263H133","S264H113","S264H133","S265H113","S265H133","S266H113","S266H133","S267H113","S267H133","S268H113","S269H113","S270H113","S270H133","S271H113","S271H133","S272H133","S273H113","S273H133","S274H113","S275H113","S276H113","S276H133","S277H113","S277H133","S278H113","S278H133","S279H113","S27H113","S280H113","S280H133","S281H113","S282H113","S282H133","S283H113","S284H113","S284H133","S285H113","S286H113","S286H133","S287H113","S287H133","S288H133","S289H113","S291H113","S292H113","S293H113","S294H113","S295H113","S296H113","S297H113","S298H113","S299H113","S300H113","S301H113","S302H113","S304H113","S306H113","S307H113","S308H113","S309H113","S310H113","S311H113","S312H113","S314H113","S315H113","S316H113","S317H113","S318H113","S319H113","S320H113","S321H113","S322H113","S323H113","S324H113","S325H113","S326H113","S327H113","S328H113","S329H113","S330H113","S331H113","S332H113","S333H113","S334H113","S335H113","S336H113","S337H113","S338H113","S339H113","S33H113","S340H113","S341H113","S342H113","S343H113","S344H113","S345H113","S346H113","S347H113","S349H113","S350H113","S351H113","S352H113","S353H113","S354H113","S355H113","S356H113","S357H113","S358H113","S359H113","S35H113","S360H113","S361H113","S362H113","S363H113","S364H113","S365H113","S366H113","S367H113","S368H113","S369H113","S370H113","S371H113","S373H113","S374H113","S377H113","S378H113","S380H113","S381H113","S44H113","S49H113","S51H113","S57H113","S59H113","S61H113","S65H113","S6H113","S76H113","S97H113","S98H113","S99H113","S102H133","S103H133","S106H133","S108H133","S200H133","S203H133","S209H133","S215H133","S227H133","S229H133","S230H133","S233H133","S234H133","S237H133","S244H133","S194H133","S195H133","S198H133","S210H133","S220H133","S221H133","S238H133","S248H133","S281H133","S109H133","S110H133","S111H133","S113H133","S122H133","S128H133","S130H133","S131H133","S132H133","S196H133","S205H133","S211H133","S213H133","S217H133","S223H133","S224H133","S240H133","S242H133","S255H133","S261H133","S268H133","S269H133","S274H133","S283H133","S133H133","S134H133","S135H133","S137H133","S138H133","S139H133","S140H133","S141H133","S142H133","S143H133","S114H133","S115H133","S116H133","S117H133","S118H133","S119H133","S120H133","S121H133","S123H133","S124H133","S125H133","S126H133")
PA=PA[match(Sph410, gsub(".*_", "", rownames(PA))),]

PA[which(PA<50)]<-0





