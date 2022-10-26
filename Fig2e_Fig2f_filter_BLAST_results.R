#Rscript print-args.R first "$filename"
#Rscript /ebio/abt6_projects7/bacterial_strain_analysis/code/filter_BLAST_results.R "Sph0349_S337H113_BLAST_nanopore.txt" "/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/AugustusAnnotation/Sphingomonas_BLAST_results"

args <- commandArgs(trailingOnly = TRUE)
  filename <- args[1]
  working_dir <- args[2]
  

matrix_format<-function(matrix, type=""){
	as.matrix(matrix)->matrix
	if(type=="" | type=="H"){
		if(ncol(matrix)==1){
			t(matrix)->matrix
		}
	}
	if(type=="V"){}
	return(matrix)
}	

allBLAST<-read.table(file=paste(working_dir, "/", filename, sep="", collapse=""), quote="", sep="\t")
names(allBLAST)=c("SAMPLE", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

allBLAST=allBLAST[which(allBLAST$pident>=40),]
allBLAST=allBLAST[which((allBLAST$qlen/allBLAST$slen)>=0.6),]
allBLAST=allBLAST[which((allBLAST$length/allBLAST$qlen)>=0.6),]

		s=unique(allBLAST$SAMPLE)
		#save sample as new table
		allBLAST[which(allBLAST$SAMPLE==s),]->extracted_sample
		#remove that sample temporarily from the full table
		allBLAST[which(allBLAST$SAMPLE!=s),]->allBLAST
		#take each gene in extracted_sample[,"sseqid"]
		do.call(rbind, strsplit(as.matrix(extracted_sample[,"sseqid"]), split="_____"))[,1]->genelist
		for(g in unique(genelist)){
			print(g)
			 extracted_sample[grep(paste(g, "_____", sep="", collapse=""), extracted_sample[,"sseqid"]),]->extracted_sample_gene
			 extracted_sample[grep(paste(g, "_____", sep="", collapse=""), extracted_sample[,"sseqid"], invert=TRUE),]->extracted_sample
			 matrix_format(extracted_sample_gene[which(extracted_sample_gene[,"pident"]==max(extracted_sample_gene[,"pident"]))[1],])->top_sample_gene
			 rbind(extracted_sample, top_sample_gene)->extracted_sample
		}


gsub("BLAST", "blast", filename)->filename
gsub(".txt", ".condensed.txt", filename)->filename
write.table(extracted_sample, paste(working_dir, "/", filename, sep="", collapse=""), quote = FALSE, row.names=FALSE, col.names=FALSE, sep="\t")


#
#/ebio/abt6_projects7/bacterial_strain_analysis/code/filter_BLAST_results.R.sh
###############################
#reserve running with 1 CPUs for this job
#$ -pe parallel 1
#request 4GB of RAM
#$ -l h_vmem=4G
#use /bin/bash to execute this script
#$ -S /bin/bash
#run this job from current working directory
#$ -cwd


filename=$1
working_dir=$2
Rscript /ebio/abt6_projects7/bacterial_strain_analysis/code/filter_BLAST_results.R $filename $working_dir
mkdir -p BestHits
mv *blast*condensed* BestHits
