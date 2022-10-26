

## REDUCE MAPPING GENOMES TO REMOVE REDUNDANCY. GO BY MASH DISTANCE?  CHOOSE GENOMES THAT DIFFER
## BY AT LEAST 0.5% MASH DISTANCE














#######RENAME HEADERS: Before starting with the reference genome using BWA programme, we need to prepare the different genomes#######
#######Since we want to call the number of reads that map to each of the different genomes (decoy and others)#######
#######We will rename the headings of the contigs. Therefore, they contain all the information about the genome references they come from and thos databases (Refseq and decoy)#######
#######Talia's genomes as local reference genomes#######
cd /ebio/abt6_projects9/Pseudomonas_diversity/data/submitted_genomes
#######UNZIP, first of all, otherwise it does not work#######
gzip -d #file
#######Extract specific genomes from the folder (e.g. Sphingomonas from the At-LSPHERE):#######
for i in $(cat Pseud.txt)
do
echo "$i"_*
mv "$i"_* ./Psy_At-ROOT/
done
#######And sorry for the spam, here it comes the actual change (the headers)#######
#######If we don't use double quotes for the sed command reads literally, we have to use " instead of ' to read varaibles#######
for o in *.fasta
do
o=$(sed 's/.pilon.contigs_renamed.fasta//' <(echo $o))
echo $o
sed -i "s/>LOCAL165_/>nonOTU5_LOCAL165_/g" "$o".pilon.contigs_renamed.fasta
done

for o in GCF_*_genomic.fna
do
genome=$(sed 's/_genomic.fna//' <(echo $o))
o=$(sed 's/GCF_.*.1_//' <(echo $genome))
o=$(sed 's/GCF_.*.2_//' <(echo $o))
#o=$(sed 's/GCF_.*.3_//' <(echo $o))
echo $genome
echo $o
sed -i "s/>/>DECOY_At-SOIL_"$o"_/g" "$genome"_genomic.fna
done
#######CONCATENATE GENOMES: to do BWA later on#######
#######We need the reference genomes, for that, we need to concatnate them#######
#######The double (>>) means that if does not rewrite it each time, but it still follws#######
for i in *
do
echo $i
cat $i >> REFERENCE.fna
done
#######The following thow should also work to concat#######
cat /ebio/abt6_projects9/Pseudomonas_diversity/data/pseudomonas_genome_assembly/2016_07_TGAC/all_assemblies/*.fasta > talia.fasta
cat $(ls /ebio/abt6_projects9/Pseudomonas_diversity/data/pseudomonas_genome_assembly/2016_07_TGAC/all_assemblies/*.fasta) > talia.fasta
#######Compress a file, while saving the previous version#######
gzip -k

#######BWA programme application (MAPPING TOOL): Mapping reads with genome references#######
#######First, we have to create the index!#######
bwa index #reference genome
#######Load trimmed raw reads after using SKEWER programme#######
cd /ebio/abt6_projects9/microbiome_analysis/data/processed_DATA/HiSeq3000-Run0104/skewer

reference=/ebio/abt6_projects9/microbiome_analysis/data/processed_DATA/HiSeq3000-Run0103/Psy_genomes/reference_genome/REFERENCE.fna.gz
#reference=/ebio/abt6_projects9/microbiome_analysis/data/processed_DATA/HiSeq3000-Run0104/Sph_genomes/reference_genome/REFERENCE.fna.gz
output_dir=/ebio/abt6_projects9/microbiome_analysis/data/processed_DATA/HiSeq3000-Run0103/bwa_mapping
#output_dir=/ebio/abt6_projects9/microbiome_analysis/data/processed_DATA/HiSeq3000-Run0104/bwa_mapping

for o in S*_skewered-trimmed-pair1.fastq.gz
do
read1=$o
read2=$(sed 's/pair1/pair2/' <(echo $o))
S=$(sed 's/_skewered-trimmed-pair..fastq.gz//' <(echo $o))
S=$(sed 's/S//' <(echo $S))
echo $read1
echo $read2
echo $S
qsub -N J"$S" -o "$output_dir"/out.txt -e "$output_dir"/err.txt -V -cwd /ebio/abt6_projects7/small_projects/rjove/map_to_reference.sh "$read1" "$read2" "$S" "$output_dir" "$reference"
#bwa mem -t 10 "$reference" "$read1" "$read2" | ~/software/samtools-1.3.1/samtools view -q 30 -b -@ 10 - | ~/software/samtools-1.3.1/samtools  sort -@ 10 - > "$output_dir"/S"$sample_number".bam
done

#######Mapping reads with decoy genomes#######
#vi /ebio/abt6_projects7/bacterial_strain_analysis/code/bacterial_reference_map.sh
#reference=/ebio/abt6_projects9/microbiome_analysis/data/processed_DATA/HiSeq3000-Run0103/Pseudo_and_Decoy_REF/PseudoREF_LocalRefseqDecoy.fasta
#readdir=/ebio/abt6_projects9/microbiome_analysis/data/processed_DATA/HiSeq3000-Run0103/skewered
#for r in "$readdir"/*pair1.fastq.gz
#do
#read1=$r
#read2=$(sed 's/pair1/pair2/' <(echo $r))
#S=$(sed 's/.*skewered\/S//' <(echo $r))
#S=$(sed 's/_skewered-trimmed-pair1.fastq.gz//' <(echo $S))
#qsub -N job"$S" -o $(pwd)/out.txt -e $(pwd)/err.txt -V -cwd bacterial_reference_map.sh $read1 $read2 $reference $S
#done

#######PICARD programme application to remove PCR duplicates#######
#######SOLVED: it seems we have a limited amount of temporal files that can be created unlimitedly -n#######
#######In out case is 1024 in length, the MarkDuplicates has an option MAX_FILE_HANDLES=Integer defined as 1000#######
mkdir noduplicates
out_dir=/ebio/abt6_projects9/microbiome_analysis/data/processed_DATA/HiSeq3000-Run0103/bwa_mapping/noduplicates

for i in S*bam
do
i=$(sed 's/S//' <(echo $i))
sample_number=$(sed 's/.bam//' <(echo $i))
echo $sample_number
INPUT=S"$sample_number"
OUTPUT=S"$sample_number"_nodups
echo $INPUT
echo $OUTPUT
qsub -N J"$sample_number" -o "$out_dir"/out.txt -e "$out_dir"/err.txt -V -cwd /ebio/abt6_projects7/bacterial_strain_analysis/code/remove_duplicatesRoger.sh "$INPUT" "$OUTPUT" "$sample_number" "$out_dir"
done

#/ebio/abt6_projects9/microbiome_analysis/data/software/jre1.8.0_73/bin/java -jar /ebio/abt6_projects9/microbiome_analysis/data/software/picard-tools-2.2.1/picard.jar MarkDuplicates I="$INPUT".bam O="$OUTPUT".bam M=metrics"$sample_number".txt REMOVE_DUPLICATES=true MAX_FILE_HANDLES=1000

#######STATISTIC CALCULATION of each BAM file: percentage of the amount of reads that are covered across the genome references, and their quality#######
#for i in S*bam
#do
#i=$(sed 's/S//' <(echo $i))
#sample_number=$(sed 's/.bam//' <(echo $i))
#echo $sample_number
#samtools stats S"$sample_number".bam > ./stats/S"$sample_number"_stats.txt
#done

#######STATISTIC & UNIQUE COUNT#######
#######Perform the BAMSTATS together with the unique table that counts the number of mapped reads to each of the different input genomes:#######
for j in *.bam
do
j=$(sed 's/.bam//' <(echo $j))
echo $j
qsub -N job"$j" -o $(pwd)/out_2.txt -e $(pwd)/err_2.txt -V -cwd /ebio/abt6_projects7/bacterial_strain_analysis/code/samtools_stats.sh $j
#samtools view -F4 "$S".bam | cut -f3 | sort | uniq -c > ./mapping_table/"$S"unique.txt
done

#######ANALYSE STATS and UNIQUE COUNT: Create summary table stats#######
#######Create STATS directory#######
mkdir summary

#######First, take the useful infomation from the file together with the name of the file#######
for i in *.txt
do
i=$(sed 's/_stats.txt//' <(echo $i))
echo "$i"
sed -n '/# The command/{s/# The .*S/S/; p}' "$i"_stats.txt >./summary/"$i"_sum.txt
grep ^SN "$i"_stats.txt | cut -f 2- >> ./summary/"$i"_sum.txt
done

#######Before taking the second column (numbers) of each sample, we take the first coloumn with the names#######
cd ./summary
mkdir full_table
cut -f 1 S10_clean_sum.txt > ./full_table/col_labels.txt

#######We take the second column#######
mkdir ./full_table/col_2
for i in *.txt
do
i=$(sed 's/_sum.txt//' <(echo $i))
echo "$i"
cut -f 2 "$i"_sum.txt > ./full_table/col_2/"$i"_col2.txt
done

#######We join all second columns#######
cd ./full_table/col_2
paste *col2.txt > ../sum_table.txt

#######We combine both files containing labels and sum_table#######
cd ..
paste col_labels.txt sum_table.txt > full_summary.txt

#######Now, for the UNIQUE lists of contigs and the corresponding MAPPED reads#######
#######To change the column order and separate with a tab#######
for i in *.txt
do
i=$(sed 's/_cleanunique.txt//' <(echo $i))
echo "$i"
awk '{ print $2 " " $1}' "$i"_cleanunique.txt | sed 's/ /\t/' - > "$i"_cleanunique_rearranged.txt
done

#######Now, we need to combine all values with REFSEQ#######
#######It might be sub-optimal, but it will be helpful#######
mkdir summary

for i in *_cleanunique_rearranged.txt
do
i=$(sed 's/_cleanunique_rearranged.txt//' <(echo $i))
echo "$i"
awk '$1 ~ /REFSEQ/ {sum += $2} END {print "REFSEQ""\t" sum}' "$i"_cleanunique_rearranged.txt > ./summary/"$i"summary_clean.txt
#awk '$1 ~ /JULIAN/ {sum += $2} END {print "JULIAN""\t" sum}' "$i"_cleanunique_rearranged.txt >> ./summary/"$i"summary_clean.txt
awk '$1 ~ /DECOY/ {sum += $2} END {print "DECOY""\t" sum}' "$i"_cleanunique_rearranged.txt >> ./summary/"$i"summary_clean.txt
awk '$1 ~ /nonOTU5_LOCAL165/ {sum += $2} END {print "nonOTU5_LOCAL165""\t" sum}' "$i"_cleanunique_rearranged.txt >> ./summary/"$i"summary_clean.txt
awk '$1 ~ /^OTU5_LOCAL165/ {sum += $2} END {print "OTU5_LOCAL165""\t" sum}' "$i"_cleanunique_rearranged.txt >> ./summary/"$i"summary_clean.txt
awk '{s+=$2}END{print "TOTAL""\t" s}' "$i"_cleanunique_rearranged.txt >> ./summary/"$i"summary_clean.txt
#######Again, be careful with the ' or " and, in this ocassion, we have to add another \ !!, otherwise it cannot read the command#######
sed -i "1 i\\"$i"_summarytable" ./summary/"$i"summary_clean.txt
done

#######Before taking the second columm (numbers) of each sample, we take the first column with the names#######
cd ./summary
mkdir full_table
cut -f 1 S10summary_clean.txt > ./full_table/col_labels.txt
#######We take the second column#######
mkdir ./full_table/col_2
for i in *.txt
do
i=$(sed 's/summary_clean.txt//' <(echo $i))
echo "$i"
cut -f 2 "$i"summary_clean.txt > ./full_table/col_2/"$i"_col2.txt
done
#######We join all second columns#######
cd ./full_table/col_2
paste *col2.txt > ../sum_table.txt
#######We combine both files (labels and summary_table)
cd ..
paste col_labels.txt sum_table.txt > full_summary_uniquemap.txt

#######REMOVE FLAGGED READS#######
#######So, we had more reads than expected, it is because of multimapping reads, to solve it, we will remove those flagged reads#######
#######To see the multiple mapping reads -> print a list with the read IDs and the count per read#######
samtools view -F4 S10_nodups.bam | cut -f1 | sort | uniq -c
#this command line will save all those IDs of the reads that show more than 2 reads
samtools view -F4 S10_nodups.bam | cut -f1 | sort | uniq -c | awk '$1>2{print$2}' > S10list.txt
#this command is the same but without the storage and to check if really grep the ones with more than 2
samtools view -F4 S10_nodups.bam | cut -f1 | sort | uniq -c | awk '$1>2{print$2"\t"$1}'| wc -l
#once we have the list I wanted to check were do they mapped, I grep the files in the list and count for unique hits:
samtools view S10nomulti.bam | grep -f S10list.txt | cut -f 3 | sort | uniq -c > mapping_multi.txt

#######Now, the real deal comes#######
out_dir=/ebio/abt6_projects9/microbiome_analysis/data/processed_DATA/HiSeq3000-Run0104/bwa_mapping/noduplicates/clean_bam
for q in *.bam
do
q=$(sed 's/_nodups.bam//' <(echo $q))
echo $q
#######In the shell script, the following three steps can be achieved in a single command step !! if wanted#######
qsub -N J"$q" -o "$out_dir"/out3.txt -e "$out_dir"/err3.txt -V -cwd /ebio/abt6_projects7/bacterial_strain_analysis/code/samtools_clean_roger.sh "$q" "$out_dir"
done
#samtools view -h -F 4 -F 256 -F 2048 -bS -o S10final.bam S10_nodups.bam
#######Here we might want to do also STATS and Summary_table, go above again#######

#######FINAL BAM file w/o DECOY: only the good mappings#######
#######This one will give back the mapped removing DECOY ones ! and the output is also in BAM file#######
out_dir=/ebio/abt6_projects9/microbiome_analysis/data/processed_DATA/HiSeq3000-Run0104/bwa_mapping/noduplicates/clean_bam
for q in *.bam
do
q=$(sed 's/_clean.bam//' <(echo $q))
echo $q
#######Be careful with the shell script to activate/inactivate whatever you need#######
qsub -N J"$q" -o "$out_dir"/out.txt -e "$out_dir"/err.txt -V -cwd /ebio/abt6_projects7/bacterial_strain_analysis/code/samtools_clean_roger.sh "$q" "$out_dir"
done
#######This one is the correct#######
#samtools view -h S10_clean.bam | awk '$3 !~ /DECOY/' | samtools view -h -bS -o S10test.bam -
#######Remove more reads than the supposed to be#######
#samtools view -h -F 4 S10_nodups.bam | grep -v "DECOY" | samtools view -bS -o S10test.bam -

#######Samtools programme application: transforming BAM file to FASTQ file#######
out_dir=/ebio/abt6_projects9/microbiome_analysis/data/processed_DATA/HiSeq3000-Run0104/bwa_mapping/noduplicates/clean_bam
for q in *.bam
do
q=$(sed 's/_final.bam//' <(echo $q))
echo $q
qsub -N J"$q" -o "$out_dir"/out.txt -e "$out_dir"/err.txt -V -cwd /ebio/abt6_projects7/bacterial_strain_analysis/code/samtools_clean_roger.sh "$q" "$out_dir"
#samtools sort -n S10_final.bam | samtools fastq - > S10.fastq.gz
done

#######The FASTQ files are goinng to be used for MASH application: LOOK at Mash_analysis.sh
#######OTHER analyses that we can do with the CLEAN data:#######
#######We want to extract the total reads for genome from the UNIQUE files, where also the contigs are included, unfortunately#######
#######The original (below), was to be RESTRICTIVE with the name, because afterwards there was the _contig....#######
#######However, if we use the newer verion "${line}", different genomes were combined and what we did was to change the list files#######
#######He added an _ after the genome name as found in the reads, so then the problem was solved#######
while read line; do grep "\<${line}\>" testfile.txt | awk '{sum+=$2} END {print sum}' > $line.count.txt ; done < pattern.txt

#######This works well, but a lot of files were created#######
#######While read line#######
#######DO this ONE#######
#grep "${line}" S10_cleanunique_rearranged.txt | awk '{sum+=$2} END {print $1"\t"sum}' > $line.count.txt
#done < /ebio/abt6_projects9/microbiome_analysis/data/processed_DATA/HiSeq3000-Run0103/Psy_genomes/genomes_list/all.txt
#for Sph --> /ebio/abt6_projects9/microbiome_analysis/data/processed_DATA/HiSeq3000-Run0104/Sph_genomes/genomes_list/all.txt
#######Loop for all these files#######
for l in *_cleanunique_rearranged.txt
do
l=$(sed 's/_cleanunique_rearranged.txt//' <(echo $l))
while read line
do
grep "${line}" "$l"_cleanunique_rearranged.txt | awk '{sum+=$2} END {print $1"\t"sum}' >> "$l".count.txt
done < /ebio/abt6_projects9/microbiome_analysis/data/processed_DATA/HiSeq3000-Run0104/Sph_genomes/genomes_list/all.txt
sed -i "1 i\\"$l"_summarytable" "$l".count.txt
echo "$l"
done

#######Process the data#######
mkdir column2
for d in *count.txt
do
d=$(sed 's/.count.txt//' <(echo $d))
echo "$d"
cat "$d".count.txt | cut -f2 > ./column2/"$d"_column2
done

#######Combine the files of column2#######
cd column2
cp /ebio/abt6_projects9/microbiome_analysis/data/processed_DATA/HiSeq3000-Run0104/Sph_genomes/genomes_list/all.txt labels.txt
sed -i '1 i\fulltable' labels.txt # this adds to the first line the name fulltable so then when combining the files, the table has the right amount of rows

#######Paste $(cat all) # this mark makes no sense, since we added the first row and ther is not such a file with this name#######
paste *column2 > fullmatrix
paste labels.txt fullmatrix > genomecount_sample.txt # combines the name column and the numbers column

#######IF we want to go from number of reads per genome to present/absent reads per genome:#######
#######Delete the Decoy rows#######
awk '!/DECOY/' genomecount_sample.txt > genomecount_samplenodecoy.txt
#######The percentage and other paremeters determining quality of the READS will be done in R and then saved#######
#######Transform the values to 0 or 1 (above 50 reads, or 0.002 read fraction) for further analyses#######
awk 'BEGIN { FS = OFS = "\t" } { for(l=1; l<=NF; l++) if($l<0.001) $l = 0 } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = 0 } {if (NR>1) {for (C=2; C<=NF; C++) {if ($C >= 0.001) {$C=1}}}print}'
#Sort the file according to the S#######
perl -pale '
   $. == 1 and
   @I = map  { $_->[1] }
        sort { $a->[0] <=> $b->[0] }
        map  { [ $F[$_] =~ /^S(\d+)$/, $_ ] } 1..$#F;
   $_ = "@F[0, @I]";
'
#######Then, we have to change the spaces for tab in textwrangles#######
#good_Sample.psy <- c("S1","S2","S6","S8","S9","S10","S11","S12","S13","S14","S15","S16","S17","S18","S19","S20","S21","S23","S25","S26","S27","S29","S30","S31","S32","S34","S35","S39","S43","S44","S46","S47","S48","S49","S52","S53","S54","S55","S56")
#all except 3, 4, 5, 7, 22, 24, 28, 33, 36, 37, 38, 40, 41, 42, 45, 50, 51
#######Then, we want to remove certain columns#######
#######Be careful! We have to add one number to each column because we have the first column name#######
#=$3=$4=$5=$7=$22=$24=$28=$33=$36=$37=$38=$40=$41=$42=$45=$50=$51
#######the one below cut, not working becasue there are 0...#######
#cut -f 1, 2, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 25, 26, 27, 29, 30, 31, 32, 34, 35, 39, 43, 44, 46, 47, 48, 49, 52, 53, 54, 55, 56
awk 'BEGIN {IFS="\t"} {OFS="\t"} {$4=$5=$6=$8=$23=$25=$29=$34=$37=$38=$39=$41=$42=$43=$44=$51=$52="\s"; print $0}' file
#######Afterwards remove the s\t in textwrangler#######
#######Then, I was tired and did the IF command in excel to summarise all the samples from one species#######

#######To extract specific reads mapping to locations on the genomes#######
samtools view -H S18_final.bam > header.sam
samtools view S18_final.bam | awk '$3 ~ /^OTU5_LOCAL165_plate13.E2/' | cat header.sam - | samtools view -Sb - > test.bam


#######OTHER USEFUL COMMANDS#######
#######Plot BAMSTATS#######
#plot-bamstats -p "53" S53_stats.txt
samtools stats S"$sample_number".bam > ./stats/S"$sample_number"_stats.txt
cd stats
/ebio/abt6_projects9/microbiome_analysis/data/software/samtools-1.3.1/bin/plot-bamstats -p "p" S"$i""$filename"_stats.txt
echo "done"

#######Sum second column#######
awk '{s+=$2}END{print s}' summary.txt

#######Take out unmapped reads and make them as FASTA input for BLAST#######
samtools view -f 4 S1_nodups.bam | awk '{print ">"$1"\n" $10}' > unmappedS1.fasta
#######BAM file to FASTA file transformation#######
#######This one is in FASTA format, but for SKETCH, we want BAM file, not FASTA file#######
samtools view -F4 -b S10_nodups.bam | awk '!/DECOY/ {OFS="\t"; print ">"$1"_"$3"\n"$10}' > S10_noDEC.fasta
#######Count READS#######
#######To check the number of reads and it has to be the same as the supposed to be mapped reads!!#######
samtools view S10_nodups.bam | grep '^ST'| wc -l
samtools view -c
