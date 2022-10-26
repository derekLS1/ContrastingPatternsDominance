mash=/ebio/abt6_projects9/microbiome_analysis/data/software/mash/mash

#Mash
NR_directory=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Pseudomonas/NonRedundant
mkdir -p $NR_directory


#>
	#>
		#>
			#> now get all input files and run sketch module from Mash programme.  Write all the sketches to the same place.
		#>
	#>
#>
#TK165_assembly
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Pseudomonas/TK165
for input in p*fasta
do
echo $input
nice -10 $mash sketch -s 5000 -r -o $NR_directory/sketch_"$input" "$input"
done
#REFSEQ_assembly
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Pseudomonas/RefSEQ_genomes
for input in GCF*fna
do
echo $input
nice -10 $mash sketch -s 5000 -r -o $NR_directory/sketch_"$input" "$input"
done
#mine_nonOTU5
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Pseudomonas/Psy_DL/nonOTU5_pilon
for input in Psy*fasta
do
echo $input
nice -10 $mash sketch -s 5000 -r -o $NR_directory/sketch_"$input" "$input"
done
#mine_OTU5
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Pseudomonas/Psy_DL/OTU5_pilon
for input in Psy*fasta
do
echo $input
nice -10 $mash sketch -s 5000 -r -o $NR_directory/sketch_"$input" "$input"
done
#Pseudomonas in L-Sphere
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-LSPHERE/Pseudo_LSPHERE
for input in GCF*fna
do
echo $input
nice -10 $mash sketch -s 5000 -r -o $NR_directory/sketch_"$input" "$input"
done
#Pseudomonas in R-Sphere
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-ROOT/Pseudo_ROOT
for input in GCF*fna
do
echo $input
nice -10 $mash sketch -s 5000 -r -o $NR_directory/sketch_"$input" "$input"
done


#then get all msh files as an output and make a list of all msh files
cd $NR_directory
ls *msh > list_of_all_sketches.txt
#then run all msh files with list_of_all_sketches.txt to get values from all genomes
for dist in *msh
do
  echo $dist
  nice -10 $mash dist -t $dist $(cat list_of_all_sketches.txt) > distmatrix_$dist.txt
done


#now copy all values in column2
for d in distmatrix*
do
echo "$d"
cat "$d" | cut -f2 > "$d"_column2
done
paste *column2 > fullmatrix.txt
#now clean unncessary files and mv all msh files into sketches folder
rm *column2
rm distmatrix*
mkdir -p sketches
mv list_of_all_sketches.txt sketches
mv sketch_*.msh sketches


#SCP THE FULLMATRIX RESULT TO MAC
scp dlundberg@burrito.eb.local:/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Pseudomonas/NonRedundant/fullmatrix.txt /Users/dlundberg/Documents/abt6/Pratchaya/Bulk_metagenomes/NonRedundantMapping/Pseudomonas

#RUN THE R SCRIPT TO PARE DOWN BASED ON ONLY GENOMES AT LEAST 1% DIFFERENT FROM ALL OTHERS
# "remove_redundant_genomes.R", set to "Pseudomonas"

#SCP THE CONDENSED LIST BACK TO BURRITO
scp /Users/dlundberg/Documents/abt6/Pratchaya/Bulk_metagenomes/NonRedundantMapping/Pseudomonas/condensed_Pref.txt dlundberg@burrito.eb.local:/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Pseudomonas/NonRedundant/


#>
	#>
		#>
			#> #MOVE THE CONDENSED SET OF GENOMES TO A SINGLE FOLDER OF ASSEMBLIES
		#>
	#>
#>
NR_directory=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Pseudomonas/NonRedundant
cd $NR_directory 
ls /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-ROOT/nonPseudo_nonSphingo_ROOT > nonPseudo_nonSphingo_ROOT.txt
ls /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-LSPHERE/nonPseudo_nonSphingo_LSPHERE > nonPseudo_nonSphingo_LSPHERE.txt
ls /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-SOIL/nonPseudo_nonSphingo_SOIL > nonPseudo_nonSphingo_SOIL.txt

cat nonPseudo_nonSphingo_ROOT.txt >> condensed_Pref.txt
cat nonPseudo_nonSphingo_LSPHERE.txt >> condensed_Pref.txt
cat nonPseudo_nonSphingo_SOIL.txt >> condensed_Pref.txt
cat condensed_Pref.txt | sort | uniq > Pseudomonas_BulkMap_Index.txt

Pseudo_NR_REF=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Pseudomonas/NonRedundant/Pseudo_NR_REF
mkdir -p $Pseudo_NR_REF
cd $Pseudo_NR_REF

#TK165_assembly
for genome in $(cat $NR_directory/condensed_Pref.txt)
do
echo $genome
    cp /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Pseudomonas/TK165/"$genome" $Pseudo_NR_REF
    cp /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Pseudomonas/RefSEQ_genomes/"$genome" $Pseudo_NR_REF
    cp /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Pseudomonas/Psy_DL/nonOTU5_pilon/"$genome" $Pseudo_NR_REF
    cp /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Pseudomonas/Psy_DL/OTU5_pilon/"$genome" $Pseudo_NR_REF
    cp /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-LSPHERE/Pseudo_LSPHERE/"$genome" $Pseudo_NR_REF
    cp /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-ROOT/Pseudo_ROOT/"$genome" $Pseudo_NR_REF
	cp /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-ROOT/nonPseudo_nonSphingo_ROOT/"$genome" $Pseudo_NR_REF
    cp /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-LSPHERE/nonPseudo_nonSphingo_LSPHERE/"$genome" $Pseudo_NR_REF
    cp /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-SOIL/nonPseudo_nonSphingo_SOIL/"$genome" $Pseudo_NR_REF
done


###DO THE SAME FOR SPECIFIC GENE CLASSES



#>
	#>
		#>
			#> #RENAME ALL CONTIG HEADERS AFTER THE GENOME
		#>
	#>
#>
Pseudo_NR_REF=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Pseudomonas/NonRedundant/Pseudo_NR_REF
cd $Pseudo_NR_REF
rename "s/_L8_pilon.fasta/H133_pilon.fasta/g" *
rename "s/PsyRun133_//g" *

for o in *.fasta *.fna
do
    genome_name=$(sed "s/_pilon.*fasta//g" <(echo $o))
    genome_name=$(sed "s/.fasta//g" <(echo $genome_name))
    genome_name=$(sed "s/GCF_/GCF/g" <(echo $genome_name))
    genome_name=$(sed "s/_.*//g" <(echo $genome_name))
    genome_name=$(sed "s/GCF/GCF_/g" <(echo $genome_name))
    echo $genome_name
	sed -i "s/>/>"$genome_name"_/g" "$o"
done


#>
	#>
		#>
			#> #CONCATENATE GENOMES: to do BWA later on. Delete the individual ones afterwards to save space.
		#>
	#>
#>
for i in *.fasta *.fna
do
echo $i
cat $i >> PSE_NR_REFERENCE.fasta2
done

rm *.fasta *.fna
cat PSE_NR_REFERENCE.fasta2 > PSE_NR_REFERENCE.fasta
rm PSE_NR_REFERENCE.fasta2
#>
	#>
		#>
			#> #MAKE INDEX
		#>
	#>
#>
fasta=PSE_NR_REFERENCE.fasta
qsub -N "$fasta" -o $(pwd)/out.txt -e $(pwd)/err.txt -V -cwd /ebio/abt6_projects7/bacterial_strain_analysis/code/bwa_index.sh "$fasta"




#NOW RUN PRATCHAYAS MAPPING SCRIPT FOR TWO SPHINGOMONAS OUTPUT DIRECTORIES. 
#SPRING
outputDir=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/process_scrape_metagenomes/Pseudomonas_Spring
path_to_REF=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Pseudomonas/NonRedundant/Pseudo_NR_REF/PSE_NR_REFERENCE.fasta
readsDir=NA
cd $outputDir
#Iterate through the "skewered" directory to get all the sample numbers
for sample_number in $(ls /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/process_scrape_metagenomes/Pseudomonas_Spring/skewered | grep "S.*_" | sed "s/_.*//g" | sed "s/S/ /g" | sort | uniq )
do 
    echo $sample_number
    qsub -N Psp"$sample_number" -o $(pwd)/out.txt -e $(pwd)/err.txt -V -cwd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/scripts/skewer_mapping_dedup_SphPsySummer.sh "$sample_number" "$readsDir" "$outputDir" "$path_to_REF"
done

#SUMMER
outputDir=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/process_scrape_metagenomes/Pseudomonas_Summer
path_to_REF=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Pseudomonas/NonRedundant/Pseudo_NR_REF/PSE_NR_REFERENCE.fasta
readsDir=NA
cd $outputDir
#Iterate through the "skewered" directory to get all the sample numbers
for sample_number in $(ls /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/process_scrape_metagenomes/Pseudomonas_Summer/skewered | grep "S.*_" | sed "s/_.*//g" | sed "s/S/ /g" | sort | uniq )
do 
    echo $sample_number
    qsub -N Psu"$sample_number" -o $(pwd)/out.txt -e $(pwd)/err.txt -V -cwd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/scripts/skewer_mapping_dedup_SphPsySummer.sh "$sample_number" "$readsDir" "$outputDir" "$path_to_REF"
done

#>
	#>
		#>
			#> #MAKE TABLE
		#>
	#>
#>




