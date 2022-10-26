mash=/ebio/abt6_projects9/microbiome_analysis/data/software/mash/mash

#Mash
NR_directory=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Sphingomonas/NonRedundant
mkdir -p $NR_directory


#>
	#>
		#>
			#> now get all input files and run sketch module from Mash programme.  Write all the sketches to the same place.
		#>
	#>
#>
#REFSEQ_assembly
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Sphingomonas/RefSEQ_genome
for input in GCF*fna
do
echo $input
nice -10 $mash sketch -s 5000 -r -o $NR_directory/sketch_"$input" "$input"
done
#mine
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Sphingomonas/Local_Sphingomonas
for input in S*H*pilon.fasta
do
echo $input
nice -10 $mash sketch -s 5000 -r -o $NR_directory/sketch_"$input" "$input"
done
#Sphingomonas in L-Sphere
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-LSPHERE/Sphingo_LSPHERE
for input in GCF*fna
do
echo $input
nice -10 $mash sketch -s 5000 -r -o $NR_directory/sketch_"$input" "$input"
done
#Sphingomonas in R-Sphere
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-ROOT/Sphingo_ROOT
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
# "remove_redundant_genomes.R", set to "Sphingomonas"


#SCP THE CONDENSED LIST BACK TO BURRITO
scp /Users/dlundberg/Documents/abt6/Pratchaya/Bulk_metagenomes/NonRedundantMapping/Sphingomonas/condensed_Sref.txt dlundberg@burrito.eb.local:/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Sphingomonas/NonRedundant/


#>
	#>
		#>
			#> #MOVE THE CONDENSED SET OF GENOMES TO A SINGLE FOLDER OF ASSEMBLIES
			#> #AND ADD DECOYS TO THE LIST.
		#>
	#>
#>
NR_directory=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Sphingomonas/NonRedundant
cd $NR_directory 
ls /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-ROOT/nonPseudo_nonSphingo_ROOT > nonPseudo_nonSphingo_ROOT.txt
ls /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-LSPHERE/nonPseudo_nonSphingo_LSPHERE > nonPseudo_nonSphingo_LSPHERE.txt
ls /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-SOIL/nonPseudo_nonSphingo_SOIL > nonPseudo_nonSphingo_SOIL.txt

cat nonPseudo_nonSphingo_ROOT.txt >> condensed_Sref.txt
cat nonPseudo_nonSphingo_LSPHERE.txt >> condensed_Sref.txt
cat nonPseudo_nonSphingo_SOIL.txt >> condensed_Sref.txt
cat condensed_Sref.txt | sort | uniq > Sphingomonas_BulkMap_Index.txt


Sphingo_NR_REF=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Sphingomonas/NonRedundant/Sphingo_NR_REF
mkdir -p $Sphingo_NR_REF
cd $Sphingo_NR_REF


for genome in $(cat $NR_directory/Sphingomonas_BulkMap_Index.txt)
do
echo $genome
    cp /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Sphingomonas/RefSEQ_genome/"$genome" $Sphingo_NR_REF    #128
    cp /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Sphingomonas/Local_Sphingomonas/"$genome" $Sphingo_NR_REF   #148
    cp /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-LSPHERE/Sphingo_LSPHERE/"$genome" $Sphingo_NR_REF   #23
    cp /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-ROOT/Sphingo_ROOT/"$genome" $Sphingo_NR_REF    #3
    cp /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-ROOT/nonPseudo_nonSphingo_ROOT/"$genome" $Sphingo_NR_REF   #180
    cp /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-LSPHERE/nonPseudo_nonSphingo_LSPHERE/"$genome" $Sphingo_NR_REF   #167
    cp /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Decoy/At-SOIL/nonPseudo_nonSphingo_SOIL/"$genome" $Sphingo_NR_REF
done



#Do the same for interesting genes
NR_directory=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Sphingomonas/NonRedundant
cd $NR_directory 
Sphingo_GENES_REF=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Sphingomonas/NonRedundant/Sphingo_GENES_REF
mkdir -p $Sphingo_GENES_REF
cd $Sphingo_GENES_REF




#>
	#>
		#>
			#> #RENAME ALL CONTIG HEADERS AFTER THE GENOME
		#>
	#>
#>
cd $Sphingo_NR_REF
for o in *.fasta *.fna
do
    genome_name=$(sed "s/_pilon.*fasta//g" <(echo $o))
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
cat $i >> SPH_NR_REFERENCE.fasta2
done

rm *.fasta *.fna
cat SPH_NR_REFERENCE.fasta2 > SPH_NR_REFERENCE.fasta
rm SPH_NR_REFERENCE.fasta2
#>
	#>
		#>
			#> #MAKE INDEX
		#>
	#>
#>
fasta=SPH_NR_REFERENCE.fasta
qsub -N "$fasta" -o $(pwd)/out.txt -e $(pwd)/err.txt -V -cwd /ebio/abt6_projects7/bacterial_strain_analysis/code/bwa_index.sh "$fasta"




#NOW RUN PRATCHAYAS MAPPING SCRIPT FOR TWO SPHINGOMONAS OUTPUT DIRECTORIES. 
#SPRING
outputDir=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/process_scrape_metagenomes/Sphingomonas_Spring
path_to_REF=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Sphingomonas/NonRedundant/Sphingo_NR_REF/SPH_NR_REFERENCE.fasta
readsDir=NA
cd $outputDir
#Iterate through the "skewered" directory to get all the sample numbers
for sample_number in $(ls /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/process_scrape_metagenomes/Sphingomonas_Spring/skewered | grep "S.*_" | sed "s/_.*//g" | sed "s/S/ /g" | sort | uniq )
do 
    echo $sample_number
    qsub -N Ssp"$sample_number" -o $(pwd)/out.txt -e $(pwd)/err.txt -V -cwd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/scripts/skewer_mapping_dedup_SphPsySummer.sh "$sample_number" "$readsDir" "$outputDir" "$path_to_REF"
done

#SUMMER
outputDir=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/process_scrape_metagenomes/Sphingomonas_Summer
path_to_REF=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/mapping_reference_genomes/Sphingomonas/NonRedundant/Sphingo_NR_REF/SPH_NR_REFERENCE.fasta
readsDir=NA
cd $outputDir
#Iterate through the "skewered" directory to get all the sample numbers
for sample_number in $(ls /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/process_scrape_metagenomes/Sphingomonas_Summer/skewered | grep "S.*_" | sed "s/_.*//g" | sed "s/S/ /g" | sort | uniq )
do 
    echo $sample_number
    qsub -N Ssu"$sample_number" -o $(pwd)/out.txt -e $(pwd)/err.txt -V -cwd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/scripts/skewer_mapping_dedup_SphPsySummer.sh "$sample_number" "$readsDir" "$outputDir" "$path_to_REF"
done





scp dlundberg@burrito.eb.local:/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/process_scrape_metagenomes/Sphingomonas_Spring/mapping_table/Sphingomonas_Spring_full_mapping_table.txt /Users/dlundberg/Documents/abt6/Pratchaya/Bulk_metagenomes/

