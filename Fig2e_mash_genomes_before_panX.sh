mash=/ebio/abt6_projects9/microbiome_analysis/data/software/mash/mash

#Mash
mash_directory=/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/MASH/DEREK_MASH_GENOME
cd $mash_directory
mkdir -p sketch 
#now get all input files and run sketch module from Mash programme
for input in p*
do
nice -10 $mash sketch -r -o sketch_"$input" "$input"
done

#then get all msh files as an output and make lists from all msh files
ls *msh > list_of_all_sketches.txt
#then run all msh files with list_of_all_sketches.txt to get values from all genomes
for dist in *msh
do
  echo $dist
  $mash dist -t $dist $(cat list_of_all_sketches.txt) > distmatrix_$dist.txt
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





