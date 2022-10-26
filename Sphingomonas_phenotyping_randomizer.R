Sphingomonas_strains=c(1:21)
Pseudomonas_strains=c("DC3000", "OTU5", "buffer")

#3 Pseudomonas treatments each with 336 plants
#Each treatment will be done in 2 blocks of 7 plates, with 8 replicates of each Sphingomonas strain
# per block (total of 16 replicates of each Sphingomonas overall)
#This is 6 blocks of Sphingomonas in total


Sph_block1=sample(rep(Sphingomonas_strains, 8), length(rep(Sphingomonas_strains, 8)))

Sph_block1=matrix(Sph_block1, ncol=6, byrow=TRUE)
blank_space=matrix(rep(rep("", 6),3), ncol=6, byrow=TRUE)
Sph_block1=rbind(Sph_block1[1:4,], blank_space, Sph_block1[5:8,],blank_space, Sph_block1[9:12,],blank_space, Sph_block1[13:16,],blank_space, Sph_block1[17:20,], blank_space, Sph_block1[21:24,],blank_space, Sph_block1[25:28,])

write.table(Sph_block1, "/Users/dlundberg/Documents/abt6/Pratchaya/random_design_tables/Sph_block6.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

