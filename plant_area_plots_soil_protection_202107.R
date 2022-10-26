library("ggplot2")
library("reshape2")

PA<-read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/soil_protection_photos_202107/size_data.txt") 
	as.matrix(PA[1,])->colnames(PA)
	as.matrix(PA[,1])->PA[,1]
	PA[2:nrow(PA),]->PA
	
	
	sort(as.matrix(PA$plantID))
	
	
	factor(PA$sphingo, levels=c("M", "B", "F", "S"))->PA$sphingo
	factor(PA$pseudo, levels=c("none", "DC", "p25c2", "MgCl"))->PA$pseudo
	factor(PA$genotype, levels=c("WT", "eds", "coi"))->PA$genotype
	factor(PA$date, levels=c("20210702", "20210707", "20210712"))->PA$date
	as.numeric(as.matrix(PA$pixels))->PA$pixels
	as.matrix(PA$plantID)->PA$plantID
	
date2=format(Sys.Date(), format="%Y%m%d")

#REMOVE PLANTS WITH LESS THAN X PIXELS
PA[which(PA$pixels>=50),]->PA  

#keep only plants in all three pictures
PA[match(PA$plantID, names(which(table(as.matrix(PA$plantID))==3)), nomatch=0)>0,]->PA

#count
sort(table(as.matrix(PA$plantID)))

names(table(as.matrix(PA$plantID)))->plants_order




#copy PA for last day to PA_relative, and replace the pixels with % change of each plant 
#MAKES PA RELATIVE
PA[which(PA$date==20210712),]->PA_relative
for(p in plants_order){
	day0=PA$pixels[intersect(which(PA$date=="20210707"), which(PA$plantID==p))]
	day5=PA$pixels[intersect(which(PA$date=="20210712"), which(PA$plantID==p))]
	(day5-day0)/day0->PA_relative$pixels[which(PA_relative$plantID==p)]
}
	
	


options(scipen=999)


#PA->PA2
#PA[which(PA$date==20210712),]->PA2
PA_relative->PA2
#Chart to see pathogenic effect of Pseudomonas
#RAW SIZE CHART. 
PA2$pseudo=factor(PA2$pseudo,levels=c("MgCl", "DC", "p25c2"))
#convert pixel to mm^2
#mm_factor=9996/187389


write.table(PA2, "/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/soil_protection_photos_202107/Plant_size_table_for_David.txt")

pdf(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/soil_protection_photos_202107/Pseudo_by_genotype_", date, ".pdf", sep=""), width = 5, height = 6, useDingbats=FALSE)
		#print(
		ggplot(PA2, aes(fill=pseudo, y=pixels*100, x=genotype)) + 
   	 	geom_boxplot(outlier.size=0, outlier.shape=NA, alpha=0.75) + 
   	 	geom_point(aes(fill=pseudo), size = 1.2, stroke=0.1, shape = 21, position = position_jitterdodge(jitter.width=0.3)) +
   	 	theme_classic() + scale_fill_manual(values=c("#C576B0", "#6DE9FF", "#F7EC2A")) + 
    	scale_y_continuous(name ="Percent_growth", breaks = c(-100, 0, 100, 200, 300), labels = c(-100, 0, 100, 200, 300), limits = c(-100,300)) +
    	#scale_x_discrete(name ="Pseudomonas Strain") +
    	theme(axis.text=element_text(size=15, color="black"), axis.title=element_text(size=14,face="bold")) + 
    	theme(axis.text.x = element_text(angle=90, hjust=0.99, vjust=0.9)) + 
    	theme(axis.ticks = element_line(colour = 'black'))
    	#)
dev.off()


#Do the stats to compare pseudomonas roups


DC3000_coi=PA2[which(PA2$pseudo=="DC" & PA2$genotype=="coi") ,"pixels"]
DC3000_eds=PA2[which(PA2$pseudo=="DC" & PA2$genotype=="eds") ,"pixels"]
DC3000_WT=PA2[which(PA2$pseudo=="DC" & PA2$genotype=="WT") ,"pixels"]

	wilcox.test(DC3000_coi, DC3000_WT, "greater")
	#wilcox.test(DC3000_WT, DC3000_eds, "greater")
	wilcox.test(DC3000_coi, DC3000_eds, "greater")





Mg_coi=PA2[which(PA2$pseudo=="MgCl" & PA2$genotype=="coi") ,"pixels"]
p25_coi=PA2[which(PA2$pseudo=="p25c2" & PA2$genotype=="coi") ,"pixels"]
Mg_eds=PA2[which(PA2$pseudo=="MgCl" & PA2$genotype=="eds") ,"pixels"]
p25_eds=PA2[which(PA2$pseudo=="p25c2" & PA2$genotype=="eds") ,"pixels"]
Mg_WT=PA2[which(PA2$pseudo=="MgCl" & PA2$genotype=="WT") ,"pixels"]
p25_WT=PA2[which(PA2$pseudo=="p25c2" & PA2$genotype=="WT") ,"pixels"]

wilcox.test(Mg_coi, p25_coi, "greater")
wilcox.test(Mg_WT, p25_WT, "greater")
wilcox.test(Mg_eds, p25_eds, "greater")


p25_coi=PA2[which(PA2$pseudo=="p25c2" & PA2$genotype=="coi") ,"pixels"]
p25_eds=PA2[which(PA2$pseudo=="p25c2" & PA2$genotype=="eds") ,"pixels"]
p25_WT=PA2[which(PA2$pseudo=="p25c2" & PA2$genotype=="WT") ,"pixels"]

wilcox.test(p25_coi, p25_WT, "less")
wilcox.test(p25_coi, p25_eds, "less")
wilcox.test(p25_eds, p25_WT, "greater")


wilcox.test(p25_coi, p25_WT, "less")

p.adjust(c(

wilcox.test(DC3000_coi, DC3000_WT, "greater")[[3]],
wilcox.test(DC3000_coi, DC3000_eds, "greater")[[3]],
wilcox.test(Mg_WT, p25_WT, "greater")[[3]],
wilcox.test(Mg_eds, p25_eds, "greater")[[3]],
wilcox.test(Mg_coi, p25_coi, "greater")[[3]]), method = "fdr")



#magnitude of p25c2 effect stronger for coi1 

mean(Mg_coi)-mean(p25_coi)
mean(Mg_eds)-mean(p25_eds)
mean(Mg_WT)-mean(p25_WT)


DC3000_WT=PA2[which(PA2$pseudo=="DC" & PA2$genotype=="WT") ,"pixels"]

OTU5=PA2[which(PA2$Pseudomonas=="OTU5" & PA2$genotype=="eds") ,"pixels"]
BOIL=PA2[which(PA2$Pseudomonas=="BOIL" & PA2$genotype=="eds") ,"pixels"]
wilcox.test(OTU5, DC3000, "greater")
wilcox.test(BOIL, DC3000, "greater")
wilcox.test(OTU5, BOIL, "greater")




p.adjust(wilcox.test(FR1, MgCl, "greater")[[3]], method="fdr", n=5)






ggplot(PA, aes(fill=Pseudomonas, y=pixel_to_mm2*get(DATE), x=Sphingomonas)) + 
   	 	geom_boxplot(outlier.size=0, outlier.shape=NA, alpha=0.75) + 
   	 	geom_point(aes(fill=Pseudomonas), size = 1.2, stroke=0.1, shape = 21, position = position_jitterdodge(jitter.width=0.3)) +
   	 	theme_classic() + scale_fill_manual(values=c("#C576B0", "#6DE9FF", "#F7EC2A")) + 
    	scale_y_continuous(name ="Pixels", breaks = c(0, 100, 200, 300), labels = c(0, 100, 200, 300), limits = c(0,300)) +
    	scale_x_discrete(name ="Sphingomonas Strain", breaks = breaks, labels = plotting_key[,2], ) +
    	theme(axis.text=element_text(size=15, color="black"), axis.title=element_text(size=14,face="bold")) + 
    	theme(axis.text.x = element_text(angle=90, hjust=0.99, vjust=0.9) + #, color=plotting_key[,4])) + 
    	theme(axis.ticks = element_line(colour = 'black'))
    	)
	#dev.off()


##
	###
		####
			######
#To see growth promotion effect of Sphingomonas without any pseudomonas
#PA[c(which(PA$pseudo=="none"), which(PA$pseudo=="MgCl")),]->PAsphingo
#PA[c(which(PA$pseudo=="none"), which(PA$pseudo=="p25c2")),]->PAsphingo
#PA[c(which(PA$pseudo=="MgCl")),]->PAsphingo
PA->PAsphingo


PAsphingo[which(PAsphingo$genotype=="WT"),]->PAsphingo
#PAsphingo[which(PAsphingo$sphingo!="B"),]->PAsphingo


	
date2=format(Sys.Date(), format="%Y%m%d")


#copy PA for last day to PA_relative, and replace the pixels with % change of each plant 
#MAKES PA RELATIVE
PAsphingo[which(PAsphingo$date==20210712),]->PAsphingo_relative
for(p in plants_order){
	day0=PAsphingo$pixels[intersect(which(PAsphingo$date=="20210702"), which(PAsphingo$plantID==p))]
	day5=PAsphingo$pixels[intersect(which(PAsphingo$date=="20210712"), which(PAsphingo$plantID==p))]   #"20210712"
	(day5-day0)/day0->PAsphingo_relative$pixels[which(PAsphingo_relative$plantID==p)]
}
	


#RAW SIZE CHART. 
#	pdf(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/soil_protection_photos_202107/Sphingo_since_sphingo_spray_", date, ".pdf", sep=""), width = 5, height = 5, useDingbats=FALSE)
		#print(
		ggplot(PAsphingo_relative, aes(fill=pseudo, y=pixels, x=sphingo)) + 
   	 	geom_boxplot(outlier.size=0, outlier.shape=NA, alpha=0.75) + 
   	 	geom_point(aes(fill=pseudo), size = 1.2, stroke=0.1, shape = 21, position = position_jitterdodge(jitter.width=0.3)) +
   	 	theme_classic() + #scale_fill_manual(values=c("#ED6925FF","#460B5DFF","#ffff20","#FCA209FF","#67CC5CFF","#FB9606FF","#297B8EFF","#A22B62FF")) + 
    	#scale_y_continuous(name ="Leaf Area (square mm)", breaks = c(0, 100, 200, 300), labels = c(0, 100, 200, 300), limits = c(0,325)) +
    	#scale_x_discrete(name ="Pseudomonas Strain") +
    	theme(axis.text=element_text(size=15, color="black"), axis.title=element_text(size=14,face="bold")) + 
    	theme(axis.text.x = element_text(angle=90, hjust=0.99, vjust=0.9)) + 
    	theme(axis.ticks = element_line(colour = 'black'))
    	#)
#	dev.off()
		
MgCl=PAsphingo_relative[which(PAsphingo_relative$sphingo=="M" & PAsphingo_relative$genotype=="WT") ,"pixels"]
F=PAsphingo_relative[which(PAsphingo_relative$sphingo=="F" & PAsphingo_relative$genotype=="WT") ,"pixels"]
S=PAsphingo_relative[which(PAsphingo_relative$sphingo=="S" & PAsphingo_relative$genotype=="WT") ,"pixels"]
B=PAsphingo_relative[which(PAsphingo_relative$sphingo=="B" & PAsphingo_relative$genotype=="WT") ,"pixels"]

wilcox.test(MgCl, F, "greater")

	



gen="eds"
FR1=PA3[which(PA3$Sphingomonas=="FR1" & PA3$genotype==gen) ,"pixels"]
S18=PA3[which(PA3$Sphingomonas=="S18" & PA3$genotype==gen) ,"pixels"]
S131=PA3[which(PA3$Sphingomonas=="S131" & PA3$genotype==gen) ,"pixels"]
S63D7=PA3[which(PA3$Sphingomonas=="S63D7" & PA3$genotype==gen) ,"pixels"]
MgCl=PA3[which(PA3$Sphingomonas=="MgCl" & PA3$genotype==gen) ,"pixels"]
Sboil=PA3[which(PA3$Sphingomonas=="Sboil" & PA3$genotype==gen) ,"pixels"]

wilcox.test(FR1, S131, "greater")
wilcox.test(FR1, MgCl, "greater")



p.adjust(wilcox.test(FR1, MgCl, "greater")[[3]], method="fdr", n=5)

DC3000=PA2[which(PA2$Pseudomonas=="DC3000" & PA2$genotype=="eds") ,"pixels"]
OTU5=PA2[which(PA2$Pseudomonas=="OTU5" & PA2$genotype=="eds") ,"pixels"]
BOIL=PA2[which(PA2$Pseudomonas=="BOIL" & PA2$genotype=="eds") ,"pixels"]
wilcox.test(OTU5, DC3000, "greater")
wilcox.test(BOIL, DC3000, "greater")
wilcox.test(OTU5, BOIL, "greater")




S131=PA2[which(PA2$Sphingomonas=="S131" & PA2$genotype=="Ey") ,"pixels"]
S63D7=PA2[which(PA2$Sphingomonas=="S63D7" & PA2$genotype=="Ey") ,"pixels"]
MgCl=PA2[which(PA2$Sphingomonas=="MgCl" & PA2$genotype=="Ey") ,"pixels"]
Sboil=PA2[which(PA2$Sphingomonas=="Sboil" & PA2$genotype=="Ey") ,"pixels"]

#  5000*3300
#each original pixel = 0.12 mm
#each original square was 360 by 360 = 129600 pixels
#each original square was 43.2 by 43.2 = 1866.24 mm^2
#each square on Ilja's program-segmented image was 375 by 375 = 140625 pixels
#To convert from Ilja's image to the original square is 129600/140625 pixels
#To further convert from the original square pixels to mm^2 is (129600/140625) * (1866.24/129600)
#means final conversion is 0.01327104 to get pixels to mm^2 for soil

  
595*390=232050 mm^2

5000/3300
595/390

##### LOAD THE PA FROM THE EARLIER EXPERIMENT HERE, TO LINE 184
##########
##########
##########
##########
##########

date=format(Sys.Date(), format="%Y%m%d")


#REMOVE PLANTS WITH LESS THAN 500 PIXELS
PA[which(PA$day7>=500),]->PA  

plot(PA$day7)

options(scipen=999)

PA4=data.frame(Sphingomonas=PA$Sphingomonas, Pseudomonas=PA$Pseudomonas, Green_pixels=PA$day7)
colnames(PA4)=c("Sphingomonas", "Pseudomonas", "pixels")

256/90205->pixel_to_mm2

PA4$Pseudomonas=factor(PA4$Pseudomonas, levels=c("OTU5", "DC3000", "MgCl"))

#RAW SIZE CHART. 
	pdf(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/soil_protection_photos_20201007/Agar_OTU5_v_DC3000_", date, ".pdf", sep=""), width = 3.2, height = 5, useDingbats=FALSE)
		#print(
		ggplot(PA4, aes(fill=Pseudomonas, y=pixel_to_mm2*pixels, x=Pseudomonas)) + 
   	 	geom_boxplot(outlier.size=0, outlier.shape=NA, alpha=0.75) + 
   	 	geom_point(aes(fill=Pseudomonas), size = 1.2, stroke=0.1, shape = 21, position = position_jitterdodge(jitter.width=0.3)) +
   	 	theme_classic() + #scale_fill_manual(values=c("#ED6925FF","#460B5DFF","#ffff20","#FCA209FF","#67CC5CFF","#FB9606FF","#297B8EFF","#A22B62FF")) + 
    	scale_y_continuous(name ="Leaf Area (square mm)", breaks = c(0, 100, 200, 300), labels = c(0, 100, 200, 300), limits = c(0,350)) +
    	#scale_x_discrete(name ="Pseudomonas Strain") +
    	theme(axis.text=element_text(size=15, color="black"), axis.title=element_text(size=14,face="bold")) + 
    	theme(axis.text.x = element_text(angle=90, hjust=0.99, vjust=0.9)) + 
    	theme(axis.ticks = element_line(colour = 'black'))
    	#)
	dev.off()

