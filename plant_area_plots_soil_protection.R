library("ggplot2")
library("reshape2")

PA<-read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/soil_protection_photos_20201007/plant_area_pixels_20201021.txt") 
	as.matrix(PA[1,])->colnames(PA)
	as.matrix(PA[,1])->PA[,1]
	PA[2:nrow(PA),]->PA
	do.call(rbind, strsplit(PA[,1], split="_"))[,1]->PA$Sphingomonas
	do.call(rbind, strsplit(PA[,1], split="_"))[,2]->PA$genotype
	do.call(rbind, strsplit(PA[,1], split="_"))[,3]->PA$Pseudomonas
	as.numeric(as.matrix(PA$Green_pixels))->PA$Green_pixels
	
	
	#PA=rbind(c("FR1_Ey_BOIL_small", "4", 10000, "FR1", "Ey", "BOIL"), PA)
	
	factor(PA$Sphingomonas, levels=c("FR1", "MgCl", "S131", "S18", "S63D7", "Sboil"))->PA$Sphingomonas
	PA$Pseudomonas=gsub("boil", "BOIL", PA$Pseudomonas)
	factor(PA$Pseudomonas, levels=c("OTU5", "DC3000", "BOIL"))->PA$Pseudomonas
	factor(PA$genotype, levels=c("Ey", "eds"))->PA$genotype
	as.numeric(as.matrix(PA$Green_pixels))->PA$Green_pixels
	

date=format(Sys.Date(), format="%Y%m%d")



#REMOVE PLANTS WITH LESS THAN 500 PIXELS
PA[which(PA$Green_pixels>=500),]->PA  


options(scipen=999)


PA2=data.frame(Sphingomonas=PA$Sphingomonas, Pseudomonas=PA$Pseudomonas, genotype=PA$genotype, Green_pixels=PA$Green_pixels)
colnames(PA2)=c("Sphingomonas", "Pseudomonas", "genotype", "pixels")



#RAW SIZE CHART. 
	pdf(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/soil_protection_photos_20201007/Soil_OTU5_v_DC3000_", date, ".pdf", sep=""), width = 5, height = 5, useDingbats=FALSE)
		#print(
		ggplot(PA2, aes(fill=Pseudomonas, y=0.01327104*pixels, x=genotype)) + 
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



foc="Ey"
PA3=PA2[which(PA2$genotype==foc),]
PA3$Sphingomonas=factor(PA3$Sphingomonas, levels=c("FR1","S63D7", "S131",  "S18",   "Sboil", "MgCl" ))
	
PA3$Sphing_Pseud=paste(PA3[,1], PA3[,2], sep="_")

PA3$Sphing_Pseud=factor(PA3$Sphing_Pseud, levels=c( "FR1_BOIL",  "FR1_DC3000",  "FR1_OTU5",
									"S63D7_BOIL",  "S63D7_DC3000",  "S63D7_OTU5",
									"S131_BOIL",  "S131_DC3000",  "S131_OTU5",
									"S18_BOIL",  "S18_DC3000",  "S18_OTU5",
									"Sboil_BOIL",  "Sboil_DC3000",  "Sboil_OTU5",
									"MgCl_BOIL",  "MgCl_DC3000",  "MgCl_OTU5"))
									
PA3$Color=PA3$Sphing_Pseud

PA3$Pseudomonas=factor(PA3$Pseudomonas, levels=c("BOIL", "DC3000", "OTU5"))


#RAW SIZE CHART. 
	pdf(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/soil_protection_photos_20201007/Soil_Sphingo_", foc, "_", date, ".pdf", sep=""), width = 5, height = 5, useDingbats=FALSE)
		#print(
		ggplot(PA3, aes(fill=Pseudomonas, y=0.01327104*pixels, x=Sphingomonas)) + 
   	 	geom_boxplot(outlier.size=0, outlier.shape=NA, alpha=0.75) + 
   	 	geom_point(aes(fill=Pseudomonas), size = 1.2, stroke=0.1, shape = 21, position = position_jitterdodge(jitter.width=0.3)) +
   	 	theme_classic() + #scale_fill_manual(values=c("#ED6925FF","#460B5DFF","#ffff20","#FCA209FF","#67CC5CFF","#FB9606FF","#297B8EFF","#A22B62FF")) + 
    	scale_y_continuous(name ="Leaf Area (square mm)", breaks = c(0, 100, 200, 300), labels = c(0, 100, 200, 300), limits = c(0,325)) +
    	#scale_x_discrete(name ="Pseudomonas Strain") +
    	theme(axis.text=element_text(size=15, color="black"), axis.title=element_text(size=14,face="bold")) + 
    	theme(axis.text.x = element_text(angle=90, hjust=0.99, vjust=0.9)) + 
    	theme(axis.ticks = element_line(colour = 'black'))
    	#)
	dev.off()



  #Ey
PA5=PA2

gen="Ey" 
Pseudomonas="DC3000"


PA5=PA5[which(PA5$genotype==gen),]
PA5=PA5[which(PA5$Pseudomonas==Pseudomonas),]


FR1=PA5[which(PA5$Sphingomonas=="FR1") ,"pixels"]
S18=PA5[which(PA5$Sphingomonas=="S18") ,"pixels"]
S131=PA5[which(PA5$Sphingomonas=="S131") ,"pixels"]
S63D7=PA5[which(PA5$Sphingomonas=="S63D7") ,"pixels"]
MgCl=PA5[which(PA5$Sphingomonas=="MgCl") ,"pixels"]
Sboil=PA5[which(PA5$Sphingomonas=="Sboil") ,"pixels"]

wilcox.test(S18, MgCl, "greater")[[3]]
wilcox.test(S18, Sboil, "greater")[[3]]
wilcox.test(S131, MgCl, "greater")[[3]]
wilcox.test(S131, Sboil, "greater")[[3]]
wilcox.test(S63D7, MgCl, "greater")[[3]]
wilcox.test(S63D7, Sboil, "greater")[[3]]
wilcox.test(FR1, Sboil, "greater")[[3]]
wilcox.test(FR1, MgCl, "greater")[[3]]

#DC3000 
p.adjust(c(0.02308, 0.002322, 0.0003898, 0.02541), method="fdr")

#OTU5
p.adjust(c(0.00896, 0.01104, 0.00896, 0.01104), method="fdr")


p.adjust(c(0.02308, 0.002322, 0.0003898, 0.02541, 0.00896, 0.01104, 0.00896, 0.01104), method="fdr")



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

