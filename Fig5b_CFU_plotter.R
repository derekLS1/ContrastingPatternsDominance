source("/Users/dlundberg/Documents/abt6/scripts/R/functions/microbiome_custom_functions.R")

CFU<-read.table("/Users/dlundberg/Documents/abt6/Pratchaya/Bulk_metagenomes/Final_CFU_counts_SpringSummer.txt", sep="\t")
as.matrix(CFU[1,])->names(CFU)
CFU[2:nrow(CFU),]->CFU

as.numeric(as.matrix(CFU$Pseudo1))->CFU$Pseudo1
as.numeric(as.matrix(CFU$Pseudo2))->CFU$Pseudo2
as.numeric(as.matrix(CFU$LB))->CFU$LB
as.numeric(as.matrix(CFU$Sphingo))->CFU$Sphingo

rep("#000000", nrow(CFU))->CFU$PlantColor

for(clr in 1:nrow(plant_color_list)){
	print(clr)
	number_matches=length(grep(plant_color_list[clr,1], CFU$Plant))
	rep(plant_color_list[clr, 2], number_matches)-> CFU$PlantColor[grep(plant_color_list[clr,1], CFU$Plant)]
}


library(ggplot2)
library(tidyr)
library(reshape2)
melt(CFU)->CFU2
c("Plant","PlantColor", "Media","cfu")->names(CFU2)

CFU2$cfu+1->CFU2$cfu


gsub("_S.*", "_S", CFU2$Plant)->CFU2$Plant
gsub("_Og.*", "_Og", CFU2$Plant)->CFU2$Plant

CFU2[which(CFU2$Media!="Pseudo1"),]->CFU2


reorder_rows=c(
	grep("Athaliana", CFU2$Plant),
	grep("Draba", CFU2$Plant),
	grep("Cardamine", CFU2$Plant),
	grep("Dandelion", CFU2$Plant),
	grep("Thistle", CFU2$Plant),	
	grep("Grass", CFU2$Plant),
	grep("NeckarKing", CFU2$Plant),
	grep("NK", CFU2$Plant),
	grep("Clover", CFU2$Plant),
	grep("Fuzzy", CFU2$Plant),
	grep("Plantago", CFU2$Plant),
	grep("SpringKraut", CFU2$Plant),
	grep("Moss", CFU2$Plant),
	grep("Soil", CFU2$Plant), 
	grep("ShortTree", CFU2$Plant)
)	
CFU2[reorder_rows,]->CFU2
reorder_rows=c(
	grep("_Og", CFU2$Plant),
	grep("_S", CFU2$Plant)
)
CFU2[reorder_rows,]->CFU2

#REMOVE SHORT TREE
CFU2[grep("ShortTree", CFU2$Plant, invert=TRUE),]->CFU2
CFU2[grep("BigClover", CFU2$Plant, invert=TRUE),]->CFU2

factor(CFU2$Plant, levels=unique(as.character(as.matrix(CFU2$Plant))))->CFU2$Plant



factor(CFU2$Media, levels=c("LB", "Sphingo", "Pseudo2"))->CFU2$Media


#adjust CFU for grams of plant tissue, given that only 50uL of glycerol stock were plated
#CFU2$cfu*220.5882353->CFU2$cfu


CFU2[grep("_Og", CFU2$Plant),]->CFU2spring
CFU2[grep("_S", CFU2$Plant),]->CFU2summer



pdf(file="/Users/dlundberg/Documents/abt6/Pratchaya/Bulk_metagenomes/spring_cfu.pdf", width = 8, height = 2.5, useDingbats=FALSE)
ggplot(CFU2spring, aes(x = Plant, y = log10(cfu), fill=Media)) + 
                   geom_boxplot(outlier.size=0) + 
                   geom_point(aes(fill=Media, alpha=0.5), size = 2, shape = 21, position = position_jitterdodge()) + 
                  theme_classic() + scale_fill_manual(values=c("#E3E3E3", "#ED2124", "#3382BE")) +
                  scale_x_discrete(name = "Plant Species") + 
                  scale_y_continuous(limits=c(0, 7), 
                  	breaks = c(0, 1, 2, 3, 4, 5, 6, 7), 
                  	labels=10^c(0, 1, 2, 3, 4, 5, 6, 7)) +
                  guides(fill=guide_legend(title="Treatments")) + 
                	theme(axis.text.x = element_text(angle = 0, vjust=1, hjust = 1)) + 
                	theme(axis.text=element_text(size=15, color="black"), axis.title=element_text(size=14,face="bold"))          
dev.off()

pdf(file="/Users/dlundberg/Documents/abt6/Pratchaya/Bulk_metagenomes/summer_cfu.pdf", width = 8, height = 2.5, useDingbats=FALSE)
ggplot(CFU2summer, aes(x = Plant, y = log10(cfu), fill=Media)) + 
                   geom_boxplot(outlier.size=0) + 
                   geom_point(aes(fill=Media, alpha=0.5), size = 2, shape = 21, position = position_jitterdodge()) + 
                  theme_classic() + scale_fill_manual(values=c("#E3E3E3", "#ED2124", "#3382BE")) +
                  scale_x_discrete(name = "Plant Species") + 
                  scale_y_continuous(limits=c(0, 7), 
                  	breaks = c(0, 1, 2, 3, 4, 5, 6, 7), 
                  	labels=10^c(0, 1, 2, 3, 4, 5, 6, 7)) +
                  guides(fill=guide_legend(title="Treatments")) + 
                	theme(axis.text.x = element_text(angle = 0, vjust=1, hjust = 1)) + 
                	theme(axis.text=element_text(size=15, color="black"), axis.title=element_text(size=14,face="bold"))          
dev.off()       



   labels = scales::trans_format("log10", scales::math_format(10^.x)),
                                  name="log10(cfu+1)") +  





new <- total_cfu[71:147,]

library(FSA)
library(rcompanion)
kruskal.test(cfu ~ Plant, data = new)
DT=dunnTest(cfu ~ Plant,
            data=new,
            method="bh")
DT
PT = DT$res

PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.01)







#Anova no pq no es parametric....
#Statistical analysis
#ANOVA (if both conditions are met)
res.aov <- aov(cfu ~ Plant*Treatment, data = total_cfu)
summary(res.aov)
TukeyHSD(res.aov)

#Another way to do the Tuckey
library(multcomp)
summary(glht(res.aov, linfct = mcp(Plant = "Tukey")))


 #Normal distribution?
   # Extract the residuals
   aov_residuals <- residuals(object = res.aov )
   # Run Shapiro-Wilk test
   shapiro.test(x = aov_residuals )
   plot(res.aov, 2)
 
 #Homogeneity test
   library(car)
   leveneTest(cfu ~ Plant*Treatment, data = total_cfu)
 
 
