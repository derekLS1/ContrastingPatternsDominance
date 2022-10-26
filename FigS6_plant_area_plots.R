install.packages("ape")
library(ape)
#expDate=20190506 
expDate=20190611


if(expDate==20190506){
	PA<-read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/random_design_tables/plant_area_20190506.txt") 
	as.matrix(PA[1,])->colnames(PA)
	PA[2:nrow(PA),]->PA
	factor(PA$Plate)->PA$Plate
	factor(PA$Pot)->PA$Pot
	factor(PA$PlatePot)->PA$PlatePot
	factor(PA$Sphingomonas, levels=c(1:21))->PA$Sphingomonas
	factor(PA$Pseudomonas, levels=c("MgCl", "DC3000","OTU5"))->PA$Pseudomonas
	day0=PA$d20190429=as.numeric(as.matrix(PA$d20190429))
	day2=PA$d20190501=as.numeric(as.matrix(PA$d20190501))
	day4=PA$d20190503=as.numeric(as.matrix(PA$d20190503))
	day7=PA$d20190506=as.numeric(as.matrix(PA$d20190506))
	c("day0","day2","day4","day7")->names(PA)[c(6,7,8,9)]
}
if(expDate==20190611){
	PA<-read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/random_design_tables/plant_area_20190611.txt") 
	as.matrix(PA[1,])->colnames(PA)
	PA[2:nrow(PA),]->PA
	factor(PA$Plate)->PA$Plate
	factor(PA$Pot)->PA$Pot
	factor(PA$PlatePot)->PA$PlatePot
	factor(PA$Sphingomonas, levels=c(1:6))->PA$Sphingomonas
	PA$Pseudomonas=gsub("p25.c2", "OTU5", PA$Pseudomonas)
	factor(PA$Pseudomonas, levels=c("MgCl", "DC3000","OTU5"))->PA$Pseudomonas
	day0=PA$d20190603=as.numeric(as.matrix(PA$d20190603))
	day2=PA$d20190605=as.numeric(as.matrix(PA$d20190605))
	day4=PA$d20190607=as.numeric(as.matrix(PA$d20190607))
	day7=PA$d20190611=as.numeric(as.matrix(PA$d20190611))
	c("day0","day2","day4","day7")->names(PA)[c(6,7,8,9)]
}
	
date=format(Sys.Date(), format="%Y%m%d")


#load core genome tree so you can order stuff
	basefolder="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS/"
	core=read.tree(file = paste(basefolder, "strain_tree.nwk", sep="", collapse=""))


#REMOVE PLANTS WITH LESS THAN 4000 PIXELS
PA[which(day0>=4000),]->PA  


sort(table(paste(PA[,4], PA[,5], sep="_")))

options(scipen=999)

if(expDate==20190506){ 
	samples=1:21
	Snames=c("S230H133","S213H113","S133H113","S111H133","S194H133","S139H133","S18H113","S237H113","S127H113","S132H113","S190H113","S216H113","S230H113","FR1","S246H133","S380H113","S136H113","S101H133","S98H133","S337H113","BLANK")
}
if(expDate==20190611){ 
	samples=1:6
	Snames=c("S18H113","S246H133", "FR1", "S139H133", "S127H113", "BLANK")
}

testvalues3=testvalues2=testvalues1=samples
for(s in samples){
	print(s)
	Subset=PA[which(PA$Sphingomonas==s),]
		Subset[is.na(Subset$day7)==FALSE,]->Subset
		SubsetOTU5=Subset[which(Subset$Pseudomonas=="OTU5"),]
		SubsetMgCl=Subset[c(which(Subset$Pseudomonas=="DC3000"), which(Subset$Pseudomonas=="MgCl")),]
	CtrlSubset=PA[which(PA$Sphingomonas==length(samples)),]
		CtrlSubset[is.na(CtrlSubset$day7)==FALSE,]->CtrlSubset
		CtrlSubsetOTU5=CtrlSubset[which(CtrlSubset$Pseudomonas=="OTU5"),]
		CtrlSubsetMgCl=CtrlSubset[which(CtrlSubset$Pseudomonas=="MgCl"),]
		
		subject=abs(SubsetOTU5$day7-median(SubsetMgCl$day7))
		query=abs(CtrlSubsetOTU5$day7-median(CtrlSubsetMgCl$day7))
	
		testvalues1[s]=wilcox.test(subject, query, alternative="less", paired=FALSE)[[3]]
		testvalues2[s]=wilcox.test(SubsetOTU5$day0, CtrlSubsetOTU5$day0, alternative="less", paired=FALSE)[[3]]
		testvalues3[s]=wilcox.test(SubsetOTU5$day7, SubsetMgCl$day7, alternative="less", paired=FALSE)[[3]]   #2021 redo
	
	names(testvalues1)[s]=names(testvalues2)[s]=names(testvalues3)[s]=Snames[s]
}
adj_testvalues1=p.adjust(testvalues1,method="fdr")
testvalues1
adj_testvalues1[which(adj_testvalues1<=0.05)]
adj_testvalues2=p.adjust(testvalues2,method="fdr")
testvalues2
adj_testvalues2[which(adj_testvalues2<=0.05)]
adj_testvalues3=p.adjust(testvalues3, method="fdr")
adj_testvalues3
adj_testvalues3[which(adj_testvalues3>=0.05)]


#uniquenames=c("S216H113", "S127H113", "FR1", "S139H133", "S246H133", "BLANK")



PA$PlatePot[intersect(which(PA$Sphingomonas==16), which(PA$Pseudomonas=="OTU5"))]
PA$PlatePot[intersect(which(PA$Sphingomonas==16), which(PA$Pseudomonas=="MgCl"))]


#CALCULATE THE RATIO CHANGE AS OPPOSED TO ONLY RAW PIXELS
PA$difference_day0= 100*(PA$day0 - PA$day0)/ PA$day0 
PA$difference_day2= 100*(PA$day2 - PA$day0)/ PA$day0 
PA$difference_day4= 100*(PA$day4 - PA$day0)/ PA$day0 
PA$difference_day7= 100*(PA$day7 - PA$day0)/ PA$day0 

if(expDate==20190506){ 
	samples=1:21
	Snames=c("S230H133","S213H113","S133H113","S111H133","S194H133","S139H133","S18H113","S237H113","S127H113","S132H113","S190H113","S216H113","S230H113","FR1","S246H133","S380H113","S136H113","S101H133","S98H133","S337H113","BLANK")
}
if(expDate==20190611){ 
	samples=1:6
	Snames=c("S18H113","S246H133", "FR1", "S139H133", "S127H113", "BLANK")
}
testvalues=samples
for(s in samples){
	print(s)
	Subset=PA[which(PA$Sphingomonas==s),]
		Subset[is.na(Subset$difference_day7)==FALSE,]->Subset
		SubsetOTU5=Subset[which(Subset$Pseudomonas=="OTU5"),]
		SubsetMgCl=Subset[which(Subset$Pseudomonas=="MgCl"),]
	CtrlSubset=PA[which(PA$Sphingomonas==length(samples)),]
		CtrlSubset[is.na(CtrlSubset$difference_day7)==FALSE,]->CtrlSubset
		CtrlSubsetOTU5=CtrlSubset[which(CtrlSubset$Pseudomonas=="OTU5"),]
		CtrlSubsetMgCl=CtrlSubset[which(CtrlSubset$Pseudomonas=="MgCl"),]
		
		subject=abs(SubsetOTU5$difference_day7-median(SubsetMgCl$difference_day7))
		query=abs(CtrlSubsetOTU5$difference_day7-median(CtrlSubsetMgCl$difference_day7))
		#testvalues[s]=wilcox.test(SubsetOTU5$difference_day7, SubsetOTU5$difference_day0, alternative="greater", paired=TRUE)[[3]]
	
		testvalues[s]=wilcox.test(subject, query, alternative="less", paired=FALSE)[[3]]
	
	names(testvalues)[s]=Snames[s]
}
adj_testvalues=p.adjust(testvalues,method="fdr")
testvalues
adj_testvalues
adj_testvalues[which(adj_testvalues<=0.05)]


PA[PA$Pseudomonas=="MgCl",]->PAanova
	Snames[as.numeric(as.matrix(PAanova$Sphingomonas))]->PAanova$Sphingomonas
	factor(PAanova$Sphingomonas, levels=Snames)->PAanova$Sphingomonas
	#growth.promotion = lm(Sphingomonas, data = PAanova)
 	#summary(growth.promotion)
	P=aov(day7~Sphingomonas, data = PAanova)
	tuki=TukeyHSD(P)[[1]]
	tuki[order(tuki[,4], decreasing=TRUE),]->tuki
	tuki[grep("BLANK", rownames(tuki)),]->tuki
	tuki[which(tuki[,4]<0.05),]
PA[PA$Pseudomonas=="OTU5",]->PAanova
	Snames[as.numeric(as.matrix(PAanova$Sphingomonas))]->PAanova$Sphingomonas
	factor(PAanova$Sphingomonas, levels=Snames)->PAanova$Sphingomonas
	#growth.promotion = lm(Sphingomonas, data = PAanova)
 	#summary(growth.promotion)
	P=aov(difference_day7~Sphingomonas, data = PAanova)
	tuki=TukeyHSD(P)[[1]]
	tuki[order(tuki[,4], decreasing=TRUE),]->tuki
	tuki[grep("BLANK", rownames(tuki)),]->tuki
	tuki[which(tuki[,4]<0.05),]




library("reshape2")
PA2=data.frame(Sphingomonas=PA$Sphingomonas, Pseudomonas=PA$Pseudomonas, difference_day2=PA$difference_day2, difference_day4=PA$difference_day4, difference_day7=PA$difference_day7)
PA2=melt(PA2)
colnames(PA2)=c("Sphingomonas", "Pseudomonas", "date", "pixels")
#paste(PA2$Pseudomonas, PA2$date, sep="")->PA2$PseudoDate

library("ggplot2")
#colors=c("#000000", "#082FE0","#B81155", "#F66817") #Sphingomonas
#colors=c("#757575", "#0AB2D0","#68d600", "#439255") #Pseudomonas

#PA[PA$Pseudomonas=="MgCl",]->PA_MgCl
#PA[PA$Pseudomonas=="OTU5",]->PA_OTU5
#PA[PA$Pseudomonas=="DC3000",]->PA_DC3000

#convert numbers 1-21 into the actual names of the Sphingomonas strains
if(expDate==20190506){ 
	uniquenames=c("S230H133","S213H113","S133H113","S111H133","S194H133","S139H133","S18H113","S237H113","S127H113","S132H113","S190H113","S216H113","S230H113","FR1","S246H133","S380H113","S136H113","S101H133","S98H133","S337H113","BLANK")
	uniquenames=cbind(1:21, uniquenames)
}
if(expDate==20190611){ 
	uniquenames=c("S18H113","S246H133", "FR1", "S139H133", "S127H113", "BLANK")
	uniquenames=cbind(1:6, uniquenames)
}






#RECOLOR SAMPLE NAMES BASED ON ASV
culture_fasta=read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/KarasovAlmario16S/all_cultured_16S_20190410_reoriented_trimmed.fa")
as.matrix(culture_fasta)->culture_fasta
	culture_fasta[sort(c(grep("Sphingo", culture_fasta), (grep("Sphingo", culture_fasta)+1))),]->Sphingo16S_culture
		#remove contaminants
		contaminants=c("S173H113","S20H113","S375H113","S26H113","S104H113","S106H113","S107H113","S112H113","S114H113","S115H113","S260H113","S303H113","S313H113","S1H113","S305H113","S149H113","S379H113","S113H113","S159H113","S162H113")
		kill=vector(length=0)
		for(c in contaminants){
			c(kill, grep(c, Sphingo16S_culture))->kill
		}
		kill=sort(c(kill, kill+1))
		Sphingo16S_culture[setdiff(1:length(Sphingo16S_culture), kill)]->Sphingo16S_culture
		split_FASTA(Sphingo16S_culture, add_integer_column=FALSE, leave_carrot=FALSE, reverse=FALSE)->Sphingo16S_table
		gsub("H133.Sphingo_16S", "H133" , Sphingo16S_table[,1])->Sphingo16S_table[,1]
		gsub("H113.Sphingo_16S", "H113" , Sphingo16S_table[,1])->Sphingo16S_table[,1]
		#Turn tip labels into the actual V3V4 ASV sequences
		Sphingo16S_table[match(uniquenames[,2], Sphingo16S_table[,1]),2]->associated_seqs_of_labels
		sphingo_16Sseqs_to_color=as.matrix(read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/sphingo_16Sseqs_to_color.txt", comment.char=""))
		associated_colors_of_labels=sphingo_16Sseqs_to_color[match(associated_seqs_of_labels, sphingo_16Sseqs_to_color[,1]),4]
		associated_colors_of_labels[is.na(associated_colors_of_labels)]<-"#000000"

#ADD SEQ AND COLOR TO THE TABLE OF SPHINGOMONAS STRAINS PHENOTYPED
plotting_key=cbind(uniquenames, associated_seqs_of_labels, associated_colors_of_labels)



#order samples based on a txt file representing the order of the core genome tree
#treeorder<-as.matrix(read.table("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/PanX/Sphingo_runs/Sph_busco85_wONT_lessREFSEQ_noIBVSS/strain_tree_order.txt", sep="\t", quote=""))
#uniquenames[match(gsub(".*-S", "S", treeorder),uniquenames[,2], nomatch=0),]->uniquenames
#setdiff(uniquenames[,2], gsub(".*-S", "S", treeorder))

if(expDate==20190506){ 
	goodorder=c("S230H133","S216H113","S133H113","S127H113","S18H113", "S194H133","S213H113","S190H113","S230H113","S237H113","S380H113","S136H113","S132H113", "S337H113","S111H133","S246H133","S139H133" ,"S98H133", "S101H133", "FR1",  "BLANK")
	breaks=c(1:21)
	wid=10
}
if(expDate==20190611){ 
	goodorder=c("S127H113", "S18H113", "S246H133", "S139H133", "FR1", "BLANK")
	breaks=c(1:6)
	wid=5
}

#match(goodorder, core$tip.label)

uniquenames[match(goodorder, uniquenames[,2]),1]

xaxisorder=uniquenames[match(goodorder, uniquenames[,2]),1]
factor(PA$Sphingomonas, levels=xaxisorder)->PA$Sphingomonas


#on day7, well diameter = 305 pixels or 16mm.
#means 305^2 = 93025 pixels = 16^2 = 256mm2.  So 256/90205 is conversion

256/90205->pixel_to_mm2

100000*pixel_to_mm2

#RAW SIZE CHART. 
for(DATE in c("day0", "day2", "day4", "day7")){
	#pdf(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/random_design_tables/raw_size_", expDate, "_", DATE, ".pdf", sep=""), width = wid, height = 5, useDingbats=FALSE)
	#	print(
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
}


#PERCENT SIZE CHART
for(date in c("difference_day0", "difference_day2", "difference_day4", "difference_day7")){
	pdf(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/random_design_tables/percent_size_", expDate, "_", date, ".pdf", sep=""), width = wid, height = 5, useDingbats=FALSE)
		print(
		ggplot(PA, aes(fill=Pseudomonas, y=get(date), x=Sphingomonas)) + 
   	 	geom_boxplot(outlier.size=0, outlier.shape=NA, alpha=0.75) + 
   	 	geom_point(aes(fill=Pseudomonas), size = 1.2, stroke=0.3, shape = 21, position = position_jitterdodge(jitter.width=0.3)) +
    	theme_classic() + scale_fill_manual(values=c("#C576B0", "#6DE9FF", "#F7EC2A")) + 
    	scale_y_continuous(name ="Percent Change", breaks = c(-100, 0, 100, 200, 300, 400), labels = c(-100, 0, 100, 200, 300, 400), limits = c(-100,400)) +
    	scale_x_discrete(name ="Sphingomonas Strain", breaks = breaks, labels = plotting_key[,2], ) +
    	theme(axis.text=element_text(size=15, color="black"), axis.title=element_text(size=14,face="bold")) + 
    	theme(axis.text.x = element_text(angle=90, hjust=0.99, vjust=0.9, color=plotting_key[,4])) + 
    	geom_hline(yintercept=0, linetype="dashed", color="red", size=1) + 
    	theme(axis.ticks = element_line(colour = 'black'))
		)
	dev.off()
}





#### PHOTO MONTAGES::::

#scp dlundberg@burrito.eb.local:/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/Phenotyping_Sphingo/Experiment1_21_strains/photos_20190506_new/processed/p15_P1040696/seg/p15_P1040696_002_seg.png /Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/random_design_tables/day6_photo_montage 

    

#6,   14,   16,   21
strain=14
date=20190506_new
Pseudo=DC3000
if [ $strain == 2 ]
then
potsOTU5=$(echo pl1po11 pl3po13 pl4po9 pl5po8 pl5po14 pl6po19 pl7po8 pl7po12 pl8po18 pl9po2 pl9po11 pl11po8 pl11po19 pl12po22 pl14po8 pl14po14 pl15po5 pl15po6 pl16po21 pl19po0 pl19po15 pl20po18 pl21po13 pl21po23)
potsMgCl=$(echo pl22po1 pl23po9 pl24po6 pl24po22 pl25po7 pl25po18 pl26po21 pl27po9 pl29po1 pl29po18 pl32po22 pl33po16 pl33po21 pl34po3 pl35po10 pl35po15)
fi
if [ $strain == 3 ]
then
potsOTU5=$(echo pl1po22   pl2po11   pl3po0    pl3po3    pl4po16   pl4po18   pl4po20   pl7po22   pl8po19   pl9po17   pl10po3   pl11po16  pl11po17  pl14po2   pl14po15  pl14po23 pl16po3   pl17po3   pl17po12  pl18po17  pl19po11  pl19po12  pl20po20  pl20po23  pl312po21)
potsMgCl=$(echo pl25po13 pl26po8  pl26po11 pl26po13 pl26po17 pl27po7  pl28po6  pl28po9  pl29po11 pl31po21 pl32po11 pl34po7  pl34po16 pl34po23 pl35po6  pl35po8)
potsDC3000=$(echo pl36po2 pl37po7 pl37po16 pl38po7 pl39po7 pl39po23 pl41po1 pl42po21)
fi
if [ $strain == 4 ]
then
potsOTU5=$(echo pl1po4   pl1po16  pl3po5   pl3po18  pl4po15  pl6po4   pl6po18  pl7po13  pl10po7  pl10po22 pl11po15 pl12po8  pl12po15 pl13po2  pl14po12 pl14po18 pl15po17 pl16po10 pl16po12 pl18po23 pl21po2  pl21po3  pl21po14 pl21po21)
potsMgCl=$(echo pl23po7  pl23po8  pl24po10 pl25po19 pl27po11 pl27po23 pl28po2  pl28po11 pl29po3  pl29po4  pl29po23 pl30po11 pl32po16 pl33po19 pl35po19 pl35po21)
fi
if [ $strain == 5 ]
potsOTU5=$(echo pl1po1   pl2po0   pl3po10  pl3po21  pl4po7   pl4po19  pl5po7   pl7po4   pl8po1   pl8po3   pl9po15  pl11po21 pl12po14 pl12po23 pl13po17 pl13po22 pl15po7  pl16po4 pl16po11 pl16po13 pl19po1  pl19po4  pl19po7  pl19po23 pl312po7)
potsMgCl=$(echo pl22po2  pl23po22 pl24po2  pl24po8  pl25po1  pl26po3  pl27po2  pl28po15 pl30po10 pl30po12 pl30po19 pl31po7  pl33po0  pl34po0  pl34po17 pl34po20)
fi
if [ $strain == 6 ]   #means S63D7
then
potsOTU5=$(echo pl1po23   pl2po13   pl4po04    pl5po00    pl5po05    pl5po10   pl7po05    pl7po23   pl11po13  pl11po20  pl12po11  pl13po03   pl13po15  pl14po03  pl14po09   pl14po16  pl16po18  pl16po20  pl18po19  pl19po05   pl20po01   pl20po04   pl21po01   pl21po06   pl312po14) #S6
potsMgCl=$(echo pl22po17 pl23po02  pl23po04  pl23po19 pl24po11 pl25po22 pl26po16 pl27po20 pl29po17 pl30po02  pl31po14 pl32po08  pl32po15 pl33po18 pl34po14)
potsDC3000=$(echo pl36po9 pl37po3 pl37po17 pl38po3 pl39po19 pl42po3 pl42po15 pl42po20)
fi
if [ $strain == 7 ]
then
potsOTU5=$(echo  pl2po7    pl3po23   pl4po6    pl4po23   pl5po20   pl6po7    pl6po20   pl7po21   pl8po8    pl8po9    pl12po5   pl12po18  pl13po1   pl13po7   pl13po8   pl13po9  pl15po2   pl15po4   pl16po17  pl18po8   pl18po9   pl19po21  pl20po22  pl21po4   pl312po13 pl312po23  )
potsMgCl=$(echo pl22po18 pl24po15 pl25po15 pl25po23 pl27po4  pl28po4  pl28po13 pl28po23 pl29po21 pl30po1  pl30po16 pl31po13 pl31po23 pl32po3  pl33po23 pl35po23  )
fi
if [ $strain == 8 ]
then
potsOTU5=$(echo pl1po20   pl2po19   pl3po15   pl6po10   pl6po14   pl6po16   pl6po17   pl7po11   pl8po13   pl8po21   pl9po14   pl10po23  pl11po22  pl13po0   pl13po6   pl14po21 pl15po0   pl15po1   pl15po16  pl15po19  pl16po7   pl17po13  pl19po9   pl20po12  pl312po18)
potsMgCl=$(echo pl22po16 pl23po13 pl25po8  pl26po10 pl26po19 pl27po18 pl28po17 pl28po18 pl29po7  pl30po14 pl30po18 pl30po21 pl31po18 pl33po10 pl35po9  pl35po14   )
fi
if [ $strain == 9 ]
then
potsOTU5=$(echo  pl2po3   pl3po1   pl5po13  pl6po2   pl6po13  pl7po1   pl7po3   pl7po14  pl9po0   pl10po21 pl11po2  pl12po10 pl12po16 pl13po10 pl14po11 pl14po17 pl15po9  pl16po8  pl16po19 pl17po6  pl17po17 pl17po18 pl18po21 pl21po20 )
potsMgCl=$(echo pl22po22 pl23po3  pl23po18 pl23po23 pl24po4  pl24po14 pl24po20 pl26po18 pl30po17 pl32po5  pl32po14 pl32po21 pl33po20 pl34po1  pl34po19 pl35po18  )
fi
if [ $strain == 10 ]
then
potsOTU5=$(echo pl2po14   pl3po4    pl3po16   pl4po11   pl4po12   pl5po3    pl6po21   pl7po17   pl9po9    pl10po6   pl10po9   pl10po16  pl10po18  pl12po19  pl13po20  pl13po21  pl15po14  pl16po2   pl16po14  pl17po8   pl17po20  pl18po5   pl18po16  pl19po22  pl312po2  pl312po15  )
potsMgCl=$(echo pl22po15 pl25po2  pl25po11 pl25po14 pl25po20 pl26po2  pl27po8  pl27po21 pl30po0  pl31po2  pl31po15 pl33po1  pl33po7  pl34po5  pl34po9  pl35po12  )
fi
if [ $strain == 11 ]
then
potsOTU5=$(echo pl1po8   pl2po23  pl3po19  pl3po20  pl4po21  pl6po9   pl6po23  pl7po16  pl9po5   pl9po8   pl11po5  pl11po9  pl12po4  pl13po12 pl14po5  pl14po20 pl15po10 pl16po1  pl16po5  pl18po7  pl19po3  pl20po6  pl21po5  pl21po15 pl312po9  )
potsMgCl=$(echo pl24po12 pl25po5  pl26po4  pl26po6  pl26po9  pl26po14 pl26po20 pl28po22 pl31po9  pl32po0  pl32po9  pl32po12 pl32po20 pl33po2  pl33po5  pl33po14  )
fi
if [ $strain == 12 ]
then
potsOTU5=$(echo  pl1po0    pl2po22   pl3po14   pl3po17   pl4po0    pl6po1    pl7po0    pl7po18   pl8po6    pl8po14   pl9po16   pl9po19   pl10po17  pl11po1   pl11po6   pl13po23  pl15po13  pl15po20  pl18po22  pl19po8   pl19po14  pl19po19  pl20po11  pl20po14  pl312po5  pl312po16 )
potsMgCl=$(echo  pl22po6  pl23po11 pl24po0  pl24po17 pl26po1  pl27po3  pl27po13 pl27po15 pl30po15 pl30po22 pl31po5  pl31po16 pl32po4  pl33po3  pl34po6  pl35po7  )
fi
if [ $strain == 13 ]
then
potsOTU5=$(echo  pl1po3   pl2po16  pl2po17  pl3po6   pl5po18  pl5po19  pl5po21  pl6po8   pl8po10  pl9po1   pl9po3   pl9po6   pl9po13  pl10po20 pl11po10 pl12po7  pl15po12 pl17po19 pl19po16 pl20po5  pl20po17 pl21po9  pl21po12 pl21po19 pl312po4)
potsMgCl=$(echo pl24po18 pl25po0  pl25po10 pl27po1  pl27po5  pl28po1  pl28po10 pl28po21 pl30po7  pl31po4  pl32po6  pl32po10 pl33po12 pl34po15 pl35po4  pl35po17  )
fi
if [ $strain == 14 ]   #means FR1
then
potsOTU5=$(echo pl1po02   pl1po15  pl2po21  pl8po23  pl9po04   pl10po00  pl13po04  pl14po19 pl17po05  pl17po10 pl18po13 pl19po13 pl19po18)  #S21
potsMgCl=$(echo pl22po03  pl22po11 pl22po13 pl27po14 pl30po04  pl31po17 pl35po11 pl35po22)  #S21
potsDC3000=$(echo pl37po12 pl38po9 pl38po17 pl39po1 pl40po1 pl41po0 pl41po20 pl42po8)
fi
if [ $strain == 15 ]
then
potsOTU5=$(echo pl1po6    pl1po21   pl2po2    pl2po8    pl5po9    pl5po16   pl6po12   pl7po19   pl8po15   pl8po17   pl8po20   pl10po4   pl10po11  pl10po13  pl11po14  pl14po7  pl15po3   pl16po15  pl17po1   pl17po7   pl17po23  pl18po12  pl19po10  pl20po3   pl312po1  pl312po3  pl312po8  pl312po11)
potsMgCl=$(echo pl22po7  pl22po9  pl22po21 pl23po1  pl23po5  pl25po16 pl28po3  pl28po20 pl29po0  pl31po1  pl31po3  pl31po8  pl31po11 pl32po18 pl35po5  pl35po13)
fi
if [ $strain == 16 ]
then
potsOTU5=$(echo pl1po10  pl3po08   pl4po10  pl4po22  pl6po03   pl6po06   pl7po02   pl7po09   pl8po05   pl9po07   pl9po18  pl9po21  pl12po13 pl12po21 pl13po18 pl14po04  pl16po09  pl17po04  pl17po09  pl18po04  pl18po10 pl18po11 pl19po17)  #S16
potsMgCl=$(echo pl22po8  pl23po14 pl23po16 pl24po7  pl26po0  pl26po12 pl27po19 pl28po12 pl29po19 pl30po8  pl30po20 pl30po23 pl32po7  pl32po23 pl33po6  pl33po17)  #S16
fi
if [ $strain == 17 ]
then
potsOTU5=$(echo pl1po7    pl1po9    pl2po6    pl2po10   pl2po18   pl4po8    pl5po11   pl6po22   pl9po10   pl10po19  pl11po7   pl12po1   pl12po6   pl13po19  pl14po1   pl14po10  pl15po22  pl16po22  pl17po2   pl18po3   pl18po18  pl20po15  pl21po17  pl21po18  pl312po19)
potsMgCl=$(echo pl22po5  pl25po4  pl25po17 pl25po21 pl27po6  pl27po17 pl28po8  pl28po14 pl29po16 pl29po22 pl31po19 pl33po4  pl33po9  pl33po13 pl33po22 pl34po2)
fi
if [ $strain == 18 ]
then
potsOTU5=$(echo pl1po13  pl1po18  pl2po9   pl3po12  pl4po3   pl5po15  pl5po17  pl7po10  pl8po4   pl9po20  pl9po23  pl10po8  pl11po4  pl12po0  pl12po9  pl12po20 pl16po0  pl17po14 pl17po22 pl18po1  pl20po2  pl20po10 pl20po16 pl21po0  pl312po6)
potsMgCl=$(echo pl22po10 pl22po23 pl23po0  pl23po6  pl23po20 pl24po3  pl24po5  pl27po12 pl29po5  pl29po6  pl29po20 pl30po3  pl31po6  pl34po8  pl34po18 pl35po16)
fi
if [ $strain == 19 ]
then
potsOTU5=$(echo pl2po4   pl2po20  pl3po9   pl4po1   pl4po13  pl5po2   pl5po12  pl6po0   pl8po2   pl11po0  pl11po11 pl11po12 pl11po23 pl12po17 pl13po11 pl14po6  pl15po8  pl15po11 pl16po6  pl17po15 pl18po6  pl19po20 pl20po8  pl20po21)
potsMgCl=$(echo pl24po1  pl25po9  pl25po12 pl26po5  pl26po22 pl26po23 pl27po10 pl28po16 pl29po2  pl29po9  pl30po5  pl32po19 pl33po8  pl34po10 pl34po21 pl35po20)
fi
if [ $strain == 20 ]
then
potsOTU5=$(echo pl1po14  pl2po1   pl4po2   pl4po5   pl4po14  pl5po6   pl6po5   pl6po11  pl8po0   pl8po16  pl8po22  pl9po12  pl9po22  pl11po3  pl13po14 pl13po16 pl16po16 pl17po16 pl18po15 pl20po9  pl20po19 pl21po7  pl21po16 pl21po22 )
potsMgCl=$(echo pl22po19 pl22po20 pl23po15 pl23po21 pl25po6  pl27po0  pl28po0  pl28po7  pl29po13 pl29po14 pl32po13 pl33po15 pl34po11 pl34po13 pl35po1  pl35po2 )
fi
if [ $strain == 21 ]  #DC3000
then
potsOTU5=$(echo pl1po12   pl1po17   pl2po05    pl3po07    pl3po11   pl3po22   pl4po17   pl5po23   pl8po11   pl8po12   pl10po01   pl10po02   pl10po14  pl11po18 pl12po03   pl14po00   pl15po15  pl16po23  pl17po21  pl18po00   pl20po13  pl21po10  pl21po11  pl312po10 pl312po22)  #S21
potsMgCl=$(echo pl22po12 pl24po09  pl24po19 pl24po23 pl27po16 pl27po22 pl28po05  pl28po19 pl29po08  pl29po12 pl31po10 pl31po22 pl32po02  pl34po12 pl35po00)  #S21
potsDC3000=$(echo pl38po4 pl39po11 pl39po20 pl40po15 pl41po3 pl41po8 pl41po16 pl42po2)
fi



if [ $Pseudo == MgCl ]
then
pots=$(echo $potsMgCl)
fi
if [ $Pseudo == OTU5 ]
then
pots=$(echo $potsOTU5)
fi
if [ $Pseudo == DC3000 ]
then
pots=$(echo $potsDC3000)
fi

#SEG or image
version=seg
cd /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/Phenotyping_Sphingo/Experiment1_21_strains/photos_"$date"
mkdir -p Sph"$strain"_"$Pseudo"_"$version"
rm ./Sph"$strain"_"$Pseudo"_"$version"/*
for i in $(echo $pots)
do
plate=$(sed "s/po.*//g" <(echo $i))
plate=$(sed "s/pl//g" <(echo $plate))
pot=$(sed "s/.*po//g" <(echo $i))
echo $plate $pot
#image version
    #cp /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/Phenotyping_Sphingo/Experiment1_21_strains/photos_"$date"/processed/p"$plate"_*/image/*_0\{1,2\}"$pot"_image.png Sph"$strain"_"$Pseudo"
#seg version
    cp /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/Phenotyping_Sphingo/Experiment1_21_strains/photos_"$date"/processed/p"$plate"_*/"$version"/*_0"$pot"_"$version".png Sph"$strain"_"$Pseudo"_"$version"
	cp /ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/Phenotyping_Sphingo/Experiment1_21_strains/photos_"$date"/processed/p"$plate"_*/"$version"/*_00"$pot"_"$version".png Sph"$strain"_"$Pseudo"_"$version"
done


cd Sph"$strain"_"$Pseudo"_"$version"
montage \
    -geometry 80x80+0+0 \
    -label '%t' \
    -pointsize 8 \
    *.png \
    -frame 0 \
    -tile 8x \
    Sph"$strain"_"$Pseudo"_"$version"_overview.png


#scp dlundberg@burrito.eb.local:/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/Phenotyping_Sphingo/Experiment1_21_strains/photos_20190506_new/S*_SEG/*overview.png /Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/random_design_tables/day6_photo_montage 
scp dlundberg@burrito.eb.local:/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/Phenotyping_Sphingo/Experiment1_21_strains/photos_20190506_new/S*/*overview.png /Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/random_design_tables/day6_photo_montage 

scp dlundberg@burrito.eb.local:/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/Phenotyping_Sphingo/Experiment1_21_strains/photos_20190506_new/S*/*DC3000*overview.png /Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/random_design_tables/day6_photo_montage 


scp dlundberg@burrito.eb.local:/ebio/abt6_projects7/bacterial_strain_analysis/Pratchaya/Phenotyping_Sphingo/Experiment1_21_strains/photos_20190506_new/Sph6_OTU5_seg/*overview.png /Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Pratchaya/random_design_tables/day6_photo_montage

