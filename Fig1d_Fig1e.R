#correlation


FOCUS=c("Sphingomonas", "Pseudomonas", "Chloroplast", "all")[3]

source("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/scripts/R/functions/microbiome_custom_functions.R")
library(vegan)
library(reshape)
library(RColorBrewer)
library(pheatmap)
library(viridis)
library(gplots)
library(dendextend)
library(magrittr)

#install.packages("dendextend")
#install.packages("magrittr")

date=format(Sys.Date(), format="%Y%m%d")


#16S
otutab="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Roger/lysates_16S/MiSeq144/MiSeq144_all_otutab_raw.txt"
tax="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Roger/lysates_16S/MiSeq144/MiSeq144_all.tax"
			#tax="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Roger/lysates_16S/MiSeq144/MiSeq144_all_no799.tax"
metadata="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Roger/lysates_16S/MiSeq144/metadataEyachMultSpec.txt"
#load OTU table
	read.countTable(file=otutab)->otureads
#load metadata
	as.matrix(read.table(metadata, sep="\t"))->metadata
	apply(cbind(metadata[,5], metadata[,7]), 1, paste, collapse="_")->codespecies
	cbind(metadata, codespecies)->metadata
#convert otu_sample names to Species Nums
   
    
    original_names=metadata[match(colnames(otureads), metadata[,1]),1]
	as.matrix(metadata[match(colnames(otureads), metadata[,1]),9])->colnames(otureads)
#load taxonomy
	as.matrix(read.table(file=tax, sep="\t"))->taxonomy

#save_FOCUS, if not all
if(FOCUS!="all"){
	remove_contaminant_taxa(otureads, taxonomy, keywords=FOCUS, invert=TRUE)->plastids
}

#remove contaminants from OTU table
	remove_contaminant_taxa(otureads, taxonomy, keywords=c("Chloroplast", "Mitochondria", "Archaea"))->otureads


otureads[,grep("Dandelion", colnames(otureads))]->Dand
otureads[,grep("Thistle", colnames(otureads))]->This

plastids[,grep("NK", colnames(otureads))]->NK

#Taraxacum sp. = Dandelion
#Trifolium sp. = clover
#Cardamine hirsuta = 
#Arabidopsis thaliana = 
#Draba verna = 
#Plantago lanceolata = 
#Impatiens glandulifera = "Himalayan balsam"




				#taxonomy[match(rownames(otureads), taxonomy[,1]),4]->tabletax
				#sum(otureads[grep("BacteriaLOST", tabletax),])/sum(otureads)
				#sum(otureads[grep("BacteriaLOST.+Pseudomonadaceae", tabletax),])/sum(otureads[grep("Bacteria.+Pseudomonadaceae", tabletax),])
				#sum(otureads[grep("BacteriaLOST.+Sphingomonadaceae", tabletax),])/sum(otureads[grep("Bacteria.+Sphingomonadaceae", tabletax),])
	
#save abundances in read #'s for OTU10 (Talia's OTU5)
otureads[c("Otu10", "Otu56", "Otu99", "Otu159"),]->Pseudos

#filter by taxonomy to ___ level
otureads<-filter_by_taxonomic_level(countTable=otureads, taxonomy=taxonomy, keywords=c("f"), last_filter=TRUE)

	
#IF YOU WANT TO INCLUDE THE FAMILY READS FROM MISEQ200, RUN "ANALYZE_MISEQ200.R" AND USE THE ASTERACEAE TABLE
#colnames(Asteraceae)=gsub("S75_Soil", "S99_Soil", colnames(Asteraceae))
#otureads=combine_taxa_tables(otureads, Asteraceae, return="both")

	
	
#order rows by average bacterial abundance
otureads[order(rowSums(otureads)),]->otureads

	
#for each sample, subtract the sum of Pseudo OTUs from Pseudomonadaceae family
split_Pseudomonas=FALSE
	if(split_Pseudomonas==TRUE){
	remaining_Pseudomonadaceae=vector(length=ncol(otureads))
	for(c in 1:ncol(otureads)){
		otureads["Pseudomonadaceae",c]-sum(Pseudos[,c])->remaining_Pseudomonadaceae[c]
		#otureads["Gammaproteobacteria",c]-sum(Pseudos[,c])->remaining_Gammaproteobacteria[c]
	}           
	otureads[which(rownames(otureads)!="Pseudomonadaceae"),]->otureads
	rbind(otureads, remaining_Pseudomonadaceae, Pseudos)->otureads
	}
	
	
	
cbind(original_names, colnames(otureads))->ENA_helper



#threshold empty samples
	otureads=otureads[,grep("empty", colnames(otureads), invert=TRUE)]
		otureads=otureads[,grep("blank", colnames(otureads), invert=TRUE)]		
	plastids=plastids[,grep("empty", colnames(plastids), invert=TRUE)]	
		plastids=plastids[,grep("blank", colnames(plastids), invert=TRUE)]	
#threshold samples 
	if(FOCUS=="Sphingomonas"){thres=100}
	if(FOCUS=="Pseudomonas"){thres=100}
	if(FOCUS=="Chloroplast"){thres=100}
	if(FOCUS=="all"){thres=300}
	otureads=otureads[,which(colSums(otureads)>=thres)]
	plastids=plastids[,which(colSums(plastids)>=thres)]
#normalize table by dividing by column sums
	normalize100(otureads)->otureads
		#t(rrarefy(t(otureads), sample=1000))->otureads
		#taxonomy[match(rownames(otureads), taxonomy[,1]),4]->fulltaxa
		
		#sum(otureads[grep("^Otu10$", rownames(otureads)),])->a
		#sum(otureads[grep("Pseudomonas", fulltaxa),])->b
		#sum(otureads[grep("Pseudomonadaceae", fulltaxa),])->c
		#a/sum(otureads)
		#a/b
		#a/c
		#b/c
		
		#sum(otureads[grep("^Otu3$", rownames(otureads)),])->A
		#sum(otureads[grep("Sphingomonas", fulltaxa),])->B
		#sum(otureads[grep("Sphingomonadaceae", fulltaxa),])->C
		#A/sum(otureads)
		#A/B
		#A/C
		#B/C
		
	normalize100(plastids)->plastids
		#t(rrarefy(t(plastids), sample=100))->plastids
		
	
if(FOCUS=="all"){otureads->plastids}
if(FOCUS!="all"){plastids->plastids}

#kill the short tree and big clover
plastids[,grep("ShortTree", colnames(plastids), invert=TRUE)]->plastids
plastids[,grep("BigClover", colnames(plastids), invert=TRUE)]->plastids

#kill these samples which are mistakes
plastids[,grep("S17_Dandelion|S43_Thistle", colnames(plastids), invert=TRUE)]->plastids

#if the focus is plastids, kill the soil
if(FOCUS=="Chloroplast"){
	plastids[,grep("Soil", colnames(plastids), invert=TRUE)]->plastids
}



ID=plastids[,grep("Moss", colnames(plastids))]
ID[which(rowSums(ID)>1),]

cbind( colnames(plastids), ENA_helper[match(colnames(plastids), ENA_helper[,2]),1] )



dendrogram=TRUE
if(dendrogram==TRUE){
vegdist(t(plastids))->plastidsdist

hClustering <- hclust(plastidsdist, method = 'complete')        
hcd = as.dendrogram(hClustering)  

DerekPalette<-c("#000000","#ffff00","#f922b9","#32fff4","#1F78B4",
				"#FB9A99","#FDBF6F","#E31A1C","#6A3D9A","#CAB2D6",
				"#FF7F00","#FFFF99","#A6CEE3", "#32fff4","#1F78B4",
				"#000000","#ffff00","#f922b9","#32fff4","#1F78B4")

gsub("Og.._", "", colnames(plastids))->plastid_simplenames
gsub("Og19.2_", "", plastid_simplenames)->plastid_simplenames
gsub("S.._", "", plastid_simplenames)->plastid_simplenames
colnames(plastids)->plastid_shapes
plastid_shapes[grep("Og.._", plastid_shapes)]<-15
plastid_shapes[grep("Og19.2_", plastid_shapes)]<-15
plastid_shapes[grep("S.._", plastid_shapes)]<-15  #15 temporarily for different display, if summer will be denoted differently
plastid_shapes[grep("blank", plastid_shapes)]<-20
as.numeric(plastid_shapes)->plastid_shapes



plant_color_list[,2]->colorCode
plant_color_list[,1]->names(colorCode)


gsub("$", " ", plastid_simplenames)->plastid_simplenames2
gsub(".+lover.+", "Trifolium", plastid_simplenames2)->plastid_simplenames2
gsub(".+Kraut.+", "Impatiens", plastid_simplenames2)->plastid_simplenames2
gsub("Dandelion", "Taraxacum", plastid_simplenames2)->plastid_simplenames2
gsub("Thistle.+", "Sonchus", plastid_simplenames2)->plastid_simplenames2
gsub(".+Tree.+", "Short_Tree", plastid_simplenames2)->plastid_simplenames2
gsub("NK", "Neckar_King", plastid_simplenames2)->plastid_simplenames2
gsub(".+thaliana.+", "Arabidopsis", plastid_simplenames2)->plastid_simplenames2

hcd %>% set("leaves_pch", plastid_shapes[order.dendrogram(hcd)]) %>%  # node point type
  set("leaves_cex", 1) %>%  # node point size
  set("leaves_col", colorCode[plastid_simplenames][order.dendrogram(hcd)]) %>% # node point color
  set("labels_cex", .5) %>%
  set("labels", plastid_simplenames2[order.dendrogram(hcd)]) %>%
  set("labels_col", colorCode[plastid_simplenames][order.dendrogram(hcd)]) %>%
  plot(main = paste("Dissimilarity of ", FOCUS, " populations based on samples with >=", thres," ", FOCUS, " seqs", sep="", collapse=""))

}

pdf(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Roger/lysates_16S/", FOCUS, "_tree_", date.2(), ".pdf", sep="", collapse=""), width = 8, height = 3, useDingbats=FALSE)
print(
hcd %>% set("leaves_pch", plastid_shapes[order.dendrogram(hcd)]) %>%  # node point type
  set("leaves_cex", .5) %>%  # node point size
  set("leaves_col", colorCode[plastid_simplenames][order.dendrogram(hcd)]) %>% # node point color
  set("labels_cex", .3) %>%
  set("labels", plastid_simplenames2[order.dendrogram(hcd)]) %>%
  set("labels_col", colorCode[plastid_simplenames][order.dendrogram(hcd)]) %>%
  plot(main = FOCUS, xaxt="n", yaxt="n")
)
dev.off()


#correct order of samples
correct_column_order=colnames(plastids)[order.dendrogram(hcd)]
plastids[,order.dendrogram(hcd)]->plastids2


#MAKE top 10 OTU HEATMAP TO GO UNDER THE DENDROGRAM
plastids2[which(rowSums(plastids2)>20),]->HEATplastids
HEATplastids[order(rowSums(HEATplastids), decreasing=TRUE),]->HEATplastids

if(FOCUS=="Chloroplast"){
	#HEATplastids[order(apply(HEATplastids, 1, max)), decreasing=TRUE),]->HEATplastids
	HEATplastids[1:100,]->HEATplastids 
}
if(FOCUS!="Chloroplast"){HEATplastids[1:10,]->HEATplastids }


#heatmap.2(sqrt(sqrt(plastids2)), trace="none", key=FALSE, dendrogram="column", margins=c(8,6), Colv = as.dendrogram(hcd), col=rev(brewer.pal(11,"Spectral")))

#OR: change this name in the code below. This is your count matrix. HEATplastids
quantile_breaks <- function(xs, n = 200) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}


#mat_breaks <- quantile_breaks(HEATplastids, n = 200)
names_mat_breaks=c("0.00%","47.24%","47.74%","48.24%","48.74%","49.25%","49.75%","50.25%","50.75%","51.26%","51.76%","52.26%","52.76%","53.27%","53.77%","54.27%","54.77%","55.28%","55.78%","56.28%","56.78%","57.29%","57.79%","58.29%","58.79%","59.30%","59.80%","60.30%","60.80%","61.31%","61.81%","62.31%","62.81%","63.32%","63.82%","64.32%","64.82%","65.33%","65.83%","66.33%","66.83%","67.34%","67.84%","68.34%","68.84%","69.35%","69.85%","70.35%","70.85%","71.36%","71.86%","72.36%","72.86%","73.37%","73.87%","74.37%","74.87%","75.38%","75.88%","76.38%","76.88%","77.39%","77.89%","78.39%","78.89%","79.40%","79.90%","80.40%","80.90%","81.41%","81.91%","82.41%","82.91%","83.42%","83.92%","84.42%","84.92%","85.43%","85.93%","86.43%","86.93%","87.44%","87.94%","88.44%","88.94%","89.45%","89.95%","90.45%","90.95%","91.46%","91.96%","92.46%","92.96%","93.47%","93.97%","94.47%","94.97%","95.48%","95.98%","96.48%","96.98%","97.49%","97.99%","98.49%","98.99%","99.50%","100.00%")
mat_breaks=as.numeric(c("0","0.01051418","0.02548292","0.02855511","0.05344371","0.06449541","0.07473077","0.08777132","0.1030792","0.13833959","0.19088079","0.20966118","0.24881384","0.29198308","0.33337395","0.3603856","0.4245307","0.43391695","0.45506257","0.49569738","0.55775115","0.60042332","0.66886278","0.70809278","0.75087294","0.77661597","0.82808326","0.90117383","0.94339623","0.97109009","1.01507331","1.10290873","1.19047619","1.24786494","1.37195436","1.57661693","1.6960936","1.75042601","1.83861885","1.87490853","1.95116549","2.00234509","2.14384399","2.26229608","2.37874388","2.46305419","2.53098997","2.70953403","2.8287865","3.21255468","3.38362304","3.52853706","3.69329246","3.97839888","4.18967134","4.54904051","4.71290972","4.85265642","5.14549079","5.30459786","5.71428571","6.0140688","6.22574077","6.54708791","6.74043293","6.92699612","7.33212813","7.54877328","7.70097556","7.9270232","8.51472019","9.58238144","10.22974349","11.14834494","11.7161639","12.34516487","13.1293813","14.54296927","15.94994151","16.88303978","18.2019155","19.33441369","21.0295789","22.22222222","23.72881356","26.05040288","28.50747933","29.37971068","30.67765968","32.85229503","35.79467583","36.94219023","37.76131741","39.8243807","41.73501001","43.88169396","49.5217753","55.71504618","59.13742699","63.53912289","65.82129328","74.1219533","84.94131056","91.4530594","93.08221363","94.40155197","99.14666667"))
names_mat_breaks->names(mat_breaks)

pdf(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Roger/lysates_16S/", FOCUS, "_tree_heatmap_", date.2(), ".pdf", sep="", collapse=""), width = 8, height = 3, useDingbats=FALSE)
print(
a<-pheatmap(HEATplastids, color=inferno(length(mat_breaks) - 1),  border_color=NA, 
         breaks=mat_breaks, cluster_rows=TRUE, cluster_cols=FALSE, fontsize_col=8, fontsize_row=8, legend=FALSE)
)
dev.off()				
pdf(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Roger/lysates_16S/", FOCUS, "_LEGEND_", date.2(), ".pdf", sep="", collapse=""), width = 8, height = 3, useDingbats=FALSE)
print(
#a<-pheatmap(HEATplastids, color=inferno(length(mat_breaks) - 1),  border_color=NA, 
 #        breaks=mat_breaks, cluster_rows=TRUE, cluster_cols=FALSE, fontsize_col=8,
  #        fontsize_row=8, legend=TRUE, legend_breaks=c(0, 0.5, 1, 5, 10, 50, 100)))

a<-pheatmap(sqrt(HEATplastids), border_color=NA, 
         cluster_rows=TRUE, cluster_cols=FALSE, fontsize_col=8,
         #color=inferno(length(mat_breaks) - 1),  
         color=inferno(20),  
          fontsize_row=8, legend=TRUE, legend_breaks=sqrt(c(0, 0.5, 1, 5, 10, 50, 80, 100)))	
)
dev.off()				
	


#low_abundance(otureads, percent=5)->otureads

gsub("Og.._", "", colnames(plastids2))->simplenames
gsub("S.._", "", simplenames)->simplenames
plastids2[c("Otu10", "Otu56", "Otu99", "Otu159") ,order(simplenames)]->Pseudos

#plot(Pseudos["Otu10",], type="l", col="green")
#points(Pseudos["Otu56",], type="l", col="darkgreen")
#points(Pseudos["Otu99",], type="l", col="blue")
#points(Pseudos["Otu159",], type="l", col="black")

#otureads[,order(simplenames)]->otureads
bottomorder=c("Otu159", "Otu99", "Otu56", "Otu10", "remaining_Pseudomonadaceae",  "Sphingomonadaceae")
toporder=c("hi")



plastids2->currentreads

#topOrder(currentreads, toporder, bottomorder)->currentreads
#order bar plot by dendrogram
#if(FOCUS!="all"){ currentreads[,order.dendrogram(hcd)]->currentreads }
#if(FOCUS=="all"){ currentreads[,order.dendrogram(hcd)]->currentreads }
#currentreads[,order(colnames(currentreads))]->currentreads}

melt(currentreads)->histDL

mean(currentreads["Pseudomonas",])
mean(currentreads["Sphingomonas",])
mean(currentreads["Sphingomonadaceae",])
mean(currentreads["Pseudomonadaceae",])


histDL[,c(2, 1, 3)]->histDL
cbind(rep(1, nrow(histDL)))->histDL$alphas
c("sample", "taxa", "abundance", "alphas")->names(histDL)


factor(histDL$taxa, levels=rownames(currentreads))->histDL$taxa
factor(histDL$sample, levels=colnames(currentreads))->histDL$sample


#SET COLOR SCHEME
rownames(currentreads)->taxa_color_pairs
COLUMN_IN_COLORSLIST="Roger"    # 3 is colors, 4 is black / white / green / purple
ALPHA_COLUMN_IN_COLORSLIST="alpha2"
0->alpha_light_taxa
if("alphalight"=="alphalight"){
	rownames(currentreads)[1:(nrow(currentreads)-20)]->alpha_light_taxa
}
source("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Regalado2017/family_colorscheme.R")
taxa_color_pairs[match(histDL$taxa, taxa_color_pairs[,1]),2]->DerekPalette
taxa_color_pairs[match(histDL$taxa, taxa_color_pairs[,1]),3]->histDL$alphas
as.numeric(histDL$alphas)->histDL$alphas

#just set anything beyond top 15 to have alpha 0.3
toptaxa=rownames(currentreads)[(length(rownames(currentreads))-14):length(rownames(currentreads))]
histDL$alphas[match(histDL[,2], toptaxa, nomatch=0)>0]=1
histDL$alphas[match(histDL[,2], toptaxa, nomatch=0)==0]=0.5

alpharange=c(0.5, 1)
PSonly=FALSE
if(PSonly==TRUE){
	DerekPalette[match(DerekPalette,c("#de2e26", "#3182BD"), nomatch=0)==0]<-"black"
	histDL$alphas[histDL$alphas!=1]=1
}

if(FOCUS=="Sphingomonas"){	
	col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
	col_vector[17:length(col_vector)]->col_vector
	"yellow"->col_vector[11]
	rep(col_vector, 50)->DerekPalette
	alpharange=1
	currentreads=currentreads[order(rowSums(currentreads)),]
}
	

#correct color methylobacteria and sphingobacteriaceae
as.matrix(rownames(currentreads)[nrow(currentreads):1][1:15][15:1])


DerekPalette[length(DerekPalette):1][1:15]->pal1
inferno(20)->pal1
image(1:length(pal1), 1, as.matrix(1:length(pal1)), 
      col=pal1, xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")



library(ggplot2)
library(reshape2)
library(RColorBrewer)
options(stringsAsFactors = FALSE)


#PLAY
pdf(paste("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Roger/lysates_16S/", FOCUS, "_stacked_bar_", date.2(), ".pdf", sep="", collapse=""), width = 8, height = 3, useDingbats=FALSE)
ggplot(histDL,aes(x=histDL$sample,y=histDL$abundance, fill=histDL$taxa,alpha=histDL$alphas)) + 
  geom_bar(position="fill", stat="identity", width=1.0) + 
  scale_fill_manual(values=DerekPalette) + theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("") + ylab("") + guides(fill=FALSE) +
  theme_classic() + theme(legend.position="none") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
	scale_alpha_continuous(range = alpharange) + 
	theme(axis.line = element_blank()) 
dev.off()				
  

#ggplot(histDL,aes(x=histDL$sample,y=histDL$abundance, fill=histDL$taxa)) + 
 # geom_bar(position="fill", stat="identity", width=1.0) + 
  #scale_fill_manual(values=DerekPalette) + 
  #xlab("") + ylab("") + theme(legend.position="bottom") + guides(fill=guide_legend(ncol=4, byrow=TRUE)) + 
  #theme_classic() +  scale_size(range=c(5,20)) +
  #theme(legend.justification=c(1,0), legend.position=c(1,0))
  
