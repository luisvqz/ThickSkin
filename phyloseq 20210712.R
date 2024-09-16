#install.packages("vegan")
#install.packages("tidyverse")
#install.packages("ggplot2")

#R version 3.6.0
#install.packages("multcompView")
#install.packages("stringr")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("phyloseq")
#install.packages("remotes")
#BiocManager::install("DESeq2", force=T)
#BiocManager::install("genefilter", force=T)
#install.packages("remotes")
#remotes::install_github("twbattaglia/btools")library(btools)
#remotes::install_github("vmikk/metagMisc")

#####Calling library####

library(phyloseq)
library(ggplot2)
library(dplyr)
library(multcompView)
library(stringr)
library(vegan)
library(btools)
library(tidyverse)
library(metagMisc)

###seting WD#####

setwd("C:/Users/Luis/Google Drive/PHD/Long_biom/phyloseq20210607") #for lenovo
#setwd("D:/Gdrive/PHD/Long_biom/phyloseq20210607") #for desktop


#####importing into phyloseq########
bnf_map <-import_qiime_sample_data ("phyloseq_req/mapfile20210712.txt")
bnf_tree <- read_tree_greengenes ("phyloseq_req/tree138-21.nwk")
longbiome <-import_biom("phyloseq_req/long_biom_2021_taxonomy_2020.biom")
longbiome
#merge into new object
batsandflies<- merge_phyloseq (longbiome, bnf_map, bnf_tree)
batsandflies
sample_data(batsandflies)
#remove un-used IDs
batsandflies <- batsandflies %>%
  subset_samples(batid!="remove")
batsandflies
sample_sum_aa<- data.frame(sum = sample_sums(batsandflies))
sample_sum_aa
colSums(sample_sum_aa) #total number of reads prefiltered

#remove extraction and PCR controls
batsandflies<- batsandflies%>%
  subset_samples(batid!="control")

###Renaming taxonomy####
#Change label for taxonomy
tax <- data.frame(tax_table(batsandflies))
tax
tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "d__",""),
                        Phylum = str_replace(tax[,2], "p__",""),
                        Class = str_replace(tax[,3], "c__",""),
                        Order = str_replace(tax[,4], "o__",""),
                        Family = str_replace(tax[,5], "f__",""),
                        Genus = str_replace(tax[,6], "g__",""),
                        Species = str_replace(tax[,7], "s__",""),
                        stringsAsFactors = FALSE)
tax.clean[is.na(tax.clean)] <- ""
tax.clean[1:20,1:7]
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}

####### Filling the NA's in the tax table
tax.clean[is.na(tax.clean)] <- ""

#### re import as matrix into the S4 object
tax_table(batsandflies) <- as.matrix(tax.clean)

tax_table(batsandflies)[1:100,1:7]

#change "NA" for ""
batsandflies = subset_taxa(batsandflies, Phylum!="")
batsandflies

#Dataframe for ASVs counts
sample_sum_df<- data.frame(sum = sample_sums(batsandflies))
sample_sum_df

#plot for ASV counts
asvcounts<-ggplot(sample_sum_df, aes(x = sum)) + geom_histogram()
asvcounts+
  ggtitle("Distribution of sample sequencing depth")+
  xlab("Read counts")+
  ylab("Frequency")+theme_classic(base_size = 16)
rank_names(batsandflies)
#safety cutoff for crosschecking
bandb <- batsandflies %>%
  subset_samples(Type!="Mock")
####top 10 Phyla #####

phylum.sum_sub = tapply(taxa_sums(bandb), tax_table(bandb)[, "Phylum"], sum, na.rm=TRUE)
phylum.sum_sub.top10=names(sort(phylum.sum_sub, TRUE))[1:10]
phylum.sum_sub.top10
bandb #number of samples and ASvs
sample_sum_bd<- data.frame(sum = sample_sums(bandb))
sample_sum_bd
colSums(sample_sum_bd) #total of reads
summary(sample_sum_bd$sum) #mean number of reads pewr sample
sd(sample_sum_bd$sum) #SD of reads per sample

#for removing anything below 8000 reads
bandb<-prune_samples(sample_sums(bandb)>=8000,bandb)
bandb
#####Calculating Alpha diversity Measurements#######
Richness<-estimate_richness(bandb,measures=c("Observed","Shannon"))
FaithsPD<-estimate_pd(bandb)
names(Richness)
sample_data(bandb)
Alpha<-sample_data(bandb)
names(Alpha)

Alpha$Observed<-Richness$Observed
Alpha$Shannon<-Richness$Shannon
Alpha$FPD<-FaithsPD$PD
head(Alpha)

#Lets also add sequencing depth per sample:
Alpha$SequencingDepth<-sample_sums(bandb)
head(Alpha)
str(Alpha)
write.csv(Alpha, "alphadiversity.csv")
#export as data frame
Alphadf<-data.frame(sample_data(Alpha))
#reassign as factors
Alphadf$Type2<-as.factor(Alphadf$Type2)
Alphadf$Locality<-as.factor(Alphadf$Locality)
#we remove the single cave point from each locality
Alphadf<- Alphadf %>%
  filter(Type2!="Cave")
#drop the unused levels
Alphadf<- droplevels.data.frame(Alphadf)
levels(Alphadf$Type2)
#subset in localities
Alphadfchamela<- Alphadf %>%
  filter(Locality=="Chamela")
Alphadfcoquimatlan<- Alphadf %>%
  filter(Locality=="Coquimatlán")

##### Stats for Alpha####
#kruskal wallis and pairwise tests
#FPD
chamela_kw_PD<-kruskal.test(FPD~Type2, data=Alphadfchamela)
chamela_kw_PD
pair_chamela_kw_PD<-pairwise.wilcox.test(Alphadfchamela$FPD, Alphadfchamela$Type2,
                     p.adjust.method = "BH")
pair_chamela_kw_PD
coquimatlan_kw_PD<-kruskal.test(FPD~Type2, data=Alphadfcoquimatlan)
coquimatlan_kw_PD
pair_coquimatlan_kw_PD<-pairwise.wilcox.test(Alphadfcoquimatlan$FPD, Alphadfcoquimatlan$Type2,
                                             p.adjust.method = "BH")
pair_coquimatlan_kw_PD
#Shannon
chamela_kw_Shanon<-kruskal.test(Shannon~Type2, data=Alphadfchamela)
chamela_kw_Shanon
pair_chamela_kw_Shannon<-pairwise.wilcox.test(Alphadfchamela$Shannon, Alphadfchamela$Type2,
                                         p.adjust.method = "BH")
pair_chamela_kw_Shannon
coquimatlan_kw_Shanon<-kruskal.test(Shannon~Type2, data=Alphadfcoquimatlan)
coquimatlan_kw_Shanon
pair_coquimatlan_kw_Shannon<-pairwise.wilcox.test(Alphadfcoquimatlan$Shannon, Alphadfcoquimatlan$Type2,
                                                  p.adjust.method = "BH")
pair_coquimatlan_kw_Shannon
#ASVs
chamela_kw_ASV<-kruskal.test(Observed~Type2, data=Alphadfchamela)
chamela_kw_ASV
pair_chamela_kw_Observed<-pairwise.wilcox.test(Alphadfchamela$Observed, Alphadfchamela$Type2,
                                         p.adjust.method = "BH")
pair_chamela_kw_Observed
coquimatlan_kw_ASV<-kruskal.test(Observed~Type2, data=Alphadfcoquimatlan)
coquimatlan_kw_ASV
pair_coquimatlan_kw_Observed<-pairwise.wilcox.test(Alphadfcoquimatlan$Observed, Alphadfcoquimatlan$Type2,
                                                   p.adjust.method = "BH")
pair_coquimatlan_kw_Observed


###################################Colors#####################################################################

###########color pallets

colas<-c("Chamela"="#c7d42d",
         "Coquimatlán"="#0096ee",
         "positive"="light grey",
         "Bat GI"="#8a73c9",
         "Bat Skin"="#87a141",
        "Cave" ="#ca5686",
        "Nycterophilia" ="#4aac88",
         "Trichobius" ="#ca7040",
        "F"="red",
        "M"="blue")
colasbeta<-c("Bat GI"="#8a73c9",
             "Bat Skin"="#87a141",
             "Cave" ="#ca5686",
             "Nycterophilia" ="#4aac88",
             "Trichobius" ="#ca7040")
Taxacolors<-c("Anaplasmataceae"="#d48229",
             "Bacillaceae"="#7c63d1",
             "Chitinophagaceae"="#4dc461",
             "Clostridiaceae"="#bf4fb4",
             "Corynebacteriaceae"="#67aa37",
             "Enterobacteriaceae"="#da4a82",
             "Flavobacteriaceae"="#62c396",
             "Halomonadaceae"="#d2424f",
             "Morganellaceae"="#48adcf",
             "Mycobacteriaceae"="#ce4d2a",
             "Mycoplasmataceae"="#6b7cc4",
             "Pasteurellaceae"="#adb649",
             "Pseudonocardiaceae"="#d18ac7",
             "Rhizobiaceae"="#5a8237",
             "Rickettsiaceae"="#a04b6b",
             "Salinisphaeraceae"="#36845f",
             "Saprospiraceae"="#c7725a",
             "Staphylococcaceae"="#d2a652",
             "Oxalobacteraceae"="#c23917",
             "Comamonadaceae"="#937be0",
             "Others (< 5% abund.)"="#807f80")
genuscolors<-c("Bartonella"="#9d6cc1",
               "Noviherbaspirillum"="#3fadaf",
               "Wolbachia"="#bb873c",
               "Arsenophonus"="#6ea84e",
               "endosymbionts"="#cb5362",
               "Others (< 5% abund.)"="#807f80")

###THEMES for plots #####
compositiontheme<-theme_bw()+theme(axis.title.x = element_blank(),axis.title.y = element_text(face="bold", size=14),
   axis.text.x=element_text(colour="black", face="bold", size=14, hjust = 0.5), 
   axis.text.y=element_text(colour="black", size = 14),legend.title = element_blank(),legend.text = element_text(size = 14),legend.position='right',
   axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())
themebloc<-theme_minimal()+theme(axis.title.x = element_text(color = "black", face="bold", size = 16),axis.title.y = element_text(color="black",face="bold", size=16),  
   axis.text.x=element_text(colour="black", face="bold", size=14,angle = 0, hjust = 0.5),axis.text.y=element_text(colour="black", size = 14),
   title = element_text(color = "black", face="bold", size = 16),legend.text = element_text(size=12),legend.position='none',strip.text.y=element_text(colour = "Black", face="bold", size = 14),
   axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())
themebetaplots<-theme_bw() +theme(axis.title.x = element_text(face="bold", size=14),axis.title.y = element_text(face="bold", size=14),axis.text.x=element_text(colour="black", face="bold", size=14, hjust = 0.5),
  axis.text.y=element_text(colour="black", size = 14),legend.title = element_blank(),legend.text = element_text(size = 10),
  legend.position='top',axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())

#####Plotting Alpha diversity####

Alphacaveless<- Alpha%>%
  subset_samples(Type2!="Cave")

#Shannon

Shanbloca<-ggplot(data=Alphacaveless,aes(x=Type2, y=Shannon, fill=Type2))+ggtitle("Shannon diversity by locality and sample type")+facet_grid(rows=vars(Locality))
Shanbloca+geom_boxplot()+geom_jitter(width = 0.25)+ylab("Shannon Index")+scale_fill_manual(values=colas)+themebloc
ggsave("Shannon boxplot all types by Locality.png",dpi = 300, units = "in", height = 5, width = 6)

#Observed ASVs

ASVbloca<-ggplot(data=Alphacaveless,aes(x=Type2, y=Observed, fill=Type2))+scale_fill_manual(values=colas)+ggtitle("Observed ASVs by locality and per type")+facet_grid(rows=vars(Locality))
ASVbloca+geom_boxplot()+geom_jitter(width = 0.25)+ylab("Number of ASVs")+themebloc

ggsave("Observed ASVs boxplot all types by Locality.png",dpi = 300, units = "in", height = 5, width = 6)

#Faith's PD

FPDbloca<-ggplot(data=Alphacaveless,aes(x=Type2, y=FPD, fill=Type2))+ggtitle("Faith's PD by locality and per type")+facet_grid(rows=vars(Locality))
FPDbloca+geom_boxplot()+geom_jitter(width = 0.25)+ylab("Faith's PD")+scale_fill_manual(values=colas)+themebloc
ggsave("Faiths PD boxplot all types by Locality.png",dpi = 300, units = "in", height = 5, width = 6)


###### Composition plots######
melt_Family_loc_type <- bandb %>%
  tax_glom(taxrank = "Family") %>%      
  psmelt() %>%                         
  arrange(Family)                     

melt_Family_loc_type_chamela<-melt_Family_loc_type %>%
  filter(Locality=="Chamela")
melt_Family_loc_type_coquimatlan<-melt_Family_loc_type %>%
  filter(Locality=="Coquimatlán")


melt_Family_loc_type_chamela2 <- aggregate(Abundance ~ Family + Type2, 
                                   data= melt_Family_loc_type_chamela, 
                                   sum)
melt_Family_loc_type_coquimatlan2 <- aggregate(Abundance ~ Family + Type2, 
                                   data= melt_Family_loc_type_coquimatlan, 
                                   sum)

melt_Family_loc_type_chamela2$Family <- as.character(melt_Family_loc_type_chamela2$Family) 
melt_Family_loc_type_coquimatlan2$Family <- as.character(melt_Family_loc_type_coquimatlan2$Family) 

melt_Family_loc_type_chamela3 <- melt_Family_loc_type_chamela2 %>% 
  dplyr::group_by(Type2) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq > 0.05)

melt_Family_loc_type_coquimatlan3 <- melt_Family_loc_type_coquimatlan2 %>% 
  dplyr::group_by(Type2) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq > 0.05)

#Get Remainers to add in graph
remainers_Family_loc_type_chamela <- melt_Family_loc_type_chamela2 %>% 
  dplyr::group_by(Type2) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq < 0.05)

remainers_Family_loc_type_chamela2 <- aggregate(rel.freq ~ Type2, 
                                                data= remainers_Family_loc_type_chamela, 
                                                sum)

remainers_Family_loc_type_chamela2$Family <- "Others (< 5% abund.)"

join_Family_loc_type_chamela <- full_join(melt_Family_loc_type_chamela3,remainers_Family_loc_type_chamela2)

join_Family_loc_type_chamela <- join_Family_loc_type_chamela %>%
  dplyr::mutate(percent = paste0(round(100 * rel.freq, 0), "%"))

remainers_Family_loc_type_coquimatlan <- melt_Family_loc_type_coquimatlan2 %>% 
  dplyr::group_by(Type2) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq < 0.05)

remainers_Family_loc_type_coquimatlan2 <- aggregate(rel.freq ~ Type2, 
                                                    data= remainers_Family_loc_type_coquimatlan, 
                                                    sum)

remainers_Family_loc_type_coquimatlan2$Family <- "Others (< 5% abund.)"

join_Family_loc_type_coquimatlan <- full_join(melt_Family_loc_type_coquimatlan3,remainers_Family_loc_type_coquimatlan2)

join_Family_loc_type_coquimatlan <- join_Family_loc_type_coquimatlan %>%
  dplyr::mutate(percent = paste0(round(100 * rel.freq, 0), "%"))

join_Family_loc_type_chamela$Family <- as.factor(join_Family_loc_type_chamela$Family)
join_Family_loc_type_chamela$Family <- reorder(join_Family_loc_type_chamela$Family, join_Family_loc_type_chamela$Abundance)
family_species_loc_type_chamela <- ggplot(join_Family_loc_type_chamela, aes(x = Type2, y = rel.freq, fill = Family)) + 
  geom_bar(stat = "identity", position = "fill") +
  ylab("Relative abundance") +
  xlab("Type2") +ggtitle("composition by sample type in Chamela")+ scale_fill_manual(values=Taxacolors)+ compositiontheme

family_species_loc_type_chamela + geom_text(aes(label=rel.freq), position=position_fill(vjust=0.5), size=7)
ggsave("Chamela Overall Type2 composition FAMILY.png",dpi = 300, units = "in", height = 5, width = 12)

join_Family_loc_type_coquimatlan$Family <- as.factor(join_Family_loc_type_coquimatlan$Family)
join_Family_loc_type_coquimatlan$Family <- reorder(join_Family_loc_type_coquimatlan$Family, join_Family_loc_type_coquimatlan$Abundance)
family_species_loc_type_coquimatlan <- ggplot(join_Family_loc_type_coquimatlan, aes(x = Type2, y = rel.freq, fill = Family)) + 
  geom_bar(stat = "identity", position = "fill") +
  ylab("Relative abundance") +
  xlab("Type2") +ggtitle("composition by sample type in Coquimatlán")+ scale_fill_manual(values=Taxacolors)+ compositiontheme
  
family_species_loc_type_coquimatlan+ geom_text(aes(label=rel.freq), position=position_fill(vjust=0.5), size=7)

ggsave("coquimatlan Overall Type2 composition FAMILY.png",dpi = 300, units = "in", height = 5, width = 12)


bandbfly<- bandb %>%
  subset_samples(BatFly=="Y")

melt_Family_loc_type_batfly <- bandbfly %>%
  tax_glom(taxrank = "Family") %>%      
  psmelt() %>%                         
  arrange(Family)                     

melt_Family_loc_type_batfly_chamela<-melt_Family_loc_type_batfly %>%
  filter(Locality=="Chamela")
melt_Family_loc_type_batfly_coquimatlan<-melt_Family_loc_type_batfly %>%
  filter(Locality=="Coquimatlán")


melt_Family_loc_type_batfly_chamela2 <- aggregate(Abundance ~ Family + Type2, 
                                                  data= melt_Family_loc_type_batfly_chamela, 
                                                  sum)
melt_Family_loc_type_batfly_coquimatlan2 <- aggregate(Abundance ~ Family + Type2, 
                                                      data= melt_Family_loc_type_batfly_coquimatlan, 
                                                      sum)

melt_Family_loc_type_batfly_chamela2$Family <- as.character(melt_Family_loc_type_batfly_chamela2$Family) 
melt_Family_loc_type_batfly_coquimatlan2$Family <- as.character(melt_Family_loc_type_batfly_coquimatlan2$Family) 

melt_Family_loc_type_batfly_chamela3 <- melt_Family_loc_type_batfly_chamela2 %>% 
  dplyr::group_by(Type2) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq > 0.05)

melt_Family_loc_type_batfly_coquimatlan3 <- melt_Family_loc_type_batfly_coquimatlan2 %>% 
  dplyr::group_by(Type2) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq > 0.05)

#Get Remainers to add in graph
remainers_Family_loc_type_batfly_chamela <- melt_Family_loc_type_batfly_chamela2 %>% 
  dplyr::group_by(Type2) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq < 0.05)

remainers_Family_loc_type_batfly_chamela2 <- aggregate(rel.freq ~ Type2, 
                                                       data= remainers_Family_loc_type_batfly_chamela, 
                                                       sum)

remainers_Family_loc_type_batfly_chamela2$Family <- "Others (< 5% abund.)"

join_Family_loc_type_batfly_chamela <- full_join(melt_Family_loc_type_batfly_chamela3,remainers_Family_loc_type_batfly_chamela2)

join_Family_loc_type_batfly_chamela <- join_Family_loc_type_batfly_chamela %>%
  dplyr::mutate(percent = paste0(round(100 * rel.freq, 0), "%"))

remainers_Family_loc_type_batfly_coquimatlan <- melt_Family_loc_type_batfly_coquimatlan2 %>% 
  dplyr::group_by(Type2) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq < 0.05)

remainers_Family_loc_type_batfly_coquimatlan2 <- aggregate(rel.freq ~ Type2, 
                                                           data= remainers_Family_loc_type_batfly_coquimatlan, 
                                                           sum)

remainers_Family_loc_type_batfly_coquimatlan2$Family <- "Others (< 5% abund.)"

join_Family_loc_type_batfly_coquimatlan <- full_join(melt_Family_loc_type_batfly_coquimatlan3,remainers_Family_loc_type_batfly_coquimatlan2)

join_Family_loc_type_batfly_coquimatlan <- join_Family_loc_type_batfly_coquimatlan %>%
  dplyr::mutate(percent = paste0(round(100 * rel.freq, 0), "%"))

join_Family_loc_type_batfly_chamela$Family <- as.factor(join_Family_loc_type_batfly_chamela$Family)
join_Family_loc_type_batfly_chamela$Family <- reorder(join_Family_loc_type_batfly_chamela$Family, join_Family_loc_type_batfly_chamela$Abundance)
Family_species_loc_type_batfly_chamela <- ggplot(join_Family_loc_type_batfly_chamela, aes(x = Type2, y = rel.freq, fill = Family)) + 
  geom_bar(stat = "identity", position = "fill") +
  ylab("Relative abundance") +
  xlab("Type2") +ggtitle("composition by batfly Family in Chamela")+ scale_fill_manual(values=Taxacolors)+ compositiontheme

Family_species_loc_type_batfly_chamela + geom_text(aes(label=rel.freq), position=position_fill(vjust=0.5), size=7)

ggsave("Chamela batfly composition Family.png",dpi = 300, units = "in", height = 5, width = 12)

join_Family_loc_type_batfly_coquimatlan$Family <- as.factor(join_Family_loc_type_batfly_coquimatlan$Family)
join_Family_loc_type_batfly_coquimatlan$Family <- reorder(join_Family_loc_type_batfly_coquimatlan$Family, join_Family_loc_type_batfly_coquimatlan$Abundance)
Family_species_loc_type_batfly_coquimatlan <- ggplot(join_Family_loc_type_batfly_coquimatlan, aes(x = Type2, y = rel.freq, fill = Family)) + 
  geom_bar(stat = "identity", position = "fill") +
  ylab("Relative abundance") +
  xlab("Type2") +ggtitle("composition by batfly Family in Coquimatlán")+ scale_fill_manual(values=Taxacolors)+ compositiontheme

Family_species_loc_type_batfly_coquimatlan+ geom_text(aes(label=rel.freq), position=position_fill(vjust=0.5), size=7)

ggsave("coquimatlan batfly composition Family.png",dpi = 300, units = "in", height = 5, width = 12)
melt_Genus_loc_type_batfly <- bandbfly %>%
  tax_glom(taxrank = "Genus") %>%      
  psmelt() %>%                         
  arrange(Genus)                     

melt_Genus_loc_type_batfly_chamela<-melt_Genus_loc_type_batfly %>%
  filter(Locality=="Chamela")
melt_Genus_loc_type_batfly_coquimatlan<-melt_Genus_loc_type_batfly %>%
  filter(Locality=="Coquimatlán")


melt_Genus_loc_type_batfly_chamela2 <- aggregate(Abundance ~ Genus + Type2, 
                                                 data= melt_Genus_loc_type_batfly_chamela, 
                                                 sum)
melt_Genus_loc_type_batfly_coquimatlan2 <- aggregate(Abundance ~ Genus + Type2, 
                                                     data= melt_Genus_loc_type_batfly_coquimatlan, 
                                                     sum)

melt_Genus_loc_type_batfly_chamela2$Genus <- as.character(melt_Genus_loc_type_batfly_chamela2$Genus) 
melt_Genus_loc_type_batfly_coquimatlan2$Genus <- as.character(melt_Genus_loc_type_batfly_coquimatlan2$Genus) 

melt_Genus_loc_type_batfly_chamela3 <- melt_Genus_loc_type_batfly_chamela2 %>% 
  dplyr::group_by(Type2) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq > 0.05)

melt_Genus_loc_type_batfly_coquimatlan3 <- melt_Genus_loc_type_batfly_coquimatlan2 %>% 
  dplyr::group_by(Type2) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq > 0.05)

#Get Remainers to add in graph
remainers_Genus_loc_type_batfly_chamela <- melt_Genus_loc_type_batfly_chamela2 %>% 
  dplyr::group_by(Type2) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq < 0.05)

remainers_Genus_loc_type_batfly_chamela2 <- aggregate(rel.freq ~ Type2, 
                                                      data= remainers_Genus_loc_type_batfly_chamela, 
                                                      sum)

remainers_Genus_loc_type_batfly_chamela2$Genus <- "Others (< 5% abund.)"

join_Genus_loc_type_batfly_chamela <- full_join(melt_Genus_loc_type_batfly_chamela3,remainers_Genus_loc_type_batfly_chamela2)

join_Genus_loc_type_batfly_chamela <- join_Genus_loc_type_batfly_chamela %>%
  dplyr::mutate(percent = paste0(round(100 * rel.freq, 0), "%"))

remainers_Genus_loc_type_batfly_coquimatlan <- melt_Genus_loc_type_batfly_coquimatlan2 %>% 
  dplyr::group_by(Type2) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq < 0.05)

remainers_Genus_loc_type_batfly_coquimatlan2 <- aggregate(rel.freq ~ Type2, 
                                                          data= remainers_Genus_loc_type_batfly_coquimatlan, 
                                                          sum)

remainers_Genus_loc_type_batfly_coquimatlan2$Genus <- "Others (< 5% abund.)"

join_Genus_loc_type_batfly_coquimatlan <- full_join(melt_Genus_loc_type_batfly_coquimatlan3,remainers_Genus_loc_type_batfly_coquimatlan2)

join_Genus_loc_type_batfly_coquimatlan <- join_Genus_loc_type_batfly_coquimatlan %>%
  dplyr::mutate(percent = paste0(round(100 * rel.freq, 0), "%"))

join_Genus_loc_type_batfly_chamela$Genus <- as.factor(join_Genus_loc_type_batfly_chamela$Genus)
join_Genus_loc_type_batfly_chamela$Genus <- reorder(join_Genus_loc_type_batfly_chamela$Genus, join_Genus_loc_type_batfly_chamela$Abundance)
Genus_species_loc_type_batfly_chamela <- ggplot(join_Genus_loc_type_batfly_chamela, aes(x = Type2, y = rel.freq, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill") +
  ylab("Relative abundance") +
  xlab("Type2") +ggtitle("composition by batfly Genus in Chamela")+ scale_fill_manual(values=genuscolors)+ compositiontheme
Genus_species_loc_type_batfly_chamela + geom_text(aes(label=rel.freq), position=position_fill(vjust=0.5), size=7)

ggsave("Chamela batfly composition Genus.png",dpi = 300, units = "in", height = 4, width = 6)

join_Genus_loc_type_batfly_coquimatlan$Genus <- as.factor(join_Genus_loc_type_batfly_coquimatlan$Genus)
join_Genus_loc_type_batfly_coquimatlan$Genus <- reorder(join_Genus_loc_type_batfly_coquimatlan$Genus, join_Genus_loc_type_batfly_coquimatlan$Abundance)
Genus_species_loc_type_batfly_coquimatlan <- ggplot(join_Genus_loc_type_batfly_coquimatlan, aes(x = Type2, y = rel.freq, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill") +
  ylab("Relative abundance") +
  xlab("Type2") +ggtitle("composition by batfly Genus in Coquimatlán")+ scale_fill_manual(values=genuscolors)+ compositiontheme

Genus_species_loc_type_batfly_coquimatlan+ geom_text(aes(label=rel.freq), position=position_fill(vjust=0.5), size=7)

ggsave("coquimatlan batfly composition Genus.png",dpi = 300, units = "in", height = 4, width = 6)

#Top10genus

genus.sum_sub = tapply(taxa_sums(bandbfly), tax_table(bandbfly)[, "Genus"], sum, na.rm=TRUE)
genus.sum_sub.top20=names(sort(genus.sum_sub, TRUE))[1:20]
genus.sum_sub.top20


#BETA DIVERSITY

##### Calculating Beta Diversity####

#remove ALL singletons
bandb1<- filter_taxa(bandb, function (x) {sum(x > 0) >1}, prune=TRUE)
sample_data(bandb1)

bandbChamela<- bandb1%>%
  subset_samples(Locality=="Chamela")

bandbCoquimatlan<- bandb1%>%
  subset_samples(Locality=="Coquimatlán")

# weighted unifrac 
DistW = distance(bandb1,method="wunifrac")
DistWcha = distance(bandbChamela,method="wunifrac")
DistWcoq = distance(bandbCoquimatlan,method="wunifrac")
#unweighted unifrac 
DistUW = distance(bandb1,method="uunifrac")
DistUWcha = distance(bandbChamela,method="uunifrac")
DistUWcoq = distance(bandbCoquimatlan,method="uunifrac")

#create ordination
ordUW = ordinate(bandb1, method = "PCoA", distance = DistUW)
ordUWcha = ordinate(bandbChamela, method = "PCoA", distance = DistUWcha)
ordUWcoq = ordinate(bandbCoquimatlan, method = "PCoA", distance = DistUWcoq)
ordW = ordinate(bandb1, method = "PCoA", distance = DistW)
ordWcha = ordinate(bandbChamela, method = "PCoA", distance = DistWcha)
ordWcoq = ordinate(bandbCoquimatlan, method = "PCoA", distance = DistWcoq)

#plot ordinations by locality
BetaWallsites<-plot_ordination(bandb1, ordW, color = "Type2",axes=c(2,1))+  stat_ellipse() +  ggtitle("Weighted UniFrac by locality")+theme_minimal() + scale_color_manual(values=colasbeta)
BetaWallsites
ggsave("PCoA Weighted UniFrac by locality.png",dpi = 300, units = c("in"), height = 5, width = 8)
BetaUWallsites<-plot_ordination(bandb1, ordUW, color = "Type2",axes=c(2,1))+  stat_ellipse() +  ggtitle("UnWeighted UniFrac by locality")+theme_minimal() + scale_color_manual(values=colasbeta)
BetaUWallsites
ggsave("PCoA Un-Weighted UniFrac by locality.png",dpi = 300, units = c("in"), height = 5, width = 8)
#plot ordinations by bat sex
BetaWbatsex<-plot_ordination(bandb1, ordW, color = "Bat.Sex",axes=c(2,1))+  stat_ellipse() +  ggtitle("Weighted UniFrac by Bat.Sex")+theme_minimal() + scale_color_manual(values=colasbeta)
BetaWbatsex
ggsave("PCoA Weighted UniFrac by Bat.Sex.png",dpi = 300, units = c("in"), height = 5, width = 8)
BetaUWbatsex<-plot_ordination(bandb1, ordUW, color = "Bat.Sex",axes=c(2,1))+  stat_ellipse() +  ggtitle("UnWeighted UniFrac by Bat.Sex")+theme_minimal() + scale_color_manual(values=colasbeta)
BetaUWbatsex
ggsave("PCoA Un-Weighted UniFrac by Bat.Sex.png",dpi = 300, units = c("in"), height = 5, width = 8)

#plot ordination for Chamela
BetaUWchamela<-plot_ordination(bandbChamela, ordUWcha, color = "Type2",axes=c(2,1))+  stat_ellipse() +  ggtitle("UnWeighted UniFrac by locality")+theme_minimal() + scale_color_manual(values=colasbeta)
BetaUWchamela
BetaWchamela<-plot_ordination(bandbChamela, ordWcha, color = "Type2",axes=c(2,1))+  stat_ellipse() +  ggtitle("Weighted UniFrac by locality")+theme_minimal() + scale_color_manual(values=colasbeta)
BetaWchamela

#plot ordination for Coquimatlán
BetaUWCoquimatlan<-plot_ordination(bandbCoquimatlan, ordUWcoq, color = "Type2",axes=c(2,1))+  stat_ellipse() +  ggtitle("UnWeighted UniFrac by locality")+theme_minimal() + scale_color_manual(values=colasbeta)
BetaUWCoquimatlan
BetaWCoquimatlan<-plot_ordination(bandbCoquimatlan, ordWcoq, color = "Type2",axes=c(2,1))+  stat_ellipse() +  ggtitle("Weighted UniFrac by locality")+theme_minimal() + scale_color_manual(values=colasbeta)
BetaWCoquimatlan


### creating database
metadata_chamela<-data.frame(sample_data(bandbChamela))
metadata_chamela$Type2<-as.factor(metadata_chamela$Type2)

####Beta Diversity Centroid plots#####

W_chamela<-data.frame(ordWcha$vectors[,1],
                      ordWcha$vectors[,2])
colnames(W_chamela)[1]<-"MDS1"
colnames(W_chamela)[2]<-"MDS2"

#Weighted Chamela
Wcentroid_chamela<-cbind(W_chamela,metadata_chamela)
centroids_W_chamela <- as.data.frame(Wcentroid_chamela %>% 
                                       dplyr::group_by(Type2) %>% # calculate functions below for each group
                                       dplyr::summarise(mean_MDS1=mean(MDS1),
                                                        mean_MDS2=mean(MDS2),
                                                        n_MDS1=length(MDS1),
                                                        n_MDS2=length(MDS2),
                                                        stdv_MDS1=sd(MDS1),
                                                        stdv_MDS2=sd(MDS2),
                                                        se_MDS1=stdv_MDS1/sqrt(n_MDS1),
                                                        se_MDS2=stdv_MDS2/sqrt(n_MDS2)))
W_unifrac_plot_by_Type2_chamela <- ggplot(data = centroids_W_chamela, aes(x=mean_MDS1, y=mean_MDS2, color=Type2))+
  geom_point(size=4) +geom_errorbar(aes(ymin=mean_MDS2-se_MDS2, 
                                        ymax=mean_MDS2+se_MDS2), width=0.03,  size=0.5)+
  geom_errorbarh(aes(xmin=mean_MDS1-se_MDS1, 
                     xmax=mean_MDS1+se_MDS1), height=0.03, size=0.5) +
  geom_point(data = Wcentroid_chamela, aes(x=MDS1, y=MDS2, color=Type2), alpha=0.35, size=4) +
  scale_color_manual(values = colasbeta)+
  labs(x = "MDS1 [46.3%]", y="MDS2 [21.7%]")+
  ggtitle("Weighted UniFrac Chamela") +themebetaplots
W_unifrac_plot_by_Type2_chamela+stat_ellipse(data=Wcentroid_chamela, aes(x=MDS1, y=MDS2, color=Type2),inherit.aes = FALSE)
ggsave("W_unifrac_plot_by_Type2r_chamela.png")

#UnWeighted Chamela
UW_chamela<-data.frame(ordUWcha$vectors[,1],
                       ordUWcha$vectors[,2])
colnames(UW_chamela)[1]<-"MDS1"
colnames(UW_chamela)[2]<-"MDS2"

UWcentroid_chamela<-cbind(UW_chamela,metadata_chamela)

centroids_UW_chamela <- as.data.frame(UWcentroid_chamela %>% 
                                        dplyr::group_by(Type2) %>% # calculate functions below for each group
                                        dplyr::summarise(mean_MDS1=mean(MDS1),
                                                         mean_MDS2=mean(MDS2),
                                                         n_MDS1=length(MDS1),
                                                         n_MDS2=length(MDS2),
                                                         stdv_MDS1=sd(MDS1),
                                                         stdv_MDS2=sd(MDS2),
                                                         se_MDS1=stdv_MDS1/sqrt(n_MDS1),
                                                         se_MDS2=stdv_MDS2/sqrt(n_MDS2)))
UW_unifrac_plot_by_Type2_chamela <- ggplot(data = centroids_UW_chamela, aes(x=mean_MDS1, y=mean_MDS2, color=Type2))+
  geom_point(size=4) +geom_errorbar(aes(ymin=mean_MDS2-se_MDS2, 
                                        ymax=mean_MDS2+se_MDS2), width=0.03,  size=0.5)+
  geom_errorbarh(aes(xmin=mean_MDS1-se_MDS1, 
                     xmax=mean_MDS1+se_MDS1), height=0.03, size=0.5) +
  geom_point(data = UWcentroid_chamela, aes(x=MDS1, y=MDS2, color=Type2), alpha=0.35, size=4) +
  scale_color_manual(values = colasbeta)+
  labs(x = "MDS1 [28%]", y="MDS2 [7.9%]")+
  ggtitle("UnWeighted UniFrac Chamela") +themebetaplots
UW_unifrac_plot_by_Type2_chamela+stat_ellipse(data=UWcentroid_chamela, aes(x=MDS1, y=MDS2, color=Type2),inherit.aes = FALSE)
ggsave("UW_unifrac_plot_by_Type2r_chamela.png")

###creating the database
metadata_coquimatlan<-data.frame(sample_data(bandbCoquimatlan))
metadata_coquimatlan$Type2<-as.factor(metadata_coquimatlan$Type2)
#Weighted Coquimatlan
W_coquimatlan<-data.frame(ordWcoq$vectors[,1],
                          ordWcoq$vectors[,2])
colnames(W_coquimatlan)[1]<-"MDS1"
colnames(W_coquimatlan)[2]<-"MDS2"

Wcentroid_coquimatlan<-cbind(W_coquimatlan,metadata_coquimatlan)

centroids_W_coquimatlan <- as.data.frame(Wcentroid_coquimatlan %>% 
                                           dplyr::group_by(Type2) %>% # calculate functions below for each group
                                           dplyr::summarise(mean_MDS1=mean(MDS1),
                                                            mean_MDS2=mean(MDS2),
                                                            n_MDS1=length(MDS1),
                                                            n_MDS2=length(MDS2),
                                                            stdv_MDS1=sd(MDS1),
                                                            stdv_MDS2=sd(MDS2),
                                                            se_MDS1=stdv_MDS1/sqrt(n_MDS1),
                                                            se_MDS2=stdv_MDS2/sqrt(n_MDS2)))
W_unifrac_plot_by_Type2_coquimatlan <- ggplot(data = centroids_W_coquimatlan, aes(x=mean_MDS1, y=mean_MDS2, color=Type2))+
  geom_point(size=4) +geom_errorbar(aes(ymin=mean_MDS2-se_MDS2, 
                                        ymax=mean_MDS2+se_MDS2), width=0.03,  size=0.5)+
  geom_errorbarh(aes(xmin=mean_MDS1-se_MDS1, 
                     xmax=mean_MDS1+se_MDS1), height=0.03, size=0.5) +
  geom_point(data = Wcentroid_coquimatlan, aes(x=MDS1, y=MDS2, color=Type2), alpha=0.35, size=4) +
  scale_color_manual(values = colasbeta)+
  labs(x = "MDS1 [86.9%]", y="MDS2 [5.5%]")+
  ggtitle("Weighted UniFrac Coquimatlán") +themebetaplots
W_unifrac_plot_by_Type2_coquimatlan+stat_ellipse(data=Wcentroid_coquimatlan, aes(x=MDS1, y=MDS2, color=Type2),inherit.aes = FALSE)
ggsave("W_unifrac_plot_by_Type2r_coquimatlan.png")

## UW Coquimatlan
UW_coquimatlan<-data.frame(ordUWcoq$vectors[,1],
                           ordUWcoq$vectors[,2])
colnames(UW_coquimatlan)[1]<-"MDS1"
colnames(UW_coquimatlan)[2]<-"MDS2"

UWcentroid_coquimatlan<-cbind(UW_coquimatlan,metadata_coquimatlan)

centroids_UW_coquimatlan <- as.data.frame(UWcentroid_coquimatlan %>% 
                                            dplyr::group_by(Type2) %>% # calculate functions below for each group
                                            dplyr::summarise(mean_MDS1=mean(MDS1),
                                                             mean_MDS2=mean(MDS2),
                                                             n_MDS1=length(MDS1),
                                                             n_MDS2=length(MDS2),
                                                             stdv_MDS1=sd(MDS1),
                                                             stdv_MDS2=sd(MDS2),
                                                             se_MDS1=stdv_MDS1/sqrt(n_MDS1),
                                                             se_MDS2=stdv_MDS2/sqrt(n_MDS2)))
UW_unifrac_plot_by_Type2_coquimatlan <- ggplot(data = centroids_UW_coquimatlan, aes(x=mean_MDS1, y=mean_MDS2, color=Type2))+
  geom_point(size=4) +geom_errorbar(aes(ymin=mean_MDS2-se_MDS2, 
                                        ymax=mean_MDS2+se_MDS2), width=0.03,  size=0.5)+
  geom_errorbarh(aes(xmin=mean_MDS1-se_MDS1, 
                     xmax=mean_MDS1+se_MDS1), height=0.03, size=0.5) +
  geom_point(data = UWcentroid_coquimatlan, aes(x=MDS1, y=MDS2, color=Type2), alpha=0.35, size=4) +
  scale_color_manual(values = colasbeta)+
  labs(x = "MDS1 [29.6%]", y="MDS2 [10.4%]")+
  ggtitle("UnWeighted UniFrac Coquimatlán") +themebetaplots
UW_unifrac_plot_by_Type2_coquimatlan+stat_ellipse(data=UWcentroid_coquimatlan, aes(x=MDS1, y=MDS2, color=Type2),inherit.aes = FALSE)
ggsave("UW_unifrac_plot_by_Type2r_coquimatlan.png")

####PERMANOVAs for Beta Diversity####
perma_W_cha<-adonis(DistWcha~Type2, data=metadata_chamela, permutations = 9999)
perma_W_cha
perma_UW_cha<-adonis(DistUWcha~Type2, data=metadata_chamela, permutations = 9999)
perma_UW_cha
perma_W_coq<-adonis(DistWcoq~Type2, data=metadata_coquimatlan, permutations = 9999)
perma_W_coq
perma_UW_coq<-adonis(DistUWcoq~Type2, data=metadata_coquimatlan, permutations = 9999)
perma_UW_coq

#### BETADISP2 for Beta Diversity####
##renaming levels
levels(metadata_chamela$Type2)<-c("BG","BS","Ca","Ny","Tc")
levels(metadata_chamela$Type2)
#coquimatlan
levels(metadata_coquimatlan$Type2)<-c("BG","BS","Ca","Ny","Tc")
levels(metadata_coquimatlan$Type2)

#W Chamela
betadips2_W_chamela <- betadisper(DistWcha, metadata_chamela$Type2, type = "centroid")
betadips2_W_chamela
permutest(betadips2_W_chamela, permutations = 9999)
boxplot(betadips2_W_chamela)
plot(betadips2_W_chamela)
TU_betadips2_W_chamela<-TukeyHSD(betadips2_W_chamela)
TU_betadips2_W_chamela
boxplot(TU_betadips2_W_chamela$group)

#UW Chamela
betadips2_UW_chamela <- betadisper(DistUWcha, metadata_chamela$Type2, type = "centroid")
betadips2_UW_chamela
permutest(betadips2_UW_chamela, permutations = 9999)
boxplot(betadips2_UW_chamela)
plot(betadips2_UW_chamela)
TU_betadips2_UW_chamela<-TukeyHSD(betadips2_UW_chamela)
TU_betadips2_UW_chamela
plot(TU_betadips2_UW_chamela)


#W Coquimatlan
betadips2_W_coquimatlan <- betadisper(DistWcoq, metadata_coquimatlan$Type2, type = "centroid")
betadips2_W_coquimatlan
permutest(betadips2_W_coquimatlan, permutations = 9999)
boxplot(betadips2_W_coquimatlan)
plot(betadips2_W_coquimatlan)
TU_betadips2_W_coquimatlan<-TukeyHSD(betadips2_W_coquimatlan)
TU_betadips2_W_coquimatlan
plot(TU_betadips2_W_coquimatlan)

#UW Coquimatlan
betadips2_UW_coquimatlan <- betadisper(DistUWcoq, metadata_coquimatlan$Type2, type = "centroid")
betadips2_UW_coquimatlan
permutest(betadips2_UW_coquimatlan, permutations = 9999)
boxplot(betadips2_UW_coquimatlan)
plot(betadips2_UW_coquimatlan)
TU_betadips2_UW_coquimatlan<-TukeyHSD(betadips2_UW_coquimatlan)
TU_betadips2_UW_coquimatlan
plot(TU_betadips2_UW_coquimatlan)



