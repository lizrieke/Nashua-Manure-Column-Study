setwd("/Users/Luby/Box Sync/column_study")

library(phyloseq)
library(ggplot2)
library(vegan)
library(RColorBrewer)
library(plyr)
library(Rmisc)

library(gridExtra)
library(DESeq2)

rare.data.taxmin5.minC85.6000.phy<-readRDS("rare.data.taxmin5.minC85.6000.phy.RDS")

# read in the sample information:
si<-read.delim("./sample_info_541samples_w_day.txt")
dim(si) # check the dimension of the sample inforamtion table
# read in the Good's estimated coverage for each sample:
C<-read.delim("data.taxmin5.sample_goods_coverage.txt")
# merge C into si:
si<-merge(si, C, by.x="sample_id", by.y="SAMPLES")
dim(si) # check merged si

# read in samples need to be resequenced:
exlude<-read.delim("samples_55_to_be_resequenced.txt")
# remove samples need to be resequenced from sample information
si.noreseq<-si[! si$id %in% exlude$unique_id, ]
dim(si.noreseq) # check new sample data's dimension

# plot the histogram of sample coverages: 
p<-ggplot(si.noreseq) + geom_histogram(aes(C, fill= sample_type, binwidth=100), position= "dodge") +facet_grid(sample_type~., scale="free_y")

p


# find the smallest C for soil column samples:
min(subset(si.noreseq, sample_type == "sc_samples")[, "C"])
# subset sample information for desired coverage:
si.noreseq.minC85 <- subset(si.noreseq, C >= "0.8489297")
dim(si.noreseq.minC85) # check dimension

#############################################
## creating phyloseq object for analyses      ##
#############################################
# read in the R ready phyloseq object:
# data.taxmin5.phy: all 541 samples, including sequencing controls, mocks, and soil column samples; otu's with sum < 5 across all samples have been removed. 
data.taxmin5.phy<-readRDS("taxsum_min5_sequence_phyloseq.RDS")
data.taxmin5.phy # check phyloseq object's content

rownames(si)<-si$sample_id # add row.names to sample information table
sample_data(data.taxmin5.phy)<-si # add sample_data to phyloseq object

# create a phyloseq object that does not contain any samples need to be resequenced or coverage less than 85%. 
data.taxmin5.minC85.phy<-subset_samples(data.taxmin5.phy, sample_id %in% si.noreseq.minC85$sample_id)
data.taxmin5.minC85.phy # check 
# remove all 0 taxa across the subsetted samples
data.taxmin5.minC85.phy<-prune_taxa(taxa_sums(data.taxmin5.minC85.phy) > 0, data.taxmin5.minC85.phy)
data.taxmin5.minC85.phy # check


#removing two bad drainage 
samplesnewPhyloObject = subset_samples(data.taxmin5.minC85.phy, id != "A020MaNtR1")
data.taxmin5.minC85.phy = subset_samples(samplesnewPhyloObject, id != "A058MaNtR2")


saveRDS(data.taxmin5.minC85.phy, "data.taxmin5.minC85.RDS")
data.taxmin5.minC85.phy<- readRDS("data.taxmin5.minC85.RDS")
sample.richness<-estimate_richness(data.taxmin5.minC85.phy)
sample.richness.matrix<-as.matrix(sample.richness)
row.names <- row.names(sample.richness.matrix)


write.table(sample.richness.matrix, "./sample.richness.matrix.txt", sep="\t", quote=F, row.names=T)




data.taxmin5.minC85.6000.phy = prune_samples(sample_sums(data.taxmin5.minC85.phy)>6000, data.taxmin5.minC85.phy)

rare.data.taxmin5.minC85.6000.phy <- rarefy_even_depth(data.taxmin5.minC85.6000.phy, sample.size=min(sample_sums(data.taxmin5.minC85.6000.phy)-1), rngseed=15879966)




##This is physloseq final object w/o manure core additions
saveRDS(rare.data.taxmin5.minC85.6000.phy, "rare.data.taxmin5.minC85.6000.phy.RDS")
rare.data.taxmin5.minC85.6000.phy<-readRDS("rare.data.taxmin5.minC85.6000.phy.RDS")


##################
########adonis
total.data.dist = phyloseq::distance(rare.data.taxmin5.minC85.6000.phy, "bray")
head(total.data.dist)
adonis(total.data.dist ~ experiment_day , perm=9999, as(sample_data(rare.data.taxmin5.minC85.6000.phy), "data.frame"))



############################comparing adonis outputs of matrix between day 24 and day 108

data.taxmin5.minC85.soil.water.phy = subset_samples(rare.data.taxmin5.minC85.6000.phy, matrix != "manure")

data.taxmin5.minC85.soil.water.manured.phy = subset_samples(data.taxmin5.minC85.soil.water.phy, manure_history == "Ma")

data.taxmin5.minC85.soil.top.water.manured.phy = subset_samples(data.taxmin5.minC85.soil.water.manured.phy, soil_layer != "B")

data.taxmin5.minC85.soil.top.water.manured.phy.d24 = subset_samples(data.taxmin5.minC85.soil.top.water.manured.phy, experiment_day == "24")

data.taxmin5.minC85.soil.top.water.manured.phy.d108 = subset_samples(data.taxmin5.minC85.soil.top.water.manured.phy, experiment_day == "108")


total.data.dist.24 = phyloseq::distance(data.taxmin5.minC85.soil.top.water.manured.phy.d24, "bray")
adonis(total.data.dist.24 ~ matrix , perm=9999, as(sample_data(data.taxmin5.minC85.soil.top.water.manured.phy.d24), "data.frame"))

total.data.dist.108 = phyloseq::distance(data.taxmin5.minC85.soil.top.water.manured.phy.d108, "bray")
adonis(total.data.dist.108 ~ matrix , perm=9999, as(sample_data(data.taxmin5.minC85.soil.top.water.manured.phy.d108), "data.frame"))



####################manure, manured soil and manured water bar chart
rare.data.taxmin5.minC85.6000.psmelt <- readRDS("rare.data.taxmin5.minC85.6000.psmelt")

rare.data.taxmin5.minC85.6000.psmelt.dply.manured<-subset(rare.data.taxmin5.minC85.6000.psmelt, manure_history !="Co")

rare.data.taxmin5.minC85.6000.psmelt.dply.manured.top.soil.water.manure<-subset(rare.data.taxmin5.minC85.6000.psmelt.dply.manured, soil_layer !="B")

rare.data.taxmin5.minC85.6000.psmelt.manured.dply <- ddply(rare.data.taxmin5.minC85.6000.psmelt.dply.manured.top.soil.water.manure, .(matrix, experiment_day, phylum), summarise, total = sum(Abundance))



write.table(rare.data.taxmin5.minC85.6000.psmelt.manured.dply, "./rare.data.taxmin5.minC85.6000.psmelt.manured.dply.txt", sep="\t", quote=F, row.names=F)

rare.data.taxmin5.minC85.6000.psmelt.manured.dply.edited<-read.delim("rare.data.taxmin5.minC85.6000.psmelt.manured.dply.edited.2.txt", header = TRUE, sep = "\t")


rare.data.taxmin5.minC85.6000.psmelt.manured.dply.edited.matrix<-subset(rare.data.taxmin5.minC85.6000.psmelt.manured.dply.edited, rare.data.taxmin5.minC85.6000.psmelt.manured.dply.edited$experiment_day!="Day 10" | rare.data.taxmin5.minC85.6000.psmelt.manured.dply.edited$experiment_day!="Day 38" |rare.data.taxmin5.minC85.6000.psmelt.manured.dply.edited$experiment_day!="Day 80")

rare.data.taxmin5.minC85.6000.psmelt.manured.dply.edited$experiment_day <- gsub("N/A", "Day 10", rare.data.taxmin5.minC85.6000.psmelt.manured.dply.edited$experiment_day)

rare.data.taxmin5.minC85.6000.psmelt.manured.dply.edited.3$experiment_day = factor(rare.data.taxmin5.minC85.6000.psmelt.manured.dply.edited.3$experiment_day, levels=c('Day 0','Day 10','Day 24','Day 38','Day 59', 'Day 80', 'Day 108'))


write.table(rare.data.taxmin5.minC85.6000.psmelt.manured.dply.edited, "./rare.data.taxmin5.minC85.6000.psmelt.manured.dply.edited.2.txt", sep="\t", quote=F, row.names=F)

rare.data.taxmin5.minC85.6000.psmelt.manured.dply.edited.2

##############version 3 combines other phyla groupings as to eliminate extra lines in the following figure
rare.data.taxmin5.minC85.6000.psmelt.manured.dply.edited.3<-read.delim("rare.data.taxmin5.minC85.6000.psmelt.manured.dply.edited.3.txt", header = TRUE, sep = "\t")
install.packages("wesanderson")
library("wesanderson")
ggplot(rare.data.taxmin5.minC85.6000.psmelt.manured.dply.edited.3, aes(x=experiment_day, y=relative_abundance, fill = Phylum))+ geom_bar(stat="identity", colour="black") + facet_grid(~matrix, scale="free_x", space="free_x")+theme_bw(20)+xlab("Experiment Day")+ylab("Relative Abundance")+theme(axis.text.x=element_text(angle = 45, vjust = 0.5))+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9","#3399FF", "#FF3300", "#CC00FF", "#66FF33", "#990000", "#003300", "#FF00CC", "#000099", "#FFFF33", "#996633", "#660066"))

ggsave("stacked.bar.chart.manured.soil.water.manure.phylum.png", width =15)



#######if needed this is how to change categories and and labels for chart
#loose.water.se$rain_event <- gsub("R1", "Day 10", nt.water.se$rain_event)
#loose.water.se$rain_event <- gsub("R2", "Day 24", nt.water.se$rain_event)
#loose.water.se$rain_event <- gsub("R3", "Day 38", nt.water.se$rain_event)
#loose.water.se$rain_event <- gsub("R4", "Day 59", nt.water.se$rain_event)
#loose.water.se$rain_event <- gsub("R5", "Day 80", nt.water.se$rain_event)
#loose.water.se$rain_event <- gsub("R6", "Day 108", nt.water.se$rain_event)
#loose.water.se$rain_event = factor(nt.water.se$rain_event, levels=c('Day 10','Day 24','Day 38','Day 59', 'Day 80', 'Day 108'))
	
	

##############nmds comparing soil man con top and bottom soil at t=0

rare.data.taxmin5.minC85.6000.phy.soil = subset_samples(rare.data.taxmin5.minC85.6000.phy, matrix == "soil")
rare.data.taxmin5.minC85.6000.phy.soil.s0 = subset_samples(rare.data.taxmin5.minC85.6000.phy.soil, soil_event == "S0")

rare.data.taxmin5.minC85.6000.phy.soil.s0.ord <- ordinate(rare.data.taxmin5.minC85.6000.phy.soil.s0, "DCA", "bray")
manure_longterm_sample_plot = plot_ordination(rare.data.taxmin5.minC85.6000.phy.soil.s0, rare.data.taxmin5.minC85.6000.phy.soil.s0.ord, type="samples", color="manure_history", shape="soil_layer")+ geom_point(size=3)+ stat_ellipse(type = "norm", linetype = 2)+theme_bw(15)
ggsave("longterm.soildca.png", width =10)

rare.data.taxmin5.minC85.6000.phy.soil.s0.psmelt <- psmelt(rare.data.taxmin5.minC85.6000.phy.soil.s0)

soil.s0.dply <- ddply(rare.data.taxmin5.minC85.6000.phy.soil.s0.psmelt, .(id, manure_history, soil_layer, phylum, class, order, family, genus), summarise, total = sum(Abundance))

ggplot(soil.s0.dply, aes(x=id, y=total, fill = phylum)) + 
	geom_bar(stat="identity", width=0.5)+facet_grid(manure_history~soil_layer)

rare.data.taxmin5.minC85.6000.phy.soil.s0.b = subset_samples(rare.data.taxmin5.minC85.6000.phy.soil.s0, soil_layer == "B")

#converting physeq object to dseq object
rare.data.taxmin5.minC85.6000.phy.soil.s0.b.dds <- phyloseq_to_deseq2(rare.data.taxmin5.minC85.6000.phy.soil.s0.b, ~ manure_history)

#run deseq
rare.data.taxmin5.minC85.6000.phy.soil.s0.b.dds.ran <- DESeq(rare.data.taxmin5.minC85.6000.phy.soil.s0.b.dds, test="Wald", fitType="parametric")


#turn deseq object into dataframe
rare.data.taxmin5.minC85.6000.phy.soil.s0.b.dds.ran.df <- data.frame(results(rare.data.taxmin5.minC85.6000.phy.soil.s0.b.dds.ran, cooksCutoff=F))
head(rare.data.taxmin5.minC85.6000.phy.soil.s0.b.dds.ran.df)
levels(rare.data.taxmin5.minC85.6000.phy.soil.s0.b.dds.ran$manure_history)

remove all OTUs w/ pvalue greater than 0.05
rare.data.taxmin5.minC85.6000.phy.soil.s0.b.dds.ran.05.df <- subset(rare.data.taxmin5.minC85.6000.phy.soil.s0.b.dds.ran.df, pvalue < 0.05)
head(rare.data.taxmin5.minC85.6000.phy.soil.s0.b.dds.ran.05.df)


#keep OTUs sig more abundant in co than ma
rar.man.taxmin5.minc85.ma.05.greater.co.df <- subset(rare.data.taxmin5.minC85.6000.phy.soil.s0.b.dds.ran.05.df, log2FoldChange < 0)

#keep OTUs sig more abundant in ma than co
rar.man.taxmin5.minc85.ma.05.greater.ma.df <- subset(rare.data.taxmin5.minC85.6000.phy.soil.s0.b.dds.ran.05.df, log2FoldChange > 0)

#phyloseq taxa table into dataframe
rare.data.taxmin5.minC85.6000.phy.soil.s0.b.tax <- data.frame(tax_table(rare.data.taxmin5.minC85.6000.phy.soil.s0.b))

#merge taxa table with significantly more abundant co OTUs to tell what otus actually are
rar.man.taxmin5.minc85.05.greater.co.df..merged.df <- merge(rar.man.taxmin5.minc85.r1.r6.ma.05.greater.co.df, rare.data.taxmin5.minC85.6000.phy.soil.s0.b.tax, by = "row.names")

#merge taxa table with significantly more abundant ma OTUs to tell what otus actually are
rar.man.taxmin5.minc85.05.greater.ma.df..merged.df <- merge(rar.man.taxmin5.minc85.r1.r6.ma.05.greater.ma.df, rare.data.taxmin5.minC85.6000.phy.soil.s0.b.tax, by = "row.names")



#adding column for those sig greater in co than ma
rare.data.taxmin5.minC85.6000.phy.soil.s0.psmelt$co_more_ma <- ifelse(rare.data.taxmin5.minC85.6000.phy.soil.s0.psmelt$OTU %in% unique(rar.man.taxmin5.minc85.05.greater.co.df..merged.df$Row.names), "co_more_ma", "none")
head(rare.data.taxmin5.minC85.6000.phy.soil.s0.psmelt)

#adding column for those sig greater in ma than rain co
rare.data.taxmin5.minC85.6000.phy.soil.s0.psmelt$ma_more_co <- ifelse(rare.data.taxmin5.minC85.6000.phy.soil.s0.psmelt$OTU %in% unique(rar.man.taxmin5.minc85.05.greater.ma.df..merged.df$Row.names), "ma_more_co", "none")
head(rare.data.taxmin5.minC85.6000.phy.soil.s0.psmelt)

saveRDS(rare.data.taxmin5.minC85.6000.phy.soil.s0.psmelt, "rare.data.taxmin5.minC85.6000.phy.soil.s0.psmelt.RDS")
rare.data.taxmin5.minC85.6000.phy.soil.s0.b.psmelt<-readRDS("rare.data.taxmin5.minC85.6000.phy.soil.s0.psmelt.RDS")


rare.data.taxmin5.minC85.6000.phy.soil.s0.b.psmelt<-subset(rare.data.taxmin5.minC85.6000.phy.soil.s0.b.psmelt, matrix=="soil")
rare.data.taxmin5.minC85.6000.phy.soil.s0.b.psmelt.s<-subset(rare.data.taxmin5.minC85.6000.phy.soil.s0.b.psmelt, soil_event=="S0")
rare.data.taxmin5.minC85.6000.phy.soil.s0.b.final.psmelt<-subset(rare.data.taxmin5.minC85.6000.phy.soil.s0.b.psmelt, soil_layer=="B")



saveRDS(rare.data.taxmin5.minC85.6000.phy.soil.s0.b.final.psmelt, "rare.data.taxmin5.minC85.6000.phy.soil.s0.b.final.psmelt.RDS")

soil.s0.b.final.dply <- ddply(rare.data.taxmin5.minC85.6000.phy.soil.s0.b.final.psmelt, .(id, manure_history, soil_layer, phylum, class, order, family, genus, ma_more_co, co_more_ma), summarise, total = sum(Abundance))




soil.s0.b.final.dply.genus <- ddply(rare.data.taxmin5.minC85.6000.phy.soil.s0.b.final.psmelt, .(genus, ma_more_co, co_more_ma), summarise, total = sum(Abundance))


soil.s0.b.final.dply.genus.more.ma<-subset(soil.s0.b.final.dply.genus, ma_more_co=="ma_more_co")
soil.s0.b.final.dply.genus.more.co<-subset(soil.s0.b.final.dply.genus, co_more_ma=="co_more_ma")


################################################################################################
#################finding which OTUs are significantly greater in manured S1 than S0 "describing short term manure affect"

rare.data.taxmin5.minC85.6000.short.term.phy<-readRDS("rare.data.taxmin5.minC85.6000.phy.RDS")

rare.data.taxmin5.minC85.6000.short.term.phy.soil<-subset_samples(rare.data.taxmin5.minC85.6000.short.term.phy, matrix== "soil")

rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured<-subset_samples(rare.data.taxmin5.minC85.6000.short.term.phy.soil, manure_history== "Ma")

rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0<-subset_samples(rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured, soil_event != "S3")

rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1<-subset_samples(rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0, soil_event != "S2")

rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top<-subset_samples(rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1, soil_layer== "T")

#converting physeq object to dseq object
rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top.dds <- phyloseq_to_deseq2(rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top, ~ soil_event)

#run deseq
rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top.dds.ran <- DESeq(rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top.dds, test="Wald", fitType="parametric")


#turn deseq object into dataframe
rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top.dds.ran.df <- data.frame(results(rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top.dds.ran, cooksCutoff=F))
head(rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top.dds.ran.df)
levels(rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top.dds.ran$soil_event)

remove all OTUs w/ pvalue greater than 0.05
rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top.dds.ran.df.ran.05.df <- subset(rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top.dds.ran.df, pvalue < 0.05)
head(rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top.dds.ran.df.ran.05.df)


#keep OTUs sig more abundant in s1 than s0
rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top.dds.ran.df.ran.05.greater.s1.df <- subset(rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top.dds.ran.df.ran.05.df, log2FoldChange > 0)

#phyloseq taxa table into dataframe
rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top.tax <- data.frame(tax_table(rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top))

#merge taxa table with significantly more abundant S1 OTUs to tell what otus actually are
rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top.dds.ran.df.ran.05.greater.s1.df..merged.df <- merge(rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top.dds.ran.df.ran.05.greater.s1.df, rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top.tax, by = "row.names")

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0 <- readRDS("rare.data.taxmin5.minC85.6000.psmelt")

#adding column for those sig greater in co than ma
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0$s1_more_s0 <- ifelse(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0$OTU %in% unique(rare.data.taxmin5.minC85.6000.short.term.phy.soil.manured.s0.s1.top.dds.ran.df.ran.05.greater.s1.df..merged.df$Row.names), "s1_more_s0", "none")
head(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0)


saveRDS(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0, "rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.RDS")

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0<-readRDS("rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.RDS")

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0, matrix=="soil")
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil, soil_layer=="T")
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top, s1_more_s0=="s1_more_s0")
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1, manure_history=="Ma")


###########manure dply
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.manure<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0, matrix=="manure")
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.manure.more.s1<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.manure, s1_more_s0=="s1_more_s0")


rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.manure.more.s1 <- ddply(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.manure.more.s1, .(matrix, phylum, class, order, family, genus, s1_more_s0), summarise, total = sum(Abundance))



#######genus level dply

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist.dply <- ddply(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist, .(manure_history, soil_event, soil_layer, phylum, class, order, family, genus, s1_more_s0), summarise, total = sum(Abundance))

write.table(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist.dply, "./rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist.dply.txt", sep="\t", quote=F, row.names=F)

short.term.soil.rel.abund<-read.delim("rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist.dply.w.relabund.txt", header = TRUE, sep = "\t")


ggplot(short.term.soil.rel.abund, aes(x=phylum, y=relative_abundance, fill = phylum)) +  ylim(-1,0.25) +
	geom_bar(stat="identity")+facet_grid(~soil_event)+coord_polar(start = 0)

#####phylum level dpy
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist.dply.phylum <- ddply(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist, .(manure_history, soil_event, soil_layer, phylum, s1_more_s0), summarise, total = sum(Abundance))

write.table(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist.dply.phylum, "./rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist.dply.phylum.txt", sep="\t", quote=F, row.names=F)

short.term.soil.rel.abund.phyla<-read.delim("rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist.dply.phylum.relabund.txt", header = TRUE, sep = "\t")

ggplot(short.term.soil.rel.abund.phyla, aes(x=Phylum, y=rel_percent, fill = Phylum)) +  ylim(-25,50) + 
	geom_bar(stat="identity")+facet_grid(~soil_event)+coord_polar(start = 0)
	
########class level dply
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist.dply.class <- ddply(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist, .(manure_history, soil_event, soil_layer, phylum, class, s1_more_s0), summarise, total = sum(Abundance))

write.table(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist.dply.class, "./rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist.dply.class.txt", sep="\t", quote=F, row.names=F)

#############
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist.dply.order <- ddply(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist, .(manure_history, soil_event, soil_layer, phylum, class, order, s1_more_s0), summarise, total = sum(Abundance))

write.table(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist.dply.order, "./rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist.dply.order.txt", sep="\t", quote=F, row.names=F)

short.term.soil.rel.abund.order<-read.delim("rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist.dply.order.rel.abund.txt", header = TRUE, sep = "\t")



short.term.soil.rel.abund.order.dply <- ddply(short.term.soil.rel.abund.order, .(soil_event, Phylum, Order, s1_more_s0), summarise, percent_relative = sum(rel_percent))

short.term.soil.rel.abund.order.dply<-short.term.soil.rel.abund.order.dply[with(short.term.soil.rel.abund.order.dply, order(Phylum)), ]

short.term.soil.rel.abund.order.dply$Order = factor(short.term.soil.rel.abund.order.dply$Order, levels=c('Actinomycetales','Bacteroidales','Flavobacteriales','Sphingobacteriales', 'Clostridiales', 'Burkholderiales', 'Caulobacterales', 'Pseudomonadales', 'Rhizobiales', 'Sphingomonadales', 'Xanthomonadales', 'Spirochaetales', 'Other'))

ggplot(short.term.soil.rel.abund.order.dply, aes(x=Order, y=percent_relative, fill = Phylum)) +  ylim(-5,22) + 
	geom_bar(stat="identity")+theme_bw() +theme(axis.text.x=element_blank())+facet_grid(~soil_event) +facet_wrap( ~ soil_event, ncol=2)+coord_polar(start = 0)
	
ggplot(short.term.soil.rel.abund.order.dply, aes(x=Order, y=percent_relative, fill = Phylum)) + 
	geom_bar(stat="identity")+theme_bw() +facet_grid(~soil_event) +facet_wrap( ~ soil_event, ncol=2)+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
	
	
###########################boxplot for order level soil impacted organisms
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist.dply.order.samples <- ddply(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist, .(id, manure_history, soil_event, soil_layer, phylum, class, order, s1_more_s0), summarise, total = sum(Abundance))

write.table(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist.dply.order.samples, "./rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist.dply.order.samples.txt", sep="\t", quote=F, row.names=F)

man.affect.sub<-rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.top.only.more.s1.manhist.dply.order.samples

man.affect.sub.done<-subset(man.affect.sub, man.affect.sub$order== "Actinomycetales" | man.affect.sub$order=="Bacteroidales" | man.affect.sub$order=="Flavobacteriales" | man.affect.sub$order=="Sphingobacteriales" | man.affect.sub$order=="Clostridiales" | man.affect.sub$order== "Burkholderiales" | man.affect.sub$order== "Caulobacterales" | man.affect.sub$order== "Pseudomonadales" | man.affect.sub$order=="Rhizobiales" | man.affect.sub$order== "Sphingomonadales" | man.affect.sub$order=="Xanthomonadales" | man.affect.sub$order== "Spirochaetales")

man.affect.sub.done$rel_abund<- (man.affect.sub.done$total)/6076

man.affect.sub.done <- summarySE(man.affect.sub.done, measurevar="rel_abund", groupvars=c("soil_event", "order"))

man.affect.sub.done$soil_event <- gsub("S0", "Day 0", man.affect.sub.done$soil_event)
man.affect.sub.done$soil_event <- gsub("S1", "Day 24", man.affect.sub.done$soil_event)
man.affect.sub.done$soil_event <- gsub("S2", "Day 59", man.affect.sub.done$soil_event)
man.affect.sub.done$soil_event <- gsub("S3", "Day 108", man.affect.sub.done$soil_event)
man.affect.sub.done$soil_event = factor(man.affect.sub.done$soil_event, levels=c('Day 0', 'Day 24','Day 59','Day 108'))

ggplot(man.affect.sub.done, aes(x=soil_event, y=rel_abund, group=order, color= order)) + 
    geom_errorbar(aes(ymin= rel_abund-se, ymax= rel_abund +se), size=0.5, width=.3) +
    geom_line(aes(linetype= order, size=0.1)) + 
    geom_point(aes(shape= order, size=0.5))+ facet_grid(~order)+facet_wrap( ~ order, ncol=3)+theme_bw(15) + xlab("Days After Manure Application") +ylab("Relative Abundance")


#################################################################################    
######################short term increases in water##############################

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0, matrix=="water")
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.more.s1<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water, s1_more_s0=="s1_more_s0")
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.more.s1.ma<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.more.s1, manure_history=="Ma")


###################class level dply
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.more.s1.ma.class <- ddply(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.more.s1.ma, .(manure_history, rain_event, phylum, class, s1_more_s0), summarise, total = sum(Abundance))

write.table(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.more.s1.ma.class, "./rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.more.s1.ma.class.txt", sep="\t", quote=F, row.names=F)

###################order level dply

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.more.s1.ma.order <- ddply(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.more.s1.ma, .(manure_history, rain_event, phylum, class, order, s1_more_s0), summarise, total = sum(Abundance))

write.table(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.more.s1.ma.order, "./rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.more.s1.ma.order.txt", sep="\t", quote=F, row.names=F)

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.more.s1.ma.order.samples <- ddply(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.more.s1.ma, .(id, manure_history, rain_event, phylum, class, order, s1_more_s0), summarise, total = sum(Abundance))

man.affect.water.sub<-rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.more.s1.ma.order.samples 

man.affect.water.sub.done<-subset(man.affect.water.sub, man.affect.water.sub$order== "Actinomycetales" | man.affect.water.sub$order=="Bacteroidales" | man.affect.water.sub$order=="Flavobacteriales" | man.affect.water.sub$order=="Sphingobacteriales" | man.affect.water.sub$order=="Clostridiales" | man.affect.water.sub$order== "Burkholderiales" | man.affect.water.sub$order== "Caulobacterales" | man.affect.water.sub$order== "Pseudomonadales" | man.affect.water.sub$order=="Rhizobiales" | man.affect.water.sub$order== "Sphingomonadales" | man.affect.water.sub$order=="Xanthomonadales" |man.affect.water.sub$order== "Spirochaetales")

man.affect.water.sub.done.r2.r4.r6<-subset(man.affect.water.sub.done, man.affect.water.sub.done$rain_event=="R2" | man.affect.water.sub.done$rain_event=="R4" | man.affect.water.sub.done$rain_event=="R6")


man.affect.water.sub.done.r2.r4.r6$rel_abund<- (man.affect.water.sub.done.r2.r4.r6$total)/6076

man.affect.water.sub.done.r2.r4.r6.se <- summarySE(man.affect.water.sub.done.r2.r4.r6, measurevar="rel_abund", groupvars=c("rain_event", "order"))


man.affect.water.sub.done.r2.r4.r6.se$rain_event <- gsub("R2", "Day 24", man.affect.water.sub.done.r2.r4.r6.se$rain_event)
man.affect.water.sub.done.r2.r4.r6.se$rain_event <- gsub("R4", "Day 59", man.affect.water.sub.done.r2.r4.r6.se$rain_event)
man.affect.water.sub.done.r2.r4.r6.se$rain_event <- gsub("R6", "Day 108", man.affect.water.sub.done.r2.r4.r6.se$rain_event)
man.affect.water.sub.done.r2.r4.r6.se$rain_event = factor(man.affect.water.sub.done.r2.r4.r6.se$rain_event, levels=c('Day 24','Day 59','Day 108'))

ggplot(man.affect.water.sub.done.r2.r4.r6.se, aes(x=rain_event, y=rel_abund, group=order, color= order)) + 
    geom_errorbar(aes(ymin= rel_abund-se, ymax= rel_abund +se), size=0.5, width=.3) +
    geom_line(aes(linetype= order, size=0.1)) + 
    geom_point(aes(shape= order, size=0.5))+ facet_grid(~order)+facet_wrap( ~ order, ncol=3)+theme_bw(15) + xlab("Days After Manure Application") +ylab("Relative Abundance")

###################################MSO figure 5
####################combining soil and water s0-s1 charts

write.table(man.affect.sub.done, "./man.affect.sub.done.txt", sep="\t", quote=F, row.names=F)

write.table(man.affect.water.sub.done.r2.r4.r6.se, "./man.affect.water.sub.done.r2.r4.r6.se.txt", sep="\t", quote=F, row.names=F)

man.affect.soil.water.se
man.affect.soil.water.se<-read.delim("man.affect.soil.water.se.txt", header = TRUE, sep = "\t")

    
man.affect.soil.water.se$soil_event = factor(man.affect.soil.water.se$soil_event, levels=c('Day 0', 'Day 24','Day 59','Day 108'))
    
ggplot(man.affect.soil.water.se, aes(x=soil_event, y=rel_abund, group=sample_type, color=sample_type))+geom_line(size=1)+facet_grid(~order)+facet_wrap( ~ order, ncol=3)+ 
    geom_errorbar(aes(ymin= rel_abund-se, ymax= rel_abund +se, width=.5))+geom_point()+ xlab("Days After Manure Application") +ylab("Relative Abundance")+theme_bw(15)

ggsave("manure_affected_orders_soil_water.png", width =11)


##################### water short term affected pie charts

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0, matrix!="soil")
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.ma<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water, manure_history=="Ma")

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.ma.r2.r4.r6<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.ma, rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.ma$rain_event != "R1" | rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.ma$rain_event!="R3"| rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.ma$rain_event !="R5")

s1.s0.water.pie.dply <- ddply(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0, .(matrix, rain_event, s1_more_s0, manure_history), summarise, total = sum(Abundance))



write.table(s1.s0.water.pie.dply, "./s1.s0.water.pie.dply.txt", sep="\t", quote=F, row.names=F)

s1.s0.water.pie.dply.w.rel.abund
s1.s0.water.pie.dply.w.rel.abund<-read.delim("s1.s0.water.pie.dply.w.rel.abund.2.txt", header = TRUE, sep = "\t")
s1.s0.water.pie.dply.w.rel.abund$rain_event = factor(s1.s0.water.pie.dply.w.rel.abund$rain_event, levels=c('Day 0','Day 24','Day 59', 'Day 108'))



ggplot(s1.s0.water.pie.dply.w.rel.abund, aes(x=matrix, y=Relative_Abundance, fill = s1_more_s0))+ geom_bar(stat="identity") + facet_grid(~rain_event, scale="free_x", space="free_x")+theme_bw(20)+xlab("Matrix")+ylab("Relative Abundance")+theme(axis.text.x=element_text(angle = 45, vjust = 0.5))



ggsave("manure_affected_water_bar.png", width =14)

########################soil short term impacted


saveRDS(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0, "rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.RDS")
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0<-readRDS("rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.RDS")

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0, matrix=="soil")
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.manured<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil, manure_history=="Ma")
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.manured.top<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.manured, soil_layer=="T")

s1.s0.soil.affected.dply <- ddply(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0, .(matrix, soil_event, s1_more_s0, manure_history, soil_layer), summarise, total = sum(Abundance))


write.table(s1.s0.soil.affected.dply, "./s1.s0.soil.affected.dply.txt", sep="\t", quote=F, row.names=F)
s1.s0.soil.affected.dply.w.rel.abund<-read.delim("s1.s0.soil.affected.dply.w.rel.abund.txt", header = TRUE, sep = "\t")

s1.s0.soil.affected.dply.w.rel.abund$soil_event = factor(s1.s0.soil.affected.dply.w.rel.abund$soil_event, levels=c('Pre-Manure Application','Day 24','Day 59', 'Day 108'))


ggplot(s1.s0.soil.affected.dply.w.rel.abund, aes(x=soil_event, y=Relative_Abundance, fill = s1_more_s0))+ geom_bar(stat="identity") + facet_grid(~matrix, scale="free_x", space="free_x")+theme_bw(20)+xlab("")+ylab("Relative Abundance")+theme(axis.text.x=element_text(angle = 45, hjust = 1)) + scale_fill_manual(values=c("#003300", "#FF9933"))

003300 FF0000
ggsave("s1.s0.soil.affected.dply.w.rel.abund.png", width =18)




###############manure core


rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0<-readRDS("rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.RDS")



rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.mcore<- subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0, rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0$OTU == "OTU_1132" | rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0$OTU == "OTU_26418" | rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0$OTU == "OTU_26484" | rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0$OTU == "OTU_26487" | rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0$OTU == "OTU_26490" | rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0$OTU == "OTU_3835" | rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0$OTU == "OTU_5213")

mcore.dply <- ddply(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.mcore, .(matrix, sample_id, OTU, genus, s1_more_s0, soil_event, rain_event, soil_layer, manure_history), summarise, total = sum(Abundance))

write.table(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.mcore, "./rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.mcore.txt", sep="\t", quote=F, row.names=F)

mcore.dply.manure<-subset(mcore.dply, matrix=="soil")
mcore.dply.manure<-subset(mcore.dply.manure, soil_event=="S0")
mcore.dply.manure<-subset(mcore.dply.manure, soil_layer=="T")
mcore.dply.manure<-subset(mcore.dply.manure, manure_history=="Co")

write.table(mcore.dply.manure, "./mcore.dply.manure.txt", sep="\t", quote=F, row.names=F)


mcore.dply.r2.r4.r6<-subset(mcore.dply, mcore.dply$rain_event=="R2" | mcore.dply$rain_event=="R4" | mcore.dply$rain_event=="R6")
mcore.dply.r2.r4.r6<-subset(mcore.dply.r2.r4.r6, manure_history=="Ma")

mcore.dply.r2.r4.r6$rel_abund<- (mcore.dply.r2.r4.r6$total)/6076

mcore.dply.r2.r4.r6.se <- summarySE(mcore.dply.r2.r4.r6, measurevar="rel_abund", groupvars=c("rain_event", "genus"))


man.affect.water.sub.done.r2.r4.r6.se$rain_event <- gsub("R2", "Day 24", man.affect.water.sub.done.r2.r4.r6.se$rain_event)
man.affect.water.sub.done.r2.r4.r6.se$rain_event <- gsub("R4", "Day 59", man.affect.water.sub.done.r2.r4.r6.se$rain_event)
man.affect.water.sub.done.r2.r4.r6.se$rain_event <- gsub("R6", "Day 108", man.affect.water.sub.done.r2.r4.r6.se$rain_event)
man.affect.water.sub.done.r2.r4.r6.se$rain_event = factor(man.affect.water.sub.done.r2.r4.r6.se$rain_event, levels=c('Day 24','Day 59','Day 108'))



ggplot(mcore.dply.r2.r4.r6.se, aes(x=rain_event, y=rel_abund, group=genus, color=genus))+geom_line(size=1)+facet_grid(~genus)+facet_wrap( ~ genus, ncol=3)+ 
    geom_errorbar(aes(ymin= rel_abund-se, ymax= rel_abund +se, width=.5))+geom_point()+ xlab("Days After Manure Application") +ylab("Relative Abundance")+theme_bw(15)
    
    
mcore.dply.soil<-subset(mcore.dply, matrix=="soil")

mcore.dply.soil.man<-subset(mcore.dply.soil, manure_history=="Ma")

mcore.dply.soil.man.top<-subset(mcore.dply.soil.man, soil_layer=="T")

mcore.dply.soil.man.top$rel_abund<- (mcore.dply.soil.man.top$total)/6076

mcore.dply.soil.man.top.se <- summarySE(mcore.dply.soil.man.top, measurevar="rel_abund", groupvars=c("soil_event", "genus"))

ggplot(mcore.dply.soil.man.top.se, aes(x=soil_event, y=rel_abund, group=genus, color=genus))+geom_line(size=1)+facet_grid(~genus)+facet_wrap( ~ genus, ncol=3)+ 
    geom_errorbar(aes(ymin= rel_abund-se, ymax= rel_abund +se, width=.5))+geom_point()+ xlab("Days After Manure Application") +ylab("Relative Abundance")+theme_bw(15)
    

write.table(mcore.dply.r2.r4.r6.se, "./mcore.dply.r2.r4.r6.se.txt", sep="\t", quote=F, row.names=F)
write.table(mcore.dply.soil.man.top.se, "./mcore.dply.soil.man.top.se.txt", sep="\t", quote=F, row.names=F)


mcore.se.combined<-read.delim("mcore.se.combined.txt", header = TRUE, sep = "\t")

mcore.se.combined$soil_event = factor(mcore.se.combined$soil_event, levels=c('Day 0', 'Day 24','Day 59','Day 108'))
    
cbPalette <- c("#003300", "#FF9933")
, "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


ggplot(mcore.se.combined, aes(x=soil_event, y=rel_abund, group=sample_type, color=sample_type))+geom_line(aes(group=1), colour= cbPalette)+facet_grid(~genus)+facet_wrap( ~ genus, ncol=3)+ 
    geom_errorbar(aes(ymin= rel_abund-se, ymax= rel_abund +se, width=.5))+geom_point()+ xlab("Days After Manure Application") +ylab("Relative Abundance")+ scale_fill_manual(values=c("#003300", "#FF9933"))+theme_bw(15)
 
ggsave("manured.soil.water.phyla.png", width =13)
003300 FF0000

######################################################################
##############making bar chart of MSOs####################
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0<-readRDS("rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.RDS")

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices <- ddply(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0, .(manure_history, soil_layer, experiment_day, rain_event, phylum, matrix, soil_layer, s1_more_s0), summarise, total = sum(Abundance))

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.water<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices, matrix=="water")

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.water.manured<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.water, manure_history=="Ma")

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.water.manured.mores1<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.water.manured, s1_more_s0=="s1_more_s0")




write.table(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.water.manured.mores1, "./rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.water.manured.mores1.txt", sep="\t", quote=F, row.names=F)


rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.water.manured.mores1.edited<-read.delim("rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.water.manured.mores1.edited.txt", header = TRUE, sep = "\t")

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.water.manured.mores1.edited$day = factor(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.water.manured.mores1.edited$day, levels=c('Day 0', 'Day 10','Day 24','Day 38','Day 59', 'Day 80', 'Day 108'))



ggplot(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.water.manured.mores1.edited, aes(x=day, y=relative_abundance, fill = Phylum))+ geom_bar(stat="identity")+ facet_grid(. ~ matrix.1, scale="free_x", space="free_x")+theme_bw(20)+xlab("Experiment Day")+ylab("Relative Abundance")+theme(axis.text.x=element_text(angle = 45, vjust = 0.5))
 
 
 
 ##soil addition
 
 rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.soil<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices, matrix=="soil")

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.soil.manured<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.soil, manure_history=="Ma")

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.soil.manured.mores1<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.soil.manured, s1_more_s0=="s1_more_s0")

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.soil.manured.mores1.t<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.soil.manured.mores1, soil_layer=="T")

write.table(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.soil.manured.mores1.t, "./rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.soil.manured.mores1.t.txt", sep="\t", quote=F, row.names=F)

############manure additions

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.manure<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices, matrix=="manure")
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.manure.mores1<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.manure, s1_more_s0=="s1_more_s0")

write.table(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.manure.mores1, "./rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.all.matrices.manure.txt", sep="\t", quote=F, row.names=F)

#############overall manure, man history soil and drainage nmds
rare.data.taxmin5.minC85.6000.phy.man = subset_samples(rare.data.taxmin5.minC85.6000.phy, manure_history != "Co")
rare.data.taxmin5.minC85.6000.phy.man.top.soil = subset_samples(rare.data.taxmin5.minC85.6000.phy.man, soil_layer != "B")


sample_data(rare.data.taxmin5.minC85.6000.phy.man.top.soil)$experiment_day<-as.character(sample_data(rare.data.taxmin5.minC85.6000.phy.man.top.soil)$experiment_day)

rare.data.taxmin5.minC85.6000.phy.man.top.soil$experiment_day 

sample_data(rare.data.taxmin5.minC85.6000.phy.man.top.soil)$experiment_day = factor(sample_data(rare.data.taxmin5.minC85.6000.phy.man.top.soil)$experiment_day, levels = c("0","10","24","38", "59", "80", "108"))


rare.data.taxmin5.minC85.6000.phy.man.ord <- ordinate(rare.data.taxmin5.minC85.6000.phy.man.top.soil, "NMDS", "bray")

manure_sample_plot = plot_ordination(rare.data.taxmin5.minC85.6000.phy.man.top.soil, rare.data.taxmin5.minC85.6000.phy.man.ord, type="samples", shape = "matrix", color="experiment_day")

manure_sample_plot + geom_point(size=5)  + stat_ellipse(type = "norm", linetype = 2,)+theme_bw(20)+ theme(legend.title = element_blank()) 



ggsave("nmds.manured.top.soil.water.manure.png", width =7)






##################################manure affected nmds

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0<-readRDS("rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.RDS")
justs1.more.s2.otus<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0, s1_more_s0=="s1_more_s0")
justs1.more.s2.otus.dply <- ddply(justs1.more.s2.otus, .(OTU), summarise, total = sum(Abundance))

rare.data.taxmin5.minC85.6000.phy.man.top.soil.subset <- subset_taxa(rare.data.taxmin5.minC85.6000.phy.man.top.soil, rownames(tax_table(rare.data.taxmin5.minC85.6000.phy.man.top.soil)) %in% justs1.more.s2.otus.dply$OTU)


rare.data.taxmin5.minC85.6000.phy.man.top.soil.subset


sample_data(rare.data.taxmin5.minC85.6000.phy.man.top.soil.subset)$experiment_day = factor(sample_data(rare.data.taxmin5.minC85.6000.phy.man.top.soil.subset)$experiment_day, levels = c("0","10","24","38", "59", "80", "108"))


rare.data.taxmin5.minC85.6000.phy.man.top.soil.subset.ord <- ordinate(rare.data.taxmin5.minC85.6000.phy.man.top.soil.subset, "NMDS", "bray")

manure_sample_plot_subset = plot_ordination(rare.data.taxmin5.minC85.6000.phy.man.top.soil.subset, rare.data.taxmin5.minC85.6000.phy.man.top.soil.subset.ord, type="samples", shape = "matrix", color= "experiment_day")

manure_sample_plot_subset + geom_point(size=5)  + stat_ellipse(type = "norm", linetype = 2)+theme_bw(20)+ theme(legend.title = element_blank())


ggsave("nmds.manure.affected.top.soil.water.manure.png", width =7)


rare.data.taxmin5.minC85.6000.phy<-readRDS("rare.data.taxmin5.minC85.6000.phy.RDS")
rare.data.taxmin5.minC85.6000.soil.top<-subset_samples(rare.data.taxmin5.minC85.6000.phy, soil_layer !="B")
rare.data.taxmin5.minC85.6000.soil.top.melt<-psmelt(rare.data.taxmin5.minC85.6000.soil.top)



rare.data.taxmin5.minC85.6000.soil.top.melt.dply <- ddply(rare.data.taxmin5.minC85.6000.soil.top.melt, .(id, matrix, soil_event, rain_event, experiment_day, soil_layer, manure_history), summarise, total = sum(Abundance))
write.table(rare.data.taxmin5.minC85.6000.soil.top.melt.dply, "./rare.data.taxmin5.minC85.6000.samples.txt", sep="\t", quote=F, row.names=F)


rare.data.taxmin5.minC85.6000.phy<-readRDS("rare.data.taxmin5.minC85.6000.phy.RDS")

soil.column.otu.table<-otu_table(rare.data.taxmin5.minC85.6000.phy)
soil.column.tax.table<-tax_table(rare.data.taxmin5.minC85.6000.phy)
soil.column.sample.data<-sample_data(rare.data.taxmin5.minC85.6000.phy)

write.table(soil.column.otu.table, "./soil.column.otu.table.txt", sep="\t", quote=F, row.names=F)
write.table(soil.column.tax.table, "./soil.column.tax.table.txt", sep="\t", quote=F, row.names=F)
write.table(soil.column.sample.data, "./soil.column.sample.data.txt", sep="\t", quote=F, row.names=F)

######################################################
#######Frontiers shifts in non-manured control columns

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0<-readRDS("rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.RDS")
head(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0)

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0, matrix=="soil" & soil_layer=="T" & manure_history=="Co")


rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.dply <- ddply(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil, .(soil_event, phylum), summarise, total = sum(Abundance))

write.table(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.dply, "./test.incubation.effect.mso.txt", sep="\t", quote=F, row.names=F)

###########################################
#######Average relative abundances and std of major phyla in soil

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0<-readRDS("rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.RDS")

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.manure<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0, matrix=="soil" & soil_layer=="T" & manure_history=="Ma")


rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.manure.major.phyla<- subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.manure, phylum=="Acidobacteria" | phylum=="Actinobacteria" | phylum=="Bacteroidetes" | phylum=="candidate division WPS-1" | phylum=="Chloroflexi" |phylum=="Cyanobacteria/Chloroplast" | phylum=="Euryarchaeota" | phylum=="Firmicutes" | phylum=="Other" | phylum=="Planctomycetes" | phylum=="Proteobacteria" | phylum=="Spirochaetes" | phylum=="unclassified_Bacteria" | phylum=="Verrucomicrobia")


rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.manure.major.phyla.dply <- ddply(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.manure.major.phyla, .(id, soil_event, phylum), summarise, total = sum(Abundance))


rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.manure.major.phyla.dply$rel_abund<-(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.manure.major.phyla.dply$total)/6076

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.manure.major.phyla.dply.se <- summarySE(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.manure.major.phyla.dply, measurevar="rel_abund", groupvars=c("soil_event", "phylum"))
 
write.table(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.manure.major.phyla.dply.se, "./soil.major.phyla.avg.stddev.txt", sep="\t", quote=F, row.names=F)


###########################################
#######Average relative abundances and std of major phyla in water

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0<-readRDS("rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.RDS")

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.manured<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0, matrix=="water" & manure_history=="Ma")


rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.manure.major.phyla<- subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.manured, phylum=="Acidobacteria" | phylum=="Actinobacteria" | phylum=="Bacteroidetes" | phylum=="candidate division WPS-1" | phylum=="Chloroflexi" |phylum=="Cyanobacteria/Chloroplast" | phylum=="Euryarchaeota" | phylum=="Firmicutes" | phylum=="Other" | phylum=="Planctomycetes" | phylum=="Proteobacteria" | phylum=="Spirochaetes" | phylum=="unclassified_Bacteria" | phylum=="Verrucomicrobia")


rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.manure.major.phyla.dply <- ddply(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.manure.major.phyla, .(id, rain_event, phylum), summarise, total = sum(Abundance))

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.manure.major.phyla.dply$rel_abund<-(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.manure.major.phyla.dply$total)/6076


rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.manure.major.phyla.dply.se <- summarySE(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.manure.major.phyla.dply, measurevar="rel_abund", groupvars=c("rain_event", "phylum"))
 
write.table(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.manure.major.phyla.dply.se, "./water.major.phyla.avg.stddev.txt", sep="\t", quote=F, row.names=F)

#################################################
############significant differences in major phyla between r1 and r6

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.manure.major.phyla.dplyr1.r6<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.manure.major.phyla.dply, rain_event=="R1" | rain_event=="R6")

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.manure.major.phyla.dplyr1.r6.single.phylum<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.manure.major.phyla.dplyr1.r6, phylum=="Verrucomicrobia")

wilcox.test(rel_abund ~ rain_event, data = rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.water.manure.major.phyla.dplyr1.r6.single.phylum)

#################################################
##################significant differences between s0 and s3

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.manure.major.phyla.dply.s0.s3<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.manure.major.phyla.dply, soil_event=="S0" 
| soil_event=="S3")

rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.manure.major.phyla.dply.s0.s3.single.phylum<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.manure.major.phyla.dply.s0.s3, phylum=="Verrucomicrobia")

wilcox.test(rel_abund ~ soil_event, data = rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.soil.manure.major.phyla.dply.s0.s3.single.phylum, exact = FALSE)

phylum=="Acidobacteria" | phylum=="Actinobacteria" | phylum=="Bacteroidetes" | phylum=="candidate division WPS-1" | phylum=="Chloroflexi" |phylum=="Cyanobacteria/Chloroplast" | phylum=="Euryarchaeota" | phylum=="Firmicutes" | phylum=="Other" | phylum=="Planctomycetes" | phylum=="Proteobacteria" | phylum=="Spirochaetes" | phylum=="unclassified_Bacteria" | phylum=="Verrucomicrobia")


###################
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0<-readRDS("rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.RDS")
justs1.more.s2.otus<-subset(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0, s1_more_s0=="s1_more_s0")
justs1.more.s2.otus.dply <- ddply(justs1.more.s2.otus, .(OTU, genus), summarise, total = sum(Abundance))

write.table(justs1.more.s2.otus.dply, "./mso.table.list.txt", sep="\t", quote=F, row.names=F)

#####################list of samples S_XX and id
rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0<-readRDS("rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0.RDS")
SXX_and_id_list <- ddply(rare.data.taxmin5.minC85.6000.psmelt.s1.greater.s0, .(Sample, id, soil), summarise, total = sum(Abundance))

write.table(SXX_and_id_list, "./SXX_and_id_list.txt", sep="\t", quote=F, row.names=F)

##############################################comparing adonis soil and drainage day 24 and day 108


saveRDS(rare.data.taxmin5.minC85.6000.phy.man.top.soil.subset, "rare.data.taxmin5.minC85.6000.phy.man.top.soil.subset.RDS")
mso.adonis<-readRDS("rare.data.taxmin5.minC85.6000.phy.man.top.soil.subset.RDS")

mso.adonis.soil.water.phy = subset_samples(mso.adonis, matrix != "manure")

mso.adonis.soil.water.manured.phy = subset_samples(mso.adonis.soil.water.phy, manure_history == "Ma")

mso.adonis.top.soil.water.manured.phy = subset_samples(mso.adonis.soil.water.manured.phy, soil_layer != "B")

mso.adonis.top.soil.water.manured.phy.d24 = subset_samples(mso.adonis.top.soil.water.manured.phy, experiment_day == "24")

mso.adonis.top.soil.water.manured.phy.d108 = subset_samples(mso.adonis.top.soil.water.manured.phy, experiment_day == "108")


mso.total.data.dist.24 = phyloseq::distance(mso.adonis.top.soil.water.manured.phy.d24, "bray")
adonis(mso.total.data.dist.24 ~ matrix , perm=9999, as(sample_data(mso.adonis.top.soil.water.manured.phy.d24), "data.frame"))

mso.total.data.dist.108 = phyloseq::distance(mso.adonis.top.soil.water.manured.phy.d108, "bray")
adonis(total.data.dist.108 ~ matrix , perm=9999, as(sample_data(mso.adonis.top.soil.water.manured.phy.d108), "data.frame"))

