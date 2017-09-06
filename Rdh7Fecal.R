#source("https://bioconductor.org/biocLite.R")
library(phyloseq)
library(ggplot2)
library(scales)
library(grid)
library(vegan)
###Set directory to where your OTU table, Mappingfile, and Tree are
setwd("/Users/kjduncan/Desktop/R/Rdh7_Fecal/")
list.files()
otufile="otu_table_filtered_min3samp"
mapfile="mappingfile__Region_fecal__.txt"
trefile="97_otus_unannotated.tree"
file<-import_biom(otufile)
map<-import_qiime_sample_data(mapfile)
treefile<-read_tree(trefile)
Analysis<-merge_phyloseq(file,map,treefile)
Analysis
colnames(tax_table(Analysis))<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
ntaxa(Analysis)
nsamples(Analysis)
sample_names(Analysis)
taxa_names(Analysis)[1:10]
sample_variables(Analysis)
###Abundance bar graphs
theme_set(theme_bw())
rank_names(Analysis)
colnames(tax_table(Analysis))<-c(k="Kingdom",p="Phylum",c="Class",o="Order",f="Family",g="Genus",s="Species")
subset=subset_taxa(Analysis, Kingdom == "k__Bacteria")
subset=transform_sample_counts(subset,function(x)x/sum(x))
plot_bar(subset,fill = "Family")
##PCoA Unweighted Unifrac
GP1=transform_sample_counts(Analysis, function(x) 1E6 * x/sum(x))
orduW = ordinate(GP1, "PCoA","unifrac",weighted=TRUE)
pW = plot_ordination(GP1, orduW, color = "Genotype", title = "Weighted Unifrac") + geom_point(size=5)
pW + theme(plot.title = element_text(size=18))+ theme(legend.text=element_text(size=14))
orduU = ordinate(GP1, "PCoA","unifrac",weighted=FALSE)
pU = plot_ordination(GP1, orduU, color = "Genotype", title = "Unweighted Unifrac") + geom_point(size=5)
pU + theme(plot.title = element_text(size=18))+ theme(legend.text=element_text(size=14))
###Permanova Weighted
set.seed(1)
analysis_unifrac_weighted<-phyloseq::distance(Analysis,method = "unifrac", weighted=TRUE)
sampledf<- data.frame(sample_data(Analysis))
adonis(analysis_unifrac_weighted ~ Genotype, data = sampledf)
beta<-betadisper(analysis_unifrac_weighted, sampledf$Genotype)
permutest(beta)
###Permanova Unweighted
analysis_unifrac_unweighted<-phyloseq::distance(Analysis,method = "unifrac", weighted=FALSE)
sampledf<- data.frame(sample_data(Analysis))
adonis(analysis_unifrac_unweighted ~ Genotype, data = sampledf)
beta<-betadisper(analysis_unifrac_unweighted, sampledf$Genotype)
permutest(beta)
                  