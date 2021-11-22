
#--------------- other things from Misha

# ANDSD => SFS for demographic analysis

# make separate files listing bams for each population (without clones and replicates)
# assume we have two populations, pop0 and pop1, 20 individuals each, with corresponding bams listed in pop0.bams and pop1.bams

# generating list of filtered SNP sites for SFS production (note: no filters that distort allele frequency!):
# sb - strand bias filter; only use for 2bRAD, GBS or WGS (not for ddRAD or RADseq)
# hetbias - detects weird heterozygotes because they have unequal representation of alleles 
# set minInd to 80-90% of all your individuals (depending on the results from quality control step)

nano sfs.sh
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd #start job in submission directory
#$ -N sfs.sh # job name, anything you want
#$ -l h_rt=24:00:00
#$ -M thenicolakriefall@gmail.com
#$ -m be
GENOME_REF="Amil.fasta"
FILTERS="-uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 91 -hwe_pval 1e-5"
TODO="-doMajorMinor 1 -doSaf 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 11"
angsd -b bams_no7 -GL 1 -anc $GENOME_REF $FILTERS $TODO -P 1 -out all.freq

# extracting and indexing list of sites to make SFS from 
# filtering out sites where heterozygote counts comprise more than 50% of all counts (likely lumped paralogs)
zcat all.freq.snpStat.gz | awk '($3+$4+$5+$6)>0' | awk '($12+$13+$14+$15)/($3+$4+$5+$6)<0.5' | cut -f 1,2  > sites2do
angsd sites index sites2do

# estimating site frequency likelihoods for each population, also saving allele frequencies (for genome scan) 
nano sfs_pops.sh
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd #start job in submission directory
#$ -N sfs_pops.sh # job name, anything you want
#$ -l h_rt=24:00:00
#$ -M thenicolakriefall@gmail.com
#$ -m be
GENOME_REF="./Amil.fasta"
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF
angsd -sites sites2do -b bams_mnwi -GL 1 -P 1 $TODO -out pop0
angsd -sites sites2do -b bams_mnwo -GL 1 -P 1 $TODO -out pop1
angsd -sites sites2do -b bams_msei -GL 1 -P 1 $TODO -out pop2
angsd -sites sites2do -b bams_mseo -GL 1 -P 1 $TODO -out pop3
angsd -sites sites2do -b bams_to -GL 1 -P 1 $TODO -out pop4
angsd -sites sites2do -b bams_ti -GL 1 -P 1 $TODO -out pop5

#exit
qsub sfs_pops.sh

# generating per-population SFS
#realSFS pop0.saf.idx >pop0.sfs
#realSFS pop1.saf.idx >pop1.sfs

#=========== Fst: global, per site, and per gene

# writing down 2d-SFS priors
realSFS pop0.saf.idx pop1.saf.idx -P 24 > p01.sfs ; realSFS fst index pop0.saf.idx pop1.saf.idx -sfs p01.sfs -fstout p01 

# global Fst between populations
realSFS fst stats p01.fst.idx

# per-site Fst
realSFS fst print p01.fst.idx > p01.fst

# extracting gene regions out of genome annotations file (gff3)
cat mygenome.gff3 | awk ' $3=="gene"' | cut -f 1,4,5,10 >gene_regions.tab

# extending to plus-minus 2 kb around the gene
awk '{print $1"\t"$2-2000"\t"$3+2000"\t"$4}' gene_regions.tab '> genes.txt

# use fstPerGene.R to compute per-gene Fst

#####	prepare the fst for easy window analysis etc 	#########

realSFS pop0.saf.idx pop1.saf.idx > p01.ml
realSFS fst index pop0.saf.idx pop1.saf.idx -sfs p01.ml -fstout p01comp


#### calculating theta #### 

echo "realSFS STJ_freq_fold.saf.idx > stj.sfs" > sfs
ls5_launcher_creator.py -j sfs -n sfs -l sfs_job -t 1:00:00 -w 1 -A tagmap 

echo "realSFS STA_freq_fold.saf.idx > sta.sfs" > sfs
ls5_launcher_creator.py -j sfs -n sfs -l sfs_job -t 1:00:00 -w 1 -A tagmap 


## if your reference genome is not the same as your study species, then use -anc for the genome ref 

TODO="-doSaf 1 -doThetas 1 -fold 1" 
echo "angsd -b ST_J -out STJ_theta -pest stj.sfs -anc $GENOME_REF $TODO -GL 1"> theta
ls5_launcher_creator.py -j theta -n theta -l thetajob -t 1:00:00 -w 1 -A tagmap 

TODO="-doSaf 1 -doThetas 1 -fold 1" 
echo "angsd -b ST_A -out STA_theta -pest sta.sfs -anc $GENOME_REF $TODO -GL 1"> theta
ls5_launcher_creator.py -j theta -n theta -l thetajob -t 1:00:00 -w 1 -A tagmap 

 
#global theta estimate for every Chromosome/scaffold
#thetaStat is an angsd subprogram, so it should be installed already with angsd 

thetaStat do_stat STA_theta.thetas.idx

## sliding window analysis of theta with 50000 kb window with 10000 kb step size if interested in regions of high/low diversity

thetaStat do_stat STJ_theta.thetas.idx -win 50000 -step 10000  -outnames STJ.theta
thetaStat do_stat STA_theta.thetas.idx -win 50000 -step 10000  -outnames STA.theta


#also this:
#-----------------------
# Computing Fst between pops

# create lists of Orpheus (starting with O) and Keppels (K) samples, names O.pop and K.pop
# based on names of .trim files that are now in your rad directory, using ls and perl -pe commands
ls K*.trim | perl -pe 's/\..+//' | perl -pe 's/2b_//g' >K.pop
ls O*.trim | perl -pe 's/\..+//' >O.pop

vcftools --vcf maxaf.vcf --weir-fst-pop K.pop --weir-fst-pop O.pop
# Weir and Cockerham weighted Fst estimate: 0.01

#### ARCHIVE - OUTDATED ####

##### MDS plot/CCA from IBS matrix #####
library(vegan)
library(adegenet) # for transp()
setwd("~/moorea_holobiont/mr_2brad/3.pop_structure")
setwd("~/moorea_holobiont/mr_2brad/3.pop_structure/tests")

#[from Misha's angsd_ibs_pca.R script]
bams=read.table("bamscl_no7")[,1] # list of bam files
bams=read.table("bams_m")[,1] # list of bam files
goods=c(1:length(bams))

## reading table of pairs of replicates (tab-delimited) - skip if there are no clones
#clonepairs=read.table("clonepairs.tab",sep="\t")
#repsa= clonepairs[,1]
#repsb= clonepairs[,2]
## removing "b" replicates
#goods=which(!(bams %in% repsb))

#--
# loading individual to population correspondences
i2p=read.table("bamscl_no7_pops_m.txt",sep="\t") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
#i2p=read.table("part3_bams_no7_year.txt") #for sequencer instead
row.names(i2p)=i2p[,1]
i2p=i2p[goods,]
site_zone=i2p[,2]
site=i2p[,3]
zone=i2p[,4]

# setting up colors for plotting
palette(rainbow(length(unique(site_zone))))
colors=as.numeric(as.factor(site_zone))
colpops=as.numeric(as.factor(sort(unique(site_zone))))

ma = as.matrix(read.table("result_m.ibsMat"))

#preliminary look
# pca(ma)
# ibs.pca <- prcomp(ma) 
# library(devtools)
# install_github("vqv/ggbiplot")
# library(ggbiplot)
# 
# ggbiplot(ibs.pca,ellipse=TRUE,groups=site)

#---
# clustering / PCoA based on identity by state (IBS) based on single read resampling
# (for low and/or uneven coverage)

# performing PCoA and CAP
conds=data.frame(cbind(site_zone,site,zone))
conds=data.frame(cbind(site_zone))

pp0=capscale(ma~1)
pp=capscale(ma~site_zone,conds)

# significance of by-site divergence
adonis(ma~site_zone,data=conds,permutations = 9999) 
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
# site        2    0.1273 0.063663  1.0937 0.01932  0.011 *
#   zone        1    0.0607 0.060695  1.0427 0.00921  0.242  
# site:zone   2    0.1174 0.058693  1.0083 0.01781  0.407  
# Residuals 108    6.2864 0.058208         0.95367         
# Total     113    6.5918                  1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#this example helped me figure out adonis: https://stackoverflow.com/questions/43736314/error-in-g-that-non-conformable-arrays
#also this: https://rdrr.io/rforge/vegan/man/adonis.html#heading-1

#beta dispersion of by-site divergence
ma.scores <- vegdist(ma)
out <- betadisper(ma.scores,group=conds$site)
anova(out)

## eigenvectors
plot(pp0$CA$eig) 

axes2plot=c(1,2)  
quartz()
cmd=pp0 #change to pp for CAP, pp0 for MDS
plot(cmd,choices=axes2plot,display="sites",type="n") # choices - axes to display
points(cmd,choices=axes2plot,pch=19,col=transp(colors,alpha=0.7))
#ordihull(cmd,choices= axes2plot,groups= conds$grp,draw="polygon",col=1+as.numeric(unique(as.factor(conds$grp))),label=T)
ordispider(cmd,choices= axes2plot,groups=conds$site,col="grey80")
ordiellipse(cmd,choices= axes2plot,groups= conds$site,draw="polygon",col=colpops,label=T)

## unscaled, to identify outliers
#plot(cmd$CA$u[,axes2plot],pch=19,col=colors)
#ordispider(cmd$CA$u[,axes2plot],groups=conds$site,col="grey80")
#ordiellipse(cmd$CA$u[,axes2plot],groups= conds$site,draw="polygon",col=colpops,label=T)
#identify(cmd$CA$u[,axes2plot],labels=colnames(ma),n=3,cex=0.7)

#### ggplot MDS plot ####
library(ggplot2)
#install.packages("remotes")
#remotes::install_github("jfq3/ggordiplots")
library("ggordiplots")
library("stringr")
library("cowplot")
library("extrafont")
loadfonts()

mds <- pp0[["CA"]][["u"]]

p1 <- ggordiplots::gg_ordiplot(ord=mds,groups=site,scaling=0.8,choices=c(1,2),conf=0.7,spiders=FALSE,pt.size=2)
#names(p1)
#spid <- p1$df_spiders
#p1$df_spiders
gg <- p1$df_ord

site <- sub("TO","TNW-F",site) #just had to make all sites have 4 characters so I could look at some other variables
site <- sub("TI","TNW-B",site)
site <- sub("I","-B",site)
site <- sub("O","-F",site)

zone <- str_sub(site,5,5)
realsite <- str_sub(site,1,3)

gg$zone <- zone
gg$reef <- realsite
quartz()
ggplot(gg,aes(x=x,y=y,color=zone,shape=reef,fill=zone))+
  geom_point(size=2)+
  # geom_segment(data=p1$df_spiders, aes(x=x, y=y, xend=cntr.x, yend=cntr.y))+
  stat_ellipse(level=0.8,aes(lty=zone),geom="polygon",alpha=0.1)+
  xlab('MDS1 (5.85%)')+
  ylab('MDS2 (5.76%)')+
  theme_cowplot()+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Backreef","Forereef"))+
scale_shape_manual(values=c(15,16,17),labels=c("MNW","MSE","TNW"))+
  scale_fill_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Backreef","Forereef"))+
  scale_linetype_manual(values=c("solid","twodash"),labels=c("Backreef","Forereef"))+
  labs(shape="Site",color="Reef zone",linetype="Reef zone",fill="Reef zone")+
  theme(text=element_text(family="Gill Sans MT"))

#### K plot from NGSadmix & beagle file ####

# assembling the input table
dir="~/Google Drive/Moorea/2brad_moorea/" # path to input files
inName="part3_don_k2.qopt" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run
pops="part3_inds2pops_no7.txt" # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.

npops=as.numeric(sub(".+(\\d+)\\..+","\\1",inName))
tbl=read.table(paste(dir,inName,sep=""),header=F)
i2p=read.table(paste(dir,pops,sep=""),header=F)

names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
row.names(tbl)=tbl$ind
head(tbl)
# putting populaitons in desired order (edit pop names as needed or skip to plot them alphabetically)
# tbl$pop=factor(tbl$pop,levels=c("O","K"))

source("~/Google Drive/Moorea/2brad_moorea/plot_admixture_v5_function.R")
ords=plotAdmixture(data=tbl,npops=npops,grouping.method="distance",vshift=0.1)

## recording cluster affiliations
#cluster.admix=apply(tbl[,1:npops],1,function(x) {return(paste(which(x>0.25),collapse=".")) })
#save(cluster.admix,file=paste(inName,"_clusters.RData",sep=""))

#### PCA from .vcf file ####

library(vcfR)
library(adegenet)
library(vegan) 
library(stringr)
#install.packages("poppr")
#library("poppr")
setwd("~/moorea_holobiont/mr_2brad/3.pop_structure")

gl=vcfR2genlight(read.vcfR("clresult_no7.vcf")) #output from angsd
#gl=vcfR2genlight(read.vcfR("~/Google Drive/Moorea/2brad_moorea/don_nolink.recode.vcf")) #output from angsd

##some visualizations
#glPlot(gl, posi="topleft")
##white blocks would indicate missing data

##allele frequency spectrum
# myFreq <- glMean(gl)
# hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies",
#      main="Distribution of (second) allele frequencies")
# temp <- density(myFreq)
# lines(temp$x, temp$y*1.8,lwd=3)

# assign populations
#pops=read.table("bamscl_no7_pops_m_noold.txt",sep="\t")
pops=read.table("bamscl_no7_pops.txt",sep="\t") #bams file with a 2nd column describing variable of interest
#pops=read.table("part3_bams_no7_year.txt",sep="\t") #another variable I looked at
pop(gl)=pops$V3 #you can have other columns with other variables in other columns, select which one for analysis here
pca=glPca(gl,nf=3,parallel=F) #make pca

#nice pca
quartz()
plot(pca$scores[,1:2],col=transp(as.numeric(as.factor(pop(gl))),0.3),pch=19)
ordispider(pca$scores[,1:2],pop(gl),col=as.numeric(as.factor(pop(gl))),label=T)
ordiellipse(pca$scores[,1:2],pops$V2,label=T,draw="polygon",col="grey90")

##I didn't do this, but useful:
## manually identifying outliers: click on outlier points in the plot.
## adjust n=3 parameter to identify more than 3 points
#outliers=identify(pca$scores[,1:2],labels=gl@ind.names,n=3,cex=0.8)

## re-making PCA without outliers
#gl2=gl[-outliers]
#pca2=glPca(gl2,nf=3,parallel=F)
#colors=as.numeric(as.numeric(as.factor(levels(pop(gl))))) # or use your own colors for populations
#s.class(pca2$scores[],pop(gl2),col=transp(colors,0.5),cstar=1,cellipse=1,clabel=1,axesell=F,grid=F,cpoint=2)

#### ggplot PCA ####
library(ggplot2)
#install.packages("ggrepel")

p1 <- ggordiplots::gg_ordiplot(ord=scores,groups=site,scaling=0.8,choices=c(1,2),conf=0.7,spiders=FALSE,pt.size=2)
#names(p1)
spid <- p1$df_spiders
#p1$df_spiders
gg <- p1$df_ord
gg$zone <- zone
gg$reef <- realsite
quartz()
ggplot(gg,aes(x=x,y=y,color=zone,shape=reef,fill=zone))+
  geom_point(size=2)+
  # geom_segment(data=p1$df_spiders, aes(x=x, y=y, xend=cntr.x, yend=cntr.y))+
  stat_ellipse(level=0.8,aes(lty=zone),geom="polygon",alpha=0.1)+
  xlab('PC1 (4.11%)')+
  ylab('PC2 (4.06%)')+
  theme_classic()+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Backreef","Forereef"))+
  scale_shape_manual(values=c(15,16,17),labels=c("MNW","MSE","TNW"))+
  scale_fill_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Backreef","Forereef"))+
  scale_linetype_manual(values=c("solid","twodash"),labels=c("Backreef","Forereef"))+
  labs(shape="Site",color="Reef zone",linetype="Reef zone",fill="Reef zone")
# geom_label(x=-0.2618,y=-0.1116,label="MNW-B",color="black",fill="white")+
# geom_label(x=-0.3347,y=-0.1352,label="MNW-F",color="black",fill="white")+
# geom_label(x=-0.1423,y=-0.9173,label="MSE-B",color="black",fill="white")+
# geom_label(x=0.6595,y=0.8860,label="MSE-F",color="black",fill="white")+
# geom_label(x=-0.0550,y=-0.1702,label="TNW-B",color="black",fill="white")+
# geom_label(x=-0.0671,y=0.4204,label="TNW-F",color="black",fill="white")

#now for sequencer instead of sites
years <- read.table("~/Google Drive/Moorea/Host/bamscl_year_no7.txt")
year <- years[,2]
year <- gsub('2015', 'HiSeq 2500',year)
year <- gsub('2017', 'HiSeq 4000',year)

quartz()
ggplot(gg,aes(x=x,y=y,color=year,shape=year,fill=year))+
  geom_point(size=2)+
  # geom_segment(data=p1$df_spiders, aes(x=x, y=y, xend=cntr.x, yend=cntr.y))+
  stat_ellipse(level=0.8,aes(lty=year),geom="polygon",alpha=0.1)+
  xlab('PC1 (4.11%)')+
  ylab('PC2 (4.06%)')+
  theme_classic()+
  labs(shape="Sequencer",color="Sequencer",linetype="Sequencer",fill="Sequencer")
