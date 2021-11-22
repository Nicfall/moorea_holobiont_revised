#mostly from https://github.com/z0on/2bRAD_denovo/blob/master/2bRAD_README.txt
#^^find instructions for downloading scripts & packages at this link^^
#edits by Nicola Kriefall thenicolakriefall(at)gmail.com

#### perform PCA on covarince matrix from pcangsd ####
library(adegenet)
library(vegan)

setwd("/Users/nicolakriefall/moorea_holobiont/mr_2brad/3.pop_structure")
C <- as.matrix(read.table("clresult_no7.cov"))
e <- eigen(C)

i2p=read.table("bamscl_no7_pops.txt",sep="\t") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
site_zone=i2p[,2]
site=i2p[,3]
zone=i2p[,4]

palette(rainbow(length(unique(site_zone))))
colors=as.numeric(as.factor(site_zone))
colpops=as.numeric(as.factor(sort(unique(site_zone))))

conds=data.frame(cbind(site_zone,site,zone))

plot(e$vectors[,1:2],xlab="PC1",ylab="PC2",col=transp(colors))
ordiellipse(e$vectors[,1:2],groups=conds$site_zone,draw="polygon",col=colpops,label=T)

adonis(C~site*zone,data=conds)
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# site        2    -1.132 -0.56590 0.99628 0.01764  0.983
# zone        1    -0.567 -0.56745 0.99901 0.00884  0.783
# site:zone   2    -1.132 -0.56582 0.99614 0.01763  0.961
# Residuals 108   -61.345 -0.56801         0.95589       
# Total     113   -64.176                  1.00000       

#eigenvalues
rbind(
  SD = sqrt(e[["values"]]),
  Proportion = e[["values"]]/sum(e[["values"]]),
  Cumulative = cumsum(e[["values"]])/sum(e[["values"]]))
#1 = 0.01302991
#2 = 0.01247939

#ggplot
library(cowplot)
library(ggplot2)

scores <- e$vectors[,1:2]
df <- cbind(scores,i2p)
df <- unname(df)
colnames(df) <- c("score1","score2","bam","site_zone","site","zone")

#site & reef zone
quartz()
ggplot(df,aes(x=score1,y=score2,color=zone,fill=zone,shape=site))+
  geom_point(size=2)+
  stat_ellipse(level=0.8,aes(lty=zone),geom="polygon",alpha=0.1)+
  xlab('PC1 (1.3%)')+
  ylab('PC2 (1.2%)')+
  theme_cowplot()+  
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("BR","FR"))+
  scale_shape_manual(values=c(8,4,9),labels=c("MNW","MSE","TNW"))+
  scale_fill_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("BR","FR"))+
  scale_linetype_manual(values=c("solid","twodash"),labels=c("BR","FR"))+
  labs(shape="Site",color="Reef zone",linetype="Reef zone",fill="Reef zone")

#site
quartz()
ggplot(df,aes(x=score1,y=score2,color=site,fill=site,shape=site))+
  geom_point(size=2)+
  stat_ellipse(level=0.8,aes(lty=site),geom="polygon",alpha=0.1)+
  xlab('PC1 (1.3%)')+
  ylab('PC2 (1.2%)')+
  theme_cowplot()+  
  scale_linetype_manual(values=c("longdash","dotted","dotdash"),labels=c("MNW","MSE","TNW"))+
  scale_color_manual(values=c("darkslategray3","darkslategray4","#000004"),labels=c("MNW","MSE","TNW"))+
  scale_shape_manual(values=c(8,4,9),labels=c("MNW","MSE","TNW"))+
  scale_fill_manual(values=c("darkslategray3","darkslategray4","#000004"),labels=c("MNW","MSE","TNW"))+
  labs(shape="Site",color="Site",linetype="Site",fill="Site")

#### K plot from .vcf ####

# primitive look at admixture data:
tbl=read.table("clresult_no7.2.Q")
barplot(t(as.matrix(tbl)), col=rainbow(5),xlab="Individual #", ylab="Ancestry", border=NA)

#---
# prettier:

# assembling the input table
dir="~/moorea_holobiont/mr_2brad/3.pop_structure/" # path to input files
inName="clresult_no7.2.Q" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run
npops=2
pops="bamscl_no7_pops.txt" # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
tbl=read.table(paste(dir,inName,sep=""),header=F)
i2p=read.table(paste(dir,pops,sep=""),header=F)
#i2p <- i2p[,1:2]
names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
row.names(tbl)=tbl$ind

head(tbl,20) # this is how the resulting dataset must look

source("~/moorea_holobiont/mr_2brad/3.pop_structure/plot_admixture_v4_function copy.R")

# putting populaitons in desired order (edit pop names as needed or skip to plot them alphabetically)
tbl$pop=factor(tbl$pop,levels=c("MNW-B","MNW-F","MSE-B","MSE-F","TNW-B","TNW-F"))

quartz()
ords=plotAdmixture(data=tbl,npops=npops,angle=0,vshift=0,hshift=0)

#### FST from .vcf file ####
#good webiste on calculating this stuff (by hand though):
#http://www.uwyo.edu/dbmcd/popecol/maylects/fst.html

library(vcfR)
#install.packages("hierfstat")
library("hierfstat")
library(adegenet)

setwd("~/")
vcf <- read.vcfR("~/moorea_holobiont/mr_2brad/3.pop_structure/clresult_no7.vcf")
genind <- vcfR2genind(vcf)
pop(genind) <- i2p$V2

basic.stats(genind)
# $overall
# Ho     Hs     Ht    Dst    Htp   Dstp    Fst   Fstp    Fis   Dest 
# 0.3167 0.3410 0.3413 0.0003 0.3414 0.0004 0.0010 0.0012 0.0712 0.0006 

fstat(genind, fstonly = FALSE, pop=NULL) #pop=null means you inherit the pops already there 
pairwise.fst(genind, pop = NULL, res.type = c("dist", "matrix"))

test.between(genind,test.lev="Locality",rand.unit="Patch")

# 1          2          3          4          5
# 2 0.01815735                                            
# 3 0.01772569 0.01547715                                 
# 4 0.01846476 0.01551543 0.01605031                      
# 5 0.01549681 0.01305472 0.01372456 0.01394897           
# 6 0.01571488 0.01373486 0.01399234 0.01395402 0.01154477

