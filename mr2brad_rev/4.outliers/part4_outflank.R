#Following instructions from:
#https://htmlpreview.github.io/?https://github.com/whitlock/OutFLANK/blob/master/inst/doc/OutFLANKAnalysis.html
#2018-04-22 version

library(devtools)
library("vcfR")
#install.packages("https://cran.r-project.org/src/contrib/Archive/raster/raster_2.5-8.tar.gz", repos = NULL, type="source")
#devtools::install_github("whitlock/OutFLANK")
library("OutFLANK")

#https://datadryad.org/bitstream/handle/10255/dryad.125768/script-outflank.R?sequence=1

###Prepare the input file in terminal
##Convert vcf file to 012 file
#vcftools --vcf mnw_hwe.vcf --012 --out mnw_hwe_out
##Replace the -1 by 9 in the 012 file

##Remove the first column
#awk '{$1=""; print substr($0,1)}' mnw_hwe_out.012 > SNPmat_mnw.txt
##Transfered SNPmat.txt to my working directory 
setwd("~/Google Drive/Moorea/outflank/")
SNPmat <- read.table("SNPmat_mnw.txt",header=FALSE)

##now for the locus names
#module load bcftools
#module load htslib
#bgzip mnw_hwe.vcf
#bcftools tabix mnw_hwe.vcf.gz  
#bcftools annotate --set-id '%CHROM\_%POS\' mnw_hwe.vcf.gz > mnw_hwe_id.vcf
#cut -f 3 mnw_hwe_id.vcf > loci_name_mnw.txt
##transferred to working directory

locusNames <- read.table("loci_name_mnw.txt",header=TRUE)

popNames <- read.table("pops_mnw.txt",header=FALSE)

fstmat <- MakeDiploidFSTMat(SNPmat,locusNames,popNames)

results <- OutFLANK(fstmat, LeftTrimFraction=0.01,
         RightTrimFraction=0.01, Hmin=0.05, 114,
         qthreshold=0.05)

OutFLANKResultsPlotter(results, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.05, binwidth = 0.005, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
