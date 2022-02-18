# library(MCMC.OTU)
# counts.nums <- 1:nrow(counts.no0) #how many samples you have (i.e. 1-93)
# 
# counts.int <- cbind(sample = 0, counts.no0)
# counts.formcmc <- cbind(X = 0, counts.int)
# 
# counts.formcmc$X <- counts.nums
# counts.formcmc$sample <- row.names(counts.no0)
# row.names(counts.formcmc) <- counts.formcmc$X
# 
# goods <- purgeOutliers(counts.formcmc,count.columns=3:11,sampleZcut=(-2.5))
