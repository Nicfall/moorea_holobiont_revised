library("ggplot2")

#setwd("~/nicfall drive/Moorea_revisions/mr2brad_rev/5.diversity")
df <- read.table("part5_het_out_depth.txt") #read in data
df$het <- df$V3/(df$V2+df$V3) #heterozygosity calculations
df$V1 <- sub("TO","TNWO",df$V1) #just renaming some sites
df$V1 <- sub("TI","TNWI",df$V1) #just renaming some sites
df$site <- substr(df$V1, 0, 4)
df$site <- as.factor(df$site)
str(df)
plot(het~site,data=df)

library("bestNormalize")
het.t <- bestNormalize(df$het)
df$het.t <- het.t$x.t
shapiro.test(df$het.t)
hist(df$het.t)

# a1 <- aov(het.t~site,data=df)
# summary(a1)
# TukeyHSD(a1)

mnw <- subset(df,site==c("MNWI","MNWO"))
wilcox.test(het ~ site, data=mnw) #W = 36, p-value = 0.3969

mse <- subset(df,site==c("MSEI","MSEO"))
wilcox.test(het ~ site, data=mse) #W = 63, p-value = 0.1564

tah <- subset(df,site==c("TNWI","TNWO"))
wilcox.test(het ~ site, data=tah) #W = 77, p-value = 0.7987

df$zone <- substr(df$V1, 4, 4)
df$zone <- as.factor(df$zone)
str(df)

df$zone <- sub("I","B",df$zone)
df$zone <- sub("O","F",df$zone)
df$site <- sub("I","-B",df$site)
df$site <- sub("O","-F",df$site)
colorz <- c("#ED7953FF","#8405A7FF")

quartz()
ggplot(df,aes(x=site,y=het,color=zone))+
  geom_boxplot()+
  theme_classic()+
  scale_color_manual(values=colorz)+
  theme(legend.position = "none",text=element_text(family="Times"))+
  ylab("Heterozygosity")+
  xlab("Site")


