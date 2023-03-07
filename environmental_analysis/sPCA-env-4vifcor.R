library(ggplot2)
library(reshape2)
library(data.table)
library(dplyr)
library(adegenet)
library(genepopedit)
library(akima) #masks lazyeval interp function
library(splancs)
library(gridExtra)
library(G)
library(zvau)
library(tidyverse)
library(dplyr)
library(ggthemes)
library(ggrepel)
library(vegan)

# Load coordinates ----

coord.ind <- read.csv("gsa-coord-ind-cart-jitter.csv",header=T)

##Jitter the geographical coordinates so that a Delaueuroy triangulation can be applied. Here we use a jitter of 1 km 
coord.ind$jLong <- jitter(coord.ind$Long,1) 
coord.ind$jLat <- jitter(coord.ind$Lat,1)


## Read in genepop files as GENIND objects ------------
gen.union.4vifcor <- read.genepop("./my_genepop_GSA_union_4vifcor.gen")
levels(pop(gen.union.4vifcor)) <- c("GSA6","GSA7",
                                 "GSA9","GSA10","GSA17",
                                 "GSA18w","GSA18e","GSA22",
                                 "GSA24")

##Add the spatial coordinates -----

env_Cart <- gen.union.4vifcor # cartesian coordieurotes
env_Cart@other$xy <- coord.ind[,c("MDS1","MDS2")]
head(env_Cart@other$xy)
tail(env_Cart@other$xy)

env_jCart <- gen.union.4vifcor # jittered Cartesian coordieurotes (Delaueuroy triangulation)
env_jCart@other$xy <- coord.ind[,c("jCartx","jCarty")]
head(env_jCart@other$xy)
tail(env_jCart@other$xy)

env_Coords <- gen.union.4vifcor # Standard lat long
env_Coords@other$xy <- coord.ind[,c("Long","Lat")]

env_jCoords <- gen.union.4vifcor # Standard lat long
env_jCoords@other$xy <- coord.ind[,c("jLong","jLat")]

####### sPCA ------

#### Delaunay (type = 1)
env_sPCA_jCart <- spca(env_jCart,ask=FALSE,type=1,scannf = T)
## positive 3 negative 3 

# Global test and local test 
randtest <- spca_randtest(env_sPCA_jCart, nperm = 999)
randtest
#"Results show that global test was significant while local test no
# $global
# Monte-Carlo test
# Call: as.randtest(sim = sims[1, ], obs = obs[1], alter = "greater")
# 
# Observation: 5.482947 
# 
# Based on 999 replicates
# Simulated p-value: 0.001 
# Alternative hypothesis: greater 
# 
# Std.Obs.pos  Expectation     Variance 
# 48.618029954  2.199377064  0.004561396 
# 
# $local
# Monte-Carlo test
# Call: as.randtest(sim = sims[2, ], obs = obs[2], alter = "greater")
# 
# Observation: 2.061353 
# 
# Based on 999 replicates
# Simulated p-value: 1 
# Alternative hypothesis: greater 
# 
# Std.Obs.neg  Expectation     Variance 
# -3.728917589  2.271775258  0.003184332 

# re-run sPCA retaing first 3 positive axes
env_sPCA_jCart <- spca(env_jCart,ask=FALSE,type=1,scannf = T)
barplot(env_sPCA_jCart$eig,col=spectral(length(env_sPCA_jCart$eig)), ylim = c(-0.05,2.5))
legend(75, 2.7, pt.bg  = spectral(2), cex = 1.2, bty = "n", xpd = T, pch = 21, y.intersp = 0.4,
       leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")

screeplot(env_sPCA_jCart)

# Interpolate sPCA to identify clines
png("env_4vifcor_SPCA_Delaunay_interpolatedMap.png",width=2500,height=2000,res=300)
x <- other(env_Cart)$xy[,1]
y <- other(env_Cart)$xy[,2]
interpX <- seq(min(x),max(x),le=200)
interpY <- seq(min(y),max(y),le=200)
temp <- interp(x, y, env_sPCA_jCart$ls[,1], xo=interpX, yo=interpY, duplicate = "mean")
myPal <- colorRampPalette(c("#8C2A1C", "#F7EC16", 
                            "azure", "#11A4C8", "#2A2771"))
annot <- function(){
  text(unique(x), unique(y), labels=gsa, pos=c(2,1,3,1,2,3,2,1,1))
  #points(x,y)
  title(xlab = "Cartesian 1", ylab = "Cartesian 2")
  #axis(side=2, at=seq(-800, 1000, by=200))
}
filled.contour(temp, color.pal=myPal, nlev=50, xlim=c(-1900, 1900), ylim=c(-800, 1100),
               key.title=title("lagged \nscore 1"), plot.title=annot())
dev.off()

plot(env_sPCA_jCart$ls[,1], coord.ind$Long)
lines(env_sPCA_jCart$ls[,1], coord.ind$Long)

A <- data.frame(Code = pop(gen.union.4vifcor), Long = coord.ind[,2], Scores =env_sPCA_jCart$ls[,1])
B <- A %>% group_by(Code) %>% mutate(meanLagS = mean(Scores))
B <- as.data.frame(B)
B$revLong <- rev(B$Long)
B$revCode <- factor(rev(B$Code), levels= rev(levels(B$Code)))
B$revMean <- rev(B$meanLagS)
bg <- c("#11A4C8", "#2A2771",
        "#396D35", "#808080",
        "#ED2224", "#ED3995","#7E277C",
        "#F8941E","#F7EC16")

ggplot(B, aes(meanLagS, Long)) + 
  geom_smooth(data = B, method = "loess", se = F, col= "black", orientation = "x", span = 0.5, size = 0.5) + 
  geom_text(data = B, aes(label = Code), size = 5, position = position_nudge(y = 1.7))+
  geom_point(data = B, aes(fill = Code), colour = "grey32", size = 6, shape = 21, show.legend = F) +
  scale_color_manual(values = bg)+
  scale_fill_manual(values = bg)+
  xlab("sPCA axis 1 lagged scores") + 
  ylab("Longitude °E")+
  #scale_y_reverse() +
  scale_x_reverse()+
  theme_base(15)

## we retained the first positive eigenvalue extracted from sPCA
locality.scores <- env_sPCA_jCart$ls[,1]

##linear regression 
locality.env <- read.csv("FourEnv-4vifcor-dbMEM.csv", header = T)
locality.env <- cbind(locality.env, coord.ind[,2:3])
locality.mod <-setNames(as.list(rep(NA,length(locality.env))), colnames(locality.env))
for (i in 1:length(locality.env)){
  lm.tmp <- lm(locality.scores ~ locality.env[,i], data = locality.env)
  locality.mod[[i]] <- lm.tmp
}
map(locality.mod, summary)
# $dbMEM1
# 
# Call:
#   lm(formula = locality.scores ~ locality.env[, i], data = locality.env)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.09911 -0.53516  0.02935  0.40169  1.95739 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.003351   0.043029   0.078    0.938    
# locality.env[, i] -27.550880   0.807295 -34.127   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.8073 on 350 degrees of freedom
# Multiple R-squared:  0.7689,	Adjusted R-squared:  0.7683 
# F-statistic:  1165 on 1 and 350 DF,  p-value: < 2.2e-16
# 
# 
# $dbMEM2
# 
# Call:
#   lm(formula = locality.scores ~ locality.env[, i], data = locality.env)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.6942 -1.7296  0.3509  1.5417  2.5785 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)       0.003351   0.089387   0.037    0.970
# locality.env[, i] 1.666090   1.677054   0.993    0.321
# 
# Residual standard error: 1.677 on 350 degrees of freedom
# Multiple R-squared:  0.002812,	Adjusted R-squared:  -3.713e-05 
# F-statistic: 0.987 on 1 and 350 DF,  p-value: 0.3212
# 
# 
# $dbMEM3
# 
# Call:
#   lm(formula = locality.scores ~ locality.env[, i], data = locality.env)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.5122 -1.9277  0.6011  1.0509  2.5631 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.003351   0.083880    0.04    0.968    
# locality.env[, i] -10.969470   1.573736   -6.97 1.58e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.574 on 350 degrees of freedom
# Multiple R-squared:  0.1219,	Adjusted R-squared:  0.1194 
# F-statistic: 48.59 on 1 and 350 DF,  p-value: 1.58e-11
# 
# 
# $dbMEM4
# 
# Call:
#   lm(formula = locality.scores ~ locality.env[, i], data = locality.env)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.6899 -1.8851  0.5581  1.4884  2.4928 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)       0.003351   0.089346   0.038    0.970
# locality.env[, i] 1.917418   1.676287   1.144    0.253
# 
# Residual standard error: 1.676 on 350 degrees of freedom
# Multiple R-squared:  0.003724,	Adjusted R-squared:  0.0008778 
# F-statistic: 1.308 on 1 and 350 DF,  p-value: 0.2535
# 
# 
# $dbMEM5
# 
# Call:
#   lm(formula = locality.scores ~ locality.env[, i], data = locality.env)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.7186 -1.8366  0.3306  1.4822  2.4888 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       0.003351   0.087262   0.038    0.969    
# locality.env[, i] 7.002316   1.637177   4.277 2.45e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.637 on 350 degrees of freedom
# Multiple R-squared:  0.04967,	Adjusted R-squared:  0.04696 
# F-statistic: 18.29 on 1 and 350 DF,  p-value: 2.446e-05
# 
# 
# $curM_Wclm_mean
# 
# Call:
#   lm(formula = locality.scores ~ locality.env[, i], data = locality.env)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -2.769 -1.522  0.257  1.442  2.306 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       -0.15716    0.08961  -1.754   0.0803 .  
# locality.env[, i] 85.16573   14.44965   5.894 8.87e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.602 on 350 degrees of freedom
# Multiple R-squared:  0.09029,	Adjusted R-squared:  0.08769 
# F-statistic: 34.74 on 1 and 350 DF,  p-value: 8.871e-09
# 
# 
# $temp_bottom_range
# 
# Call:
#   lm(formula = locality.scores ~ locality.env[, i], data = locality.env)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.6165 -1.6376  0.6616  1.3085  2.0697 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        2.02708    0.24921   8.134 7.31e-15 ***
#   locality.env[, i] -0.27531    0.03205  -8.591 2.89e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.526 on 350 degrees of freedom
# Multiple R-squared:  0.1742,	Adjusted R-squared:  0.1718 
# F-statistic: 73.81 on 1 and 350 DF,  p-value: 2.891e-16
# 
# 
# $sal_bottom_mean
# 
# Call:
#   lm(formula = locality.scores ~ locality.env[, i], data = locality.env)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.3222 -0.7602  0.3577  0.7464  1.5746 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       104.7499     4.3004   24.36   <2e-16 ***
#   locality.env[, i]  -2.7330     0.1122  -24.36   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.023 on 350 degrees of freedom
# Multiple R-squared:  0.629,	Adjusted R-squared:  0.6279 
# F-statistic: 593.4 on 1 and 350 DF,  p-value: < 2.2e-16
# 
# 
# $pco_Wclm_mean
# 
# Call:
#   lm(formula = locality.scores ~ locality.env[, i], data = locality.env)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.0939 -1.1891  0.5953  1.3916  2.2029 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       -9.090984   2.141631  -4.245 2.81e-05 ***
#   locality.env[, i]  0.024825   0.005841   4.250 2.75e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.638 on 350 degrees of freedom
# Multiple R-squared:  0.04907,	Adjusted R-squared:  0.04636 
# F-statistic: 18.06 on 1 and 350 DF,  p-value: 2.745e-05
# 
# 
# $Long
# 
# Call:
#   lm(formula = locality.scores ~ locality.env[, i], data = locality.env)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.56498 -0.30803 -0.01204  0.26143  1.38139 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        2.083491   0.047342   44.01   <2e-16 ***
#   locality.env[, i] -0.132818   0.002404  -55.26   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5386 on 350 degrees of freedom
# Multiple R-squared:  0.8972,	Adjusted R-squared:  0.8969 
# F-statistic:  3053 on 1 and 350 DF,  p-value: < 2.2e-16
# 
# 
# $Lat
# 
# Call:
#   lm(formula = locality.scores ~ locality.env[, i], data = locality.env)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.60799 -1.01831  0.02073  1.01939  2.73195 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       -15.40114    1.11120  -13.86   <2e-16 ***
#   locality.env[, i]   0.37584    0.02705   13.89   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.348 on 350 degrees of freedom
# Multiple R-squared:  0.3554,	Adjusted R-squared:  0.3536 
# F-statistic:   193 on 1 and 350 DF,  p-value: < 2.2e-16

map(locality.mod, RsquareAdj)
# $dbMEM1
# $dbMEM1$r.squared
# [1] 0.7689281
# 
# $dbMEM1$adj.r.squared
# [1] 0.7682679
# 
# 
# $dbMEM2
# $dbMEM2$r.squared
# [1] 0.002811976
# 
# $dbMEM2$adj.r.squared
# [1] -3.713233e-05
# 
# 
# $dbMEM3
# $dbMEM3$r.squared
# [1] 0.121895
# 
# $dbMEM3$adj.r.squared
# [1] 0.1193862
# 
# 
# $dbMEM4
# $dbMEM4$r.squared
# [1] 0.00372433
# 
# $dbMEM4$adj.r.squared
# [1] 0.0008778286
# 
# 
# $dbMEM5
# $dbMEM5$r.squared
# [1] 0.04967046
# 
# $dbMEM5$adj.r.squared
# [1] 0.04695523
# 
# 
# $curM_Wclm_mean
# $curM_Wclm_mean$r.squared
# [1] 0.09029194
# 
# $curM_Wclm_mean$adj.r.squared
# [1] 0.08769277
# 
# 
# $temp_bottom_range
# $temp_bottom_range$r.squared
# [1] 0.1741513
# 
# $temp_bottom_range$adj.r.squared
# [1] 0.1717917
# 
# 
# $sal_bottom_mean
# $sal_bottom_mean$r.squared
# [1] 0.6289958
# 
# $sal_bottom_mean$adj.r.squared
# [1] 0.6279358
# 
# 
# $pco_Wclm_mean
# $pco_Wclm_mean$r.squared
# [1] 0.04907422
# 
# $pco_Wclm_mean$adj.r.squared
# [1] 0.04635728
# 
# 
# $Long
# $Long$r.squared
# [1] 0.8971531
# 
# $Long$adj.r.squared
# [1] 0.8968593
# 
# 
# $Lat
# $Lat$r.squared
# [1] 0.3554181
# 
# $Lat$adj.r.squared
# [1] 0.3535764


# # Which alleles exhibit structure (squared loading)
myLoadings <- env_sPCA_jCart$c1[,1]^2
names(myLoadings) <- rownames(env_sPCA_jCart$c1)

#to identify the 5% of alleles with the greatest contributions to the first global structure in Spca result
temp <- loadingplot(myLoadings, threshold=quantile(myLoadings, 0.95),
                    xlab="Alleles",ylab="Weight of the alleles",
                    main="Contribution of alleles \n to the first sPCA axis",
                    cex.fac=0.6, srt = 45, lab.jitter = 35, adj = c(0.7, 0.3))
save.image("sPCA-Data-4vifcor.RData")
