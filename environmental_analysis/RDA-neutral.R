library(adegenet)
library(adegraphics)
library(adespatial)
library(ape)
library(car)
library(codep)
library(diveRsity)
library(purrr)
library(vegan)
library(usdm)

## Data import----
# genetic
genpop_neut <- read.genepop("./my_genepop-GSA-neutral.gen")
# environmental
Env1 <- read.csv("./env-data-paper-gsa.csv", header=TRUE)
Env <- read.csv("./env-data-paper-new-order.csv", header=TRUE) 



# # check order individuals in the file
# Order <- data.frame(ID_indiv=gsub("_", "-",indNames(genpop_solea)))
# 
# p <- map(popNames(genpop_solea), ~strsplit(.x, "_"))
# p <- map(p, ~.x[[1]][2])
# p <- unlist(p)
# 
# for (i in 1:13){
#   print(p[i])
#   ord_tmp <- identical(Order$ID_indiv[grep(p[i],Order$ID_indiv)],  Env$ID_indiv[grep(p[i],Env$ID_indiv)])
#   print(ord_tmp)
# } # individuals are in the same order in genetics and env data


## Crepare the response variables----
# extract matrix of allele frequency
afMat <- tab(genpop_neut, freq=TRUE, NA.method = "mean")
# hellinger transformation
afMatH <- decostand(afMat, method = "hellinger")
# perform a PCA on the standardized matrix
Pca <- prcomp(afMatH,center = T, scale. = F)
s <- summary(Pca)
#View(s$importance*100)
# explained variance by component (you can get this info also from "importance" in summary() object)
round((Pca$sdev^2/sum(Pca$sdev^2))*100, 2)[1:5]
# cumulative explained variance
cump <- cumsum(round((Pca$sdev^2/sum(Pca$sdev^2))*100, 2))
length(which(cump <= 80))
# extract the first Pca principal components that cumulatively explained 80% of the total variation (Selmoni2020). Use them as response variable in the RDA
Y <- Pca$x[,1:length(which(cump <= 80))]
plot(Y[,1], Y[,2])


## Create the explanatory variables----
# create Moran Eigenvector's Maps (MEMs) to be the spatial variables
Coor <- Env[,3:2]
Coorxy <- Env[,2:3]
plot(Coor, asp=1)
# Compute spatial distances among sites accounting for the earth curvature
DistSpatial <- gcd.hf(Coor)
# Compute distance-based MEM (follow Selmoni) 
dbMEM = as.data.frame(pcnm(DistSpatial)$vectors)
colnames(dbMEM) = paste0('dbMEM', 1:ncol(dbMEM))

# Select explanatory variables -----
# We tested collinearity 
#### selection of env data
v1 <- vifcor(Env[,-c(1:3)], th=0.65) 
Env <- exclude(Env, v1)
v2 <- vifcor(dbMEM, th=0.65)
#No variable from the 5 input variables has collinearity problem.
# combine dbMEM and environmental variable as an unique dataset named Y used as explanatory in the RDA
X=cbind(dbMEM, Env)



## Compute the Redundancy Analysis----
# GLOBAL model----
# RDA with all explanatory variables 
#set.seed(19215) ## with original vifcor env data
set.seed(18884)
modG = rda(Y ~ ., X, scale = F)
G = summary(modG)

# Global test of the RDA result
anova(modG) # global significance
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = Y ~ dbMEM1 + dbMEM2 + dbMEM3 + dbMEM4 + dbMEM5 + curM_Wclm_mean + sal_bottom_mean + temp_bottom_range + pco_Wclm_mean, data = X, scale = F)
# Df Variance      F Pr(>F)    
# Model      8 0.016358 5.8584  0.001 ***
#   Residual 343 0.119717                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Test of all variables
anova(modG, by="term") # significance by variable
# Permutation test for rda under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = Y ~ dbMEM1 + dbMEM2 + dbMEM3 + dbMEM4 + dbMEM5 + curM_Wclm_mean + sal_bottom_mean + temp_bottom_range + pco_Wclm_mean, data = X, scale = F)
# Df Variance       F Pr(>F)    
# dbMEM1              1 0.009035 25.8857  0.001 ***
#   dbMEM2              1 0.001694  4.8542  0.001 ***
#   dbMEM3              1 0.001493  4.2786  0.001 ***
#   dbMEM4              1 0.000842  2.4122  0.001 ***
#   dbMEM5              1 0.000888  2.5455  0.001 ***
#   curM_Wclm_mean      1 0.000550  1.5745  0.007 ** 
#   sal_bottom_mean     1 0.001200  3.4379  0.001 ***
#   temp_bottom_range   1 0.000656  1.8786  0.004 ** 
#   Residual          343 0.119717                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Tests of all canonical axes
anova(modG, by="axis") # significance by RDA
# Permutation test for rda under reduced model
# Forward tests for axes
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = Y ~ dbMEM1 + dbMEM2 + dbMEM3 + dbMEM4 + dbMEM5 + curM_Wclm_mean + sal_bottom_mean + temp_bottom_range + pco_Wclm_mean, data = X, scale = F)
# Df Variance       F Pr(>F)    
# RDA1       1 0.010912 31.2648  0.001 ***
#   RDA2       1 0.001925  5.5165  0.001 ***
#   RDA3       1 0.001015  2.9089  0.001 ***
#   RDA4       1 0.000765  2.1932  0.001 ***
#   RDA5       1 0.000570  1.6339  0.002 ** 
#   RDA6       1 0.000475  1.3622  0.049 *  
#   RDA7       1 0.000427  1.2220  0.148    
# RDA8       1 0.000267  0.7657  0.975    
# Residual 343 0.119717                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

RsquareAdj(modG) 
# $r.squared
# [1] 0.1202133
# 
# $adj.r.squared
# [1] 0.0996935


### the proportions of accumulated constrained eigenvalues 
round(G$concont$importance[3,"RDA2"]*RsquareAdj(modG)[[2]]*100,2)
# [1] 7.82



# ORDISTEP model----
# variance inflation factors in the RDA: variables X with high VIFs should generally not be manually removed before the application of a procedure of selection of variables (ordistep).
vif.cca(modG)
# dbMEM1            dbMEM2            dbMEM3            dbMEM4            dbMEM5 
# 30.778516          8.046730         16.163141          8.517340          3.450066 
# curM_Wclm_mean   sal_bottom_mean temp_bottom_range     pco_Wclm_mean 
# 1.669041         54.284101          7.645620                NA 

# Forward selection of explanatory variables to reduce the number of explanatory variables.
#set.seed(19456); 
mod0 = rda(Y ~ 1, data = X, scale = T)  
RDA = ordistep(mod0, scope = formula(modG),perm.max=200, direction = "forward")
S=summary(RDA)

anova(RDA)
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = Y ~ sal_bottom_mean + pco_Wclm_mean + dbMEM1 + dbMEM2 + temp_bottom_range + dbMEM3, data = X, scale = T)
# Df Variance      F Pr(>F)    
# Model      6    3.783 1.9043  0.001 ***
#   Residual 345  114.217                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

RDA$anova
anova(RDA, by="term")
# Permutation test for rda under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = Y ~ sal_bottom_mean + pco_Wclm_mean + dbMEM1 + dbMEM2 + temp_bottom_range + dbMEM3, data = X, scale = T)
# Df Variance      F Pr(>F)    
# sal_bottom_mean     1    0.891 2.6922  0.001 ***
#   pco_Wclm_mean       1    0.692 2.0901  0.001 ***
#   dbMEM1              1    0.591 1.7849  0.001 ***
#   dbMEM2              1    0.585 1.7656  0.001 ***
#   temp_bottom_range   1    0.577 1.7435  0.001 ***
#   dbMEM3              1    0.447 1.3495  0.001 ***
#   Residual          345  114.217                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(RDA, by="axis")

RsquareAdj(RDA) 
# $r.squared
# [1] 0.03205663
# 
# $adj.r.squared
# [1] 0.01522283


### the proportions of accumulated constrained eigenvalues
round(S$concont$importance[3,"RDA2"]*RsquareAdj(RDA)[[2]]*100,2)
# [1] 0.7


## Plot----
eco <- read.csv("./Area_coordinates.csv", header = T)
eco <- data.frame(IND = Env1$ID,
                  STRATA = rep(eco$area_name, c(40,73,36,30,40,16,11,40,66)))
eco$STRATA <- factor(eco$STRATA, levels = levels(eco$STRATA)[c(4,6,8,9,1,2,3,5,7)])

### 6 nice colors for the MPA
bg <- c("#11A4C8", "#2A2771",
        "#396D35", "#808080",
        "#ED2224", "#ED3995","#7E277C",
        "#F8941E","#F7EC16")
leg.text <- c("GSA6","GSA7",
              "GSA9","GSA10","GSA17",
              "GSA18w","GSA18e","GSA22",
              "GSA24")

#PLOT WITH ALL VARIABLES

#RDA1-RDA2 save as US LEGAL 8x14
#global----
pdf("RDA-global-neutral-4vifcor.pdf")
par(mar = c(5, 5, 1, 7) + 0.1)
plot(modG, xlab=paste0('RDA1 (',round(G$concont$importance[2,1]*100,1),'%)'), ylab=paste0('RDA2 (',round(G$concont$importance[2,2]*100,1),'%)'), type = "n", scaling=2, display=c("lc","wa"), cex.lab=1.2, cex.axis=1.2)
points(modG, display="wa", cex=1.6, scaling=2, col = "gray32" ,pch=21, bg=bg[as.numeric(eco$STRATA)]) # the pops
# points(RDA, display="lc", pch=20, cex=3, col="red", scaling=2) CENTROID 
text(modG, scaling=2, display="bp", col="darkcyan", cex=1.3, head.arrow = 0.1,
     font=2, arrow.mul = 0.55, adj=1) # the predictors
text(0.6, -0.38, expression("R"["adj"]^2*"=0.0997"), cex=1.2)
text(0.608,-.41,"p-value=0.001", cex=1.2)
#text(-0.3,0.1,"GSA6", col="#11A4C8")
# text(-0.3,0.08,"GSA7",col="#2A2771", cex=1.1)
legend(0.7,0.1, legend=leg.text, bty="n", col="gray32", pch=21, cex=1.3, pt.bg=bg,y.intersp = 0.5, x.intersp = 0.5, xpd=T)
dev.off()

# RDA1-RDA3
pdf("RDA-global-13-neutral.pdf")
par(mar = c(5, 4, 3, 7) + 0.1)
plot(modG, xlab=paste0('RDA1 (',round(G$concont$importance[2,1]*100,1),'%)'), ylab=paste0('RDA3 (',round(G$concont$importance[2,3]*100,1),'%)'), type = "n", scaling=2, choices=c(1,3), display=c("lc","wa"), cex.lab=1.2, cex.axis=1.2)
points(modG, display="wa", cex=1.8, scaling=2, col = "gray32" ,pch=21, bg=bg[as.numeric(eco$STRATA)], choices=c(1,3)) # the pops
# points(RDA, display="lc", pch=20, cex=3, col="red", scaling=2) CENTROID 
text(modG, scaling=2, display="cn", col="darkcyan", cex=1.1, choices=c(1,3), arrow.mul = 3, head.arrow = 0.1)                           # the predictors
legend(4,1.4, xpd=T, legend=leg.text, bty="n", col="gray32", pch=21, cex=1.2, pt.bg=bg)
dev.off()



# PLOT WITH VARIABLES IDENTIFIED BY ORDISTEP
#ordistep----
pdf("RDA-ordistep-neutral-4vifcor.pdf")
par(mar = c(5, 5, 1, 7))
plot(RDA, xlab=paste0('RDA1 (',round(S$concont$importance[2,1]*RsquareAdj(RDA)[[2]]*100,2),'%)'), ylab=paste0('RDA2 (',round(S$concont$importance[2,2]*RsquareAdj(RDA)[[2]]*100,2),'%)'), type="none", cex.lab=1.2, cex.axis=1.2, scaling=2, display=c("lc","wa"))
points(RDA, display="wa", cex=1.6, scaling=2, col = "gray32", bg = bg[as.numeric(eco$STRATA)],pch=21)
# points(RDA, display="lc", pch=20, cex=3, col="red", scaling=2) CENTROID 
text(RDA, scaling=2, display="bp", col="darkcyan", cex=1.2, arrow.mul = 1.5, head.arrow = 0.1, adj=0.5, font=2)
text(3, -1.98, expression("R"["adj"]^2*"=0.0152"), cex=1.2)
text(3.05,-2.13,"p-value=0.001", cex=1.2)# the predictors
legend(2.7,1, legend=leg.text, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg,y.intersp = 0.35, x.intersp = 0.3, xpd=T, pt.cex=1.5)
dev.off()

# ordistep ordisurf----
op <- ordiplot(RDA, choices=c(1,2))

pdf("ordisurf-ordistep-neutral-4vifcor-temp.pdf")
pdf("ordisurf-ordistep-neutral-4vifcor-sal.pdf")
pdf("ordisurf-ordistep-neutral-4vifcor-pco.pdf")

par(mar = c(5, 5, 1, 6) + 0.5)
plot(RDA, xlab=paste0('RDA1 (',round(S$concont$importance[2,1]*RsquareAdj(RDA)[[2]]*100,2),'%)'), ylab=paste0('RDA2 (',round(S$concont$importance[2,2]*RsquareAdj(RDA)[[2]]*100,2),'%)'), type="none", cex.lab=1.2, cex.axis=1.2, scaling=2, display=c("lc","wa"))
points(RDA, display="wa", cex=1.6, scaling=2, col = "gray32", pch=c(0:9)[as.numeric(eco$STRATA)])
text(3, -1.98, expression("R"["adj"]^2*"=0.0152"), cex=1.2)
text(3.05,-2.13,"p-value=0.001", cex=1.2)# the predictors
legend(2.5,1.7, legend=leg.text, bty="n", col="gray32", pch=0:9, cex=1,y.intersp = 0.12, x.intersp = 0.12, xpd=T, pt.cex=1.5)

##temp-----
#add gradient
temp.surface.RDA <- ordisurf(op ~ temp_bottom_range,X, add=T,
                             col = "red3", bubble=T, labcex=1.2,
                             lwd.cl=1.5)
#add arrow
temp.fit.RDA <- envfit(RDA ~ temp_bottom_range, X)
plot(temp.fit.RDA, col="darkcyan", font=2, cex=1.3, arrow.mul=1.7)

##sal----
sal.surface.RDA <- ordisurf(op ~ sal_bottom_mean,X, add=T,
                            col = "darkorchid3", labcex=1.2, bubble=T,
                            lwd.cl=1.5)
sal.fit.RDA <- envfit(RDA ~ sal_bottom_mean, X)
plot(sal.fit.RDA, col="darkcyan", font=2, cex=1.3, adj=1, arrow.mul=1.3)
##pco----
pco.surface.RDA <- ordisurf(op ~ pco_Wclm_mean,X, add=T,
                            col = "coral3", labcex=1.2, bubble=T,
                            lwd.cl=1.5)
pco.fit.RDA <- envfit(RDA ~ pco_Wclm_mean, X)
plot(pco.fit.RDA, col="darkcyan", font=2, cex=1.3, arrow.mul = 1)
dev.off()



#varpart----
RDA$call ### use ordistep variables
env.sub <- Env[,-1]
db.sub <- X[,c(1:3)]
all.sub <- cbind(env.sub, db.sub)
RDA.var <- varpart(Y,env.sub,db.sub, permutations = 1000)
showvarparts(2)
pdf("./varpart-whole.pdf", width = 5.83,height = 8.57)
plot(RDA.var, digits = 2, Xnames = c("Environmental","dbMEM"), bg = c('navy', 'tomato')) 
#save 10x12
dev.off()

## PARTIAL Redundancy Analysis----
# GLOBAL model----
# RDA with all explanatory variables 
set.seed(231); modG = rda(Y ~ curM_Wclm_mean + sal_bottom_mean + temp_bottom_range + pco_Wclm_mean + Condition(dbMEM1 + dbMEM2 + dbMEM3 + dbMEM4 + dbMEM5), data = X, scale = T)
G = summary(modG)

# Global test of the RDA result
anova(modG) # global significance
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = Y ~ curM_Wclm_mean + sal_bottom_mean + temp_bottom_range + pco_Wclm_mean + Condition(dbMEM1 + dbMEM2 + dbMEM3 + dbMEM4 + dbMEM5), data = X, scale = T)
# Df Variance      F Pr(>F)    
# Model      3    1.483 1.4948  0.001 ***
#   Residual 343  113.426                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Test of all variables
anova(modG, by="term") # significance by variable

# Tests of all canonical axes
anova(modG, by="axis") # significance by RDA

RsquareAdj(modG) 
# $r.squared
# [1] 0.01256704
# 
# $adj.r.squared
# [1] 0.004219837

# ORDISTEP model----
vif.cca(modG)
# dbMEM1            dbMEM2            dbMEM3            dbMEM4            dbMEM5 
# 30.778516          8.046730         16.163141          8.517340          3.450066 
# curM_Wclm_mean   sal_bottom_mean temp_bottom_range     pco_Wclm_mean 
# 1.669041         54.284101          7.645620                NA 
# Forward selection of explanatory variables to reduce the number of explanatory variables.
set.seed(5467); mod0 = rda(Y ~ 1, data = X, scale = T)  
RDA = ordistep(mod0, scope = formula(modG),perm.max=200, direction = "forward")
S=summary(RDA)

RDA$anova
anova(RDA)
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = Y ~ sal_bottom_mean + pco_Wclm_mean + temp_bottom_range, data = X, scale = T)
# Df Variance      F Pr(>F)    
# Model      3     2.16 2.1628  0.001 ***
#   Residual 348   115.84                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(RDA, by="term")
# Permutation test for rda under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = Y ~ sal_bottom_mean + pco_Wclm_mean + temp_bottom_range, data = X, scale = T)
# Df Variance      F Pr(>F)    
# sal_bottom_mean     1    0.891 2.6776  0.001 ***
#   pco_Wclm_mean       1    0.692 2.0788  0.001 ***
#   temp_bottom_range   1    0.577 1.7322  0.001 ***
#   Residual          348  115.840                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(RDA, by="axis")

RsquareAdj(RDA) 
# $r.squared
# [1] 0.01830387
# 
# $adj.r.squared
# [1] 0.009840973

# PLOT WITH VARIABLES IDENTIFIED BY ORDISTEP
#ordistep----
pdf("partial-RDA-ordistep-neutral-4vifcor-env.pdf")
par(mar = c(5, 5, 1, 7) + 0.5)
plot(RDA, xlab=paste0('RDA1 (',round(S$concont$importance[2,1]*RsquareAdj(RDA)[[2]]*100,2),'%)'), ylab=paste0('RDA2 (',round(S$concont$importance[2,2]*RsquareAdj(RDA)[[2]]*100,2),'%)'), type="none", cex.lab=1.2, cex.axis=1.2, scaling=2, display=c("lc","wa"))
points(RDA, display="wa", cex=1.6, scaling=2, col = "gray32", bg = bg[as.numeric(eco$STRATA)],pch=21)
# points(RDA, display="lc", pch=20, cex=3, col="red", scaling=2) CENTROID 
text(RDA, scaling=2, display="bp", col="darkcyan", cex=1.2, arrow.mul = 1.2, head.arrow = 0.1, adj=1, font=2)
text(3.3, -2.1, expression("R"["adj"]^2*"=0.0098"), cex=1.2)
text(3.35,-2.25,"p-value=0.001", cex=1.2)# the predictors
legend(3,1.7, legend=leg.text, bty="n", col="gray32", pch=21, cex=1.3, pt.bg=bg,y.intersp = 0.12, x.intersp = 0.12, xpd=T)
dev.off()

## Compute the Redundancy Analysis----
# GLOBAL model----
# RDA with all explanatory variables 
set.seed(231); modG = rda(Y ~ dbMEM1 + dbMEM2 + dbMEM3 + dbMEM4 + dbMEM5 + Condition(curM_Wclm_mean + sal_bottom_mean + temp_bottom_range + pco_Wclm_mean), data = X, scale = T)
G = summary(modG)

# Global test of the RDA result
anova(modG) # global significance
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = Y ~ dbMEM1 + dbMEM2 + dbMEM3 + dbMEM4 + dbMEM5 + Condition(curM_Wclm_mean + sal_bottom_mean + temp_bottom_range + pco_Wclm_mean), data = X, scale = T)
# Df Variance      F Pr(>F)    
# Model      4    2.024 1.5303  0.001 ***
#   Residual 343  113.426                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Test of all variables
anova(modG, by="term") # significance by variable

# Tests of all canonical axes
anova(modG, by="axis") # significance by RDA

RsquareAdj(modG) 
# $r.squared
# [1] 0.01715419
# 
# $adj.r.squared
# [1] 0.00601296

# ORDISTEP model----
# variance inflation factors in the RDA
vif.cca(modG)
# dbMEM1            dbMEM2            dbMEM3            dbMEM4            dbMEM5 
# 30.778516          8.046730         16.163141          8.517340          3.450066 
# curM_Wclm_mean   sal_bottom_mean temp_bottom_range     pco_Wclm_mean 
# 1.669041         54.284101          7.645620                NA 
# Forward selection of explanatory variables to reduce the number of explanatory variables.
set.seed(5467); mod0 = rda(Y ~ 1, data = X, scale = T)  
RDA = ordistep(mod0, scope = formula(modG),perm.max=200, direction = "forward")
S=summary(RDA)

RDA$anova
anova(RDA)
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = Y ~ dbMEM1 + dbMEM2 + dbMEM3 + dbMEM4, data = X, scale = T)
# Df Variance      F Pr(>F)    
# Model      4    2.708 2.0375  0.001 ***
#   Residual 347  115.292                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(RDA, by="term")
# Permutation test for rda under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = Y ~ dbMEM1 + dbMEM2 + dbMEM3 + dbMEM4, data = X, scale = T)
# Df Variance      F Pr(>F)    
# dbMEM1     1    0.880 2.6473  0.001 ***
#   dbMEM2     1    0.701 2.1099  0.001 ***
#   dbMEM3     1    0.584 1.7580  0.001 ***
#   dbMEM4     1    0.543 1.6350  0.001 ***
#   Residual 347  115.292                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(RDA, by="axis")

RsquareAdj(RDA) 
# $r.squared
# [1] 0.0229484
# 
# $adj.r.squared
# [1] 0.01168555

# PLOT WITH VARIABLES IDENTIFIED BY ORDISTEP
#ordistep----
pdf("partial-RDA-ordistep-neutral-4vifcor-dbmem.pdf")
par(mar = c(5, 5, 1, 7) + 0.5)
plot(RDA, xlab=paste0('RDA1 (',round(S$concont$importance[2,1]*RsquareAdj(RDA)[[2]]*100,2),'%)'), ylab=paste0('RDA2 (',round(S$concont$importance[2,2]*RsquareAdj(RDA)[[2]]*100,2),'%)'), type="none", cex.lab=1.2, cex.axis=1.2, scaling=2, display=c("lc","wa"))
points(RDA, display="wa", cex=1.6, scaling=2, col = "gray32", bg = bg[as.numeric(eco$STRATA)],pch=21)
# points(RDA, display="lc", pch=20, cex=3, col="red", scaling=2) CENTROID 
text(RDA, scaling=2, display="bp", col="darkcyan", cex=1.2, arrow.mul = 1.5, head.arrow = 0.1, adj=1, font=2)
text(3.6, -2.08, expression("R"["adj"]^2*"=0.012"), cex=1.2)
text(3.65,-2.25,"p-value=0.001", cex=1.2)# the predictors
legend(3.5,1.8, legend=leg.text, bty="n", col="gray32", pch=21, cex=1.3, pt.bg=bg,y.intersp = 0.15, x.intersp = 0.15, xpd=T)
dev.off()


