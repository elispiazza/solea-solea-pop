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
?chooseCN #this will show your options for the SPCA type

# fuction to group population in gsa 
poolSome <- function(l, p1, p2, np){
  lI <- l
  pos1 <- p1
  pos2 <- p2
  newnames <- np
  lF <- c()
  i = 1
  k = 1
  j = 1
  while (i < length(lI)){
    if (i %in% pos1){
      a <- pos1[k]
      b <- pos2[k]
      new.pop <- repool(lI[[a]], lI[[b]])
      pop(new.pop) <- factor(rep(newnames[j], nInd(new.pop)), levels = newnames[j])
      indNames(new.pop) <- paste0(newnames[j], paste0("_",1:length(indNames(new.pop))))
      lF <- c(lF, new.pop)
      i = i+2
      k = k+1
      j = j+1
    } else {
      new.pop <- lI[[i]]
      pop(new.pop) <- factor(rep(newnames[j], nInd(new.pop)), levels = newnames[j])
      indNames(new.pop) <- paste0(newnames[j], paste0("_",1:length(indNames(new.pop))))
      lF <- c(lF, new.pop)
      i = i+1
      j = j+1
    }
  }
  names(lF) <- newnames
  newgenind <- repool(lF)
  return(newgenind)
}


## Load coordinetes ------------

coord.gsa <- read.csv("gsa-coord-cart.csv",header=T)
coord <- read.csv("env-data-paper-new-order.csv", header = T)
coord$ID_indiv <- indNames(genindGSA_neut)
coord.ind <- coord[,c(1:3)]
coord.ind$MDS1 <- rep(coord.gsa$MDS1, table(pop(genindGSA_neut)))
coord.ind$MDS2 <- rep(coord.gsa$MDS2, table(pop(genindGSA_neut)))
colnames(coord.ind)[1:3] <- c("Code", "Long", "Lat")


##Jitter the species coordieurotes so that a Delaueuroy triangulation can be applied. Here we use a jitter of 2 km 
coord.ind$jCartx <- jitter(coord.ind$MDS1,1) 
coord.ind$jCarty <- jitter(coord.ind$MDS2,1)

write.csv(coord.ind, "gsa-coord-ind-cart-jitter.csv", quote = F, row.names = F)

## Read in genepop files as GENIND objects------------
# neutral 
genindGSA_neut <- read.genepop("./my_genepop-GSA-neutral.gen")

# outlier 
genindGSA_out <- read.genepop("./my_genepop-GSA-outlier.gen")

# all
genindGSA_all <- read.genepop("./my_genepop-GSA-all.gen")

# to evaluate check Inds names and 
inds<-genepop_detective("my_genepop-GSA-neutral.gen",variable = "Inds")
duplicated(inds)
sum(duplicated(inds))

##Add the spatial coordinates -----

#######Neutral loci#####

neutral_Cart <- genindGSA_neut # cartesian coordieurotes
neutral_Cart@other$xy <- coord.ind[,c("MDS1","MDS2")]
head(neutral_Cart@other$xy)
tail(neutral_Cart@other$xy)

neutral_jCart <- genindGSA_neut # jittered Cartesian coordieurotes (Delaueuroy triangulation)
neutral_jCart@other$xy <- coord.ind[,c("jCartx","jCarty")]
head(neutral_jCart@other$xy)
tail(neutral_jCart@other$xy)

neutral_Coords <- genindGSA_neut # Standard lat long
neutral_Coords@other$xy <- coord.ind[,c("Long","Lat")]


#######Outlier loci#####
outlier_Cart <- genindGSA_out # cartesian coordieurotes
outlier_Cart@other$xy <- coord.ind[,c("MDS1","MDS2")]

outlier_jCart <- genindGSA_out # jittered Cartesian coordieurotes (Delaueuroy triangulation)
outlier_jCart@other$xy <- coord.ind[,c("jCartx","jCarty")]

outlier_Coords <- genindGSA_out # Standard lat long
outlier_Coords@other$xy <- coord.ind[,c("Long","Lat")]

save.image("SPCA_Data.RData")
gc()

#######All loci#####
all_Cart <- genindGSA_all # cartesian coordieurotes
all_Cart@other$xy <- coord.ind[,c("MDS1","MDS2")]

all_jCart <- genindGSA_all # jittered Cartesian coordieurotes (Delaueuroy triangulation)
all_jCart@other$xy <- coord.ind[,c("jCartx","jCarty")]

all_Coords <- genindGSA_all # Standard lat long
all_Coords@other$xy <- coord.ind[,c("Long","Lat")]

#########Step 3. Now run the SPCA########### ------
##### Neutral loci -----
#### Delaunay (type = 1)
neutral_sPCA_jCart <- spca(neutral_jCart,ask=FALSE,type=1,scannf = FALSE)
barplot(neutral_sPCA_jCart$eig,col=spectral(length(neutral_sPCA_jCart$eig)))
plot(neutral_sPCA_jCart)
#legend("topright", fill=spectral(2),
#      leg=c("Global structures", "Local structures"))
#abline(h=0,col="grey")

save.image("SPCA_Data.RData")

png("neutral_SPCA_Delaunay_interpolatedMap.png",width=2500,height=2000,res=300)
x <- other(neutral_Cart)$xy[,1]
y <- other(neutral_Cart)$xy[,2]
interpX <- seq(min(x),max(x),le=200)
interpY <- seq(min(y),max(y),le=200)
temp <- interp(x, y, neutral_sPCA_jCart$ls[,1], xo=interpX, yo=interpY, duplicate = "mean")
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

#Global and local tests for significance on neutral 
neutral_jCart_naomit <- tab(neutral_jCart, NA.method = "mean")
neutralGtest<-global.rtest(neutral_jCart_naomit,neutral_sPCA_jCart$lw,nperm=9999)
plot(neutralGtest)
neutralLtest<-local.rtest(neutral_jCart_naomit,neutral_sPCA_jCart$lw,nperm=9999)
plot(neutralLtest)

#Save data again
save.image("SPCA_Data.RData")

#SPCA axis 1 lagged scores for neutral plotted against latitude for the sampling site demonstrating a longitudinal genetic cline

plot(neutral_sPCA_jCart$ls[,1], coord.ind$Long)
lines(neutral_sPCA_jCart$ls[,1], coord.ind$Long)

A <- data.frame(Code = pop(genindGSA_neut), Long = coord.ind[,2], Scores =neutral_sPCA_jCart$ls[,1])
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
  geom_point(data = B, aes(fill = Code), colour = "grey32", size = 6, shape = 21) +
  scale_color_manual(values = bg)+
  scale_fill_manual(values = bg)+
  xlab("sPCA axis 1 lagged scores") + 
  ylab("Longitude °E")+
  scale_y_reverse() +
  scale_x_reverse()+
  theme_base()

# Which alleles exhibit structure (squared loading)
myLoadings <- neutral_sPCA_jCart$c1[,1]^2
names(myLoadings) <- rownames(neutral_sPCA_jCart$c1)
loadingplot(myLoadings, xlab="Alleles",
            ylab="Weight of the alleles",
            main="Contribution of alleles \n to the first sPCA axis")
# AVERAGE CONTRIBUTION OF EACH MARKER, BOXPLOT
boxplot(myLoadings~neutral_jCart$loc.fac, las=3, ylab="Contribution", xlab="Marker",
        main="Contributions by markers \nto the first global score", col="grey")
# threshold value above which alleles are annotated.
#to identify the 5% of alleles with 
# the greatest contributions to the first global structure in Spca result
temp_neut <- loadingplot(myLoadings, threshold=quantile(myLoadings, 0.95),
                    xlab="Alleles",ylab="Weight of the alleles",
                    main="Contribution of alleles \n to the first sPCA axis",
                    cex.fac=0.6, srt = 45, lab.jitter = 60, adj = c(0.7, 0.3))


##### Outlier loci -----
##Delaunay (type=1)   
outlier_sPCA_jCart <- spca(outlier_jCart,ask=FALSE,type=1,scannf = FALSE)
barplot(outlier_sPCA_jCart$eig,col=spectral(length(outlier_sPCA_jCart$eig)))
legend("topright", fill=spectral(2),
       leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")

png("outlier_SPCA_Delaunay_interpolatedMap.png",width=2500,height=2000,res=300)
x <- other(outlier_Cart)$xy[,1]
y <- other(outlier_Cart)$xy[,2]
interpX <- seq(min(x),max(x),le=200)
interpY <- seq(min(y),max(y),le=200)
temp2 <- interp(x, y, outlier_sPCA_jCart$ls[,1], xo=interpX, yo=interpY, duplicate = "mean")
myPal <- colorRampPalette(c("#8C2A1C", "#F7EC16", 
                            "azure", "#11A4C8", "#2A2771"))
annot <- function(){
  text(unique(x), unique(y), labels=gsa, pos=c(2,1,3,1,2,3,2,1,1))
  #points(x,y)
  title(xlab = "Cartesian 1", ylab = "Cartesian 2")
}
filled.contour(temp2, color.pal=myPal, nlev=50, xlim=c(-1900, 1900), ylim=c(-800, 1100),
               key.title=title("lagged \nscore 1"), plot.title=annot())
dev.off()

#Global and local tests for significance on outliers
outlier_jCart_naomit <- tab(outlier_jCart, NA.method = "mean")
outGtest<-global.rtest(outlier_jCart_naomit,outlier_sPCA_jCart$lw,nperm=9999)
plot(outGtest)
outLtest<-local.rtest(outlier_jCart_naomit,outlier_sPCA_jCart$lw,nperm=9999)
plot(outLtest)

#SPCA axis 1 lagged scores for outlier plotted against latitude for the sampling site demonstrating a longitudinal genetic cline

plot(outlier_sPCA_jCart$ls[,1], coord.ind$Long)
lines(outlier_sPCA_jCart$ls[,1], coord.ind$Long)

C <- data.frame(Code = pop(genindGSA_out), Long = coord.ind[,2], Scores =outlier_sPCA_jCart$ls[,1])
D <- C %>% group_by(Code) %>% mutate(meanLagS = mean(Scores))
D <- as.data.frame(D)

ggplot(D, aes(meanLagS, Long)) + 
  geom_smooth(data = D, method = "loess", se = F, col= "black", orientation = "x", span = 0.7, size = 0.5) + 
  geom_text(data = D, aes(label = Code), size = 5, position = position_nudge(y = 1.7))+
  geom_point(data = D, aes(fill = Code), colour = "grey32", size = 6, shape = 21) +
  scale_color_manual(values = bg)+
  scale_fill_manual(values = bg)+
  xlab("sPCA axis 1 lagged scores") + 
  ylab("Longitude °E")+
  scale_y_reverse() +
  scale_x_reverse()+
  theme_base()

# Which alleles exhibit structure (squared loading)
myLoadings <- outlier_sPCA_jCart$c1[,1]^2
names(myLoadings) <- rownames(outlier_sPCA_jCart$c1)
loadingplot(myLoadings, xlab="Alleles",
            ylab="Weight of the alleles",
            main="Contribution of alleles \n to the first sPCA axis")
# AVERAGE CONTRIBUTION OF EACH MARKER, BOXPLOT
boxplot(myLoadings~outlier_jCart$loc.fac, las=3, ylab="Contribution", xlab="Marker",
        main="Contributions by markers \nto the first global score", col="grey")
# threshold value above which alleles are annotated.
#to identify the 5% of alleles with 
# the greatest contributions to the first global structure in Spca result
temp_out <- loadingplot(myLoadings, threshold=quantile(myLoadings, 0.95),
                    xlab="Alleles",ylab="Weight of the alleles",
                    main="Contribution of alleles \n to the first sPCA axis",
                    cex.fac=0.6, srt = 45, lab.jitter = 35, adj = c(0.7, 0.3))

### All loci ####
#### Delaunay (type = 1)
all_sPCA_jCart <- spca(all_jCart,ask=FALSE,type=1,scannf = FALSE)
barplot(all_sPCA_jCart$eig,col=spectral(length(all_sPCA_jCart$eig)))
plot(all_sPCA_jCart)
#legend("topright", fill=spectral(2),
#      leg=c("Global structures", "Local structures"))
#abline(h=0,col="grey")

save.image("SPCA_Data.RData")

png("all_SPCA_Delaunay_interpolatedMap.png",width=2500,height=2000,res=300)
x <- other(all_Cart)$xy[,1]
y <- other(all_Cart)$xy[,2]
interpX <- seq(min(x),max(x),le=200)
interpY <- seq(min(y),max(y),le=200)
temp4 <- interp(x, y, all_sPCA_jCart$ls[,1], xo=interpX, yo=interpY, duplicate = "mean")
myPal <- colorRampPalette(c("#8C2A1C", "#F7EC16", 
                            "azure", "#11A4C8", "#2A2771"))
annot <- function(){
  text(unique(x), unique(y), labels=gsa, pos=c(2,1,3,1,2,3,2,1,1))
  #points(x,y)
  title(xlab = "Cartesian 1", ylab = "Cartesian 2")
  #axis(side=2, at=seq(-800, 1000, by=200))
}
filled.contour(temp4, color.pal=myPal, nlev=50, xlim=c(-1900, 1900), ylim=c(-800, 1100),
               key.title=title("lagged \nscore 1"), plot.title=annot())
dev.off()

#Global and local tests for significance on neutral 

all_jCart_naomit <- tab(all_jCart, NA.method = "mean")
allGtest<-global.rtest(all_jCart_naomit,all_sPCA_jCart$lw,nperm=9999)
plot(allGtest)
allLtest<-local.rtest(all_jCart_naomit,all_sPCA_jCart$lw,nperm=9999)
plot(allLtest)

#Save data again
save.image("SPCA_Data.RData")

#SPCA axis 1 lagged scores for neutral plotted against latitude for the sampling site demonstrating a longitudinal genetic cline

plot(all_sPCA_jCart$ls[,1], coord.ind$Long)
lines(all_sPCA_jCart$ls[,1], coord.ind$Long)

E <- data.frame(Code = pop(genindGSA_all), Long = coord.ind[,2], Scores =all_sPCA_jCart$ls[,1])
G <- E %>% group_by(Code) %>% mutate(meanLagS = mean(Scores))
G <- as.data.frame(G)
bg <- c("#11A4C8", "#2A2771",
        "#396D35", "#808080",
        "#ED2224", "#ED3995","#7E277C",
        "#F8941E","#F7EC16")

ggplot(G, aes(meanLagS, Long)) + 
  geom_smooth(data = G, method = "loess", se = F, col= "black", orientation = "x", span = 0.5, size = 0.5) + 
  geom_text(data = G, aes(label = Code), size = 5, position = position_nudge(y = 1.7))+
  geom_point(data = G, aes(fill = Code), colour = "grey32", size = 6, shape = 21) +
  scale_color_manual(values = bg)+
  scale_fill_manual(values = bg)+
  xlab("sPCA axis 1 lagged scores") + 
  ylab("Longitude °E")+
  scale_y_reverse() +
  scale_x_reverse()+
  theme_base()


# Which alleles exhibit structure (squared loading)
myLoadings <- all_sPCA_jCart$c1[,1]^2
names(myLoadings) <- rownames(all_sPCA_jCart$c1)
loadingplot(myLoadings, xlab="Alleles",
            ylab="Weight of the alleles",
            main="Contribution of alleles \n to the first sPCA axis")
# AVERAGE CONTRIBUTION OF EACH MARKER, BOXPLOT
boxplot(myLoadings~all_jCart$loc.fac, las=3, ylab="Contribution", xlab="Marker",
        main="Contributions by markers \nto the first global score", col="grey")
# threshold value above which alleles are annotated.
#to identify the 5% of alleles with 
# the greatest contributions to the first global structure in Spca result
temp_all <- loadingplot(myLoadings, threshold=quantile(myLoadings, 0.95),
                    xlab="Alleles",ylab="Weight of the alleles",
                    main="Contribution of alleles \n to the first sPCA axis",
                    cex.fac=0.6, srt = 45, lab.jitter = 35, adj = c(0.7, 0.3))

unique(gsub("\\.[0-9]*", "", temp_all$var.names))[unique(gsub("\\.[0-9]*", "", temp_all$var.names)) %in% locNames(genindGSA_out)]
