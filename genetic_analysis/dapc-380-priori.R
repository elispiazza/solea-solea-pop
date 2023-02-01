load("~/Desktop/manuscript-solea/graphs/dapc-380-priori.RData")
source("~/Documents/LAVORO_BETTA/my_data/adegenet_dapc-plots.R")
library(foreach)
library(tidyverse)
library(ggthemes)

# set.seed(110)
# dapc.out2.raw <- dapc(my_o, pop(my_o), 
#                       n.pca = NULL, n.da = NULL, 
#                       scale = F, var.contrib = T, 
#                       pca.info = T) 
# set.seed(110)
# ascore.out <- a.score(dapc.out2.raw, n.sim=10)
# set.seed(110)
# opt.ascore.out <- optim.a.score(dapc.out2.raw, 
#                                 n.pca=1:ncol(dapc.out2.raw$tab),
#                                 smart = F, plot=TRUE,
#                                 n.sim=10, 
#                                 n.da= 5)

#ALPHA-SCORE GRAPH----
boxplot(opt.ascore.out$pop.score, 
        col = "gold")
points(opt.ascore.out$mean, pch = "o")
points(opt.ascore.out$best, opt.ascore.out$mean[21],
       col="black", bg = "red",
       pch= 21,   
       cex=1)

boxplot(opt.ascore.out$pop.score[1:30], 
        col = "gold")
lines(opt.ascore.out$mean[1:30], lwd = 2)
points(opt.ascore.out$best, opt.ascore.out$mean[21],
       col="black", bg= "red",
       pch= 21,   
       cex=1.5)

# 
# my_mat <- tab(my_o, NA.method="mean")
# set.seed(110)
# xval.out1 <- xvalDapc(my_mat, pop(my_o), 
#                  training.set = 0.9,
#                  result = "groupMean", 
#                  center = TRUE, scale = FALSE,
#                  n.pca = NULL, n.rep = 30, 
#                  xval.plot = TRUE)
# set.seed(110)
# xval.out2 <- xvalDapc(my_mat, pop(my_o), 
#          n.pca = 100:200,
#          training.set = 0.9, result = c("groupMean"),
#          center = TRUE, scale = FALSE,
#          n.rep = 1000, xval.plot = TRUE)
# set.seed(110)
# xval.out3 <- xvalDapc(my_mat, pop(my_o), 
#                       n.pca = 100:200,
#                       training.set = 0.9, 
#                       result = c("groupMean"),
#                       center = TRUE, scale = FALSE,
#                       n.rep = 100, xval.plot = TRUE)


#PARAM SCATTER----
# new colors insted of palette above
strong <- c("#11A4C8","#63C2C5","#1D4F9F", "#2A2771",
            "#396D35","#80C342", "#808080",
            "#ED2224", "#ED3995","#7E277C",
            "#F8941E","#F7EC16","#8C2A1C")

my.col <- strong
my.lab <- c("GSA6_2000",
                "GSA6_1999",
                "GSA7_2009",
                "GSA7_2003",
                "GSA9_2009",
                "GSA10_2000",
                "GSA10_2003",
                "GSA17_2009",
                "GSA18w_2000",
                "GSA18e_2000",
                "GSA22_2009",
                "GSA24_2009",
                "GSA24_2002"
        )



#SCATTER 380 PRIORI---- 
set.seed(110)
dapc.out2 <- dapc(my_o, pop(my_o), n.pca = 125,
                  n.da = NULL, scale = F, var.contrib = T,
                  var.loadings = T, pca.info = T) 
dapc2.ld <- as.data.frame.matrix(dapc.out2$ind.coord)
dapc2.ld$cluster <- factor(dapc.out2$grp, levels = unique(dapc.out2$grp))

p2.1 <- scat_plot(dapc2.ld, my.col, LD1, LD2) +
        scale_fill_manual(values = my.col, labels = my.lab) +
        theme(plot.background = element_rect(color = 'white'),
                legend.title = element_blank(),
              text = element_text(size = 18),
              legend.key.size = unit(0.85, "cm")
              # legend.text = element_text(size = 16),
              # axis.title = element_text(size = 16),
              # axis.text = element_text(size = 14)
        )

ggsave(filename = "dapc-380-priori.pdf", plot = p2.1, dpi = 600, width = 269, height = 219.119, units = "mm", device = "pdf")

# > dapc.out2$eig[1]/sum(dapc.out2$eig) + dapc.out2$eig[2]/sum(dapc.out2$eig)
# [1] 0.8410595
# > dapc.out2$var
# [1] 0.8351989

#CLUSTER CURRENT VS INFERRED----
km.out_k <- table(dapc.out2$grp, dapc.out2$assign)
km.out_b <- as.data.frame.matrix(km.out_k)
rownames(km.out_b) <- my.lab
colnames(km.out_b) <- my.lab
km.out_b$Current <- factor(rownames(km.out_b), levels=rownames(km.out_b))
km.out_b$Current <- ordered(km.out_b$Current, levels=rev(my.lab))
km.out_bl <- as.data.frame(pivot_longer(km.out_b, -Current, names_to = "Inferred", values_to = "Size"))
km.out_bl$Inferred <- ordered(km.out_bl$Inferred, levels=my.lab)

addCol <- function(d){
        newc <- vector()
        foreach(n=1:nrow(d)) %do% if (as.character(d[,1][n]) == as.character(d[,2][n])){
                newc <- c(newc, "grey")
        } else {
                newc <- c(newc, "white")
        }
        return(newc)
}
km.out_bl$Color <- addCol(km.out_bl)

p1.1 <- km.out_bl %>%  filter(Size>0) %>%
        ggplot(aes(Inferred, Current, size=Size, fill=Color)) +
        geom_point(shape = 21) +
        scale_size("size", range=c(0, 15), breaks=c(10, 20, 30, 40)) +
        geom_text(aes(label = Size), size=3.2) +
        scale_fill_manual(values=c("grey", "white")) +
        theme_base() +
        ylab("Current population") + xlab("Inferred population") +
        theme(plot.background = element_rect(color = 'white'),
              legend.title = element_blank(),
              legend.position = "none",
              axis.text.x = element_text(angle=90),
              text = element_text(size = 16),
              panel.grid.major = element_line(size = 0.25, linetype = 'dashed', colour = "grey"),
              panel.grid.minor = element_line(size = 0.25, linetype = 'dashed', colour = "grey")
              # legend.text = element_text(size = 16),
              # axis.title = element_text(size = 16),
              # axis.text = element_text(size = 14)
              ) 
        #guides(color=guide_legend(), size=guide_legend()) 

ggsave(filename = "dapc-380-priori-cluster-size.pdf", plot = p1.1, dpi = 600, width = 259, height = 219, units = "mm", device = "pdf")



#DAPC 380 (COLOUR) + FINDCLUSTER (SHAPE)----
dapc1.ld <- read.csv("./dapc-380-findcluster.csv", header = T)
dapc12.ld <- cbind(dapc2.ld, as.factor(dapc1.ld$cluster))
colnames(dapc12.ld)[7] <- "clusterfind"


scat_plot2 <- function(dat, clrs, labs, x_ld, y_ld){
        x_ld <- enquo(x_ld)
        y_ld <- enquo(y_ld)
        p <- ggplot(dat, aes(!! x_ld, !! y_ld, group=cluster))
        p <- p + geom_point(aes(fill=cluster, shape=clusterfind),
                            size = 4.5,
                            alpha = .8)
        #p <- p + scale_fill_manual(values = clrs, labels=labs)
        #p <- p + scale_shape_manual(values=c(21, 24, 22))
        p <- p + theme_base()
}


dapc12.plot <- scat_plot2(dapc12.ld, my.col, my.lab, LD1, LD2) +
        scale_shape_manual(values=c(21, 24, 22),
                           labels = c("Cluster 1",
                                      "Cluster 2",
                                      "Cluster 3")
        ) +
        scale_fill_manual(values = my.col, labels = my.lab) +
        theme(plot.background = element_rect(color = 'white'),
              legend.title = element_blank(),
              text = element_text(size = 18),
              legend.key.size = unit(0.85, "cm")
              # legend.text = element_text(size = 16),
              # axis.title = element_text(size = 16),
              # axis.text = element_text(size = 14)
        ) +
        guides(fill = guide_legend(override.aes = list(shape = 21)))

ggsave(filename = "dapc-380-findcluster&priori.pdf", plot = dapc12.plot, dpi = 600, width = 269, height = 219.119, units = "mm", device = "pdf")




# 
dapc.out2_info <- summary(dapc.out2)
assignplot(dapc.out2, subset=1:25)

compoplot(dapc.out2, 
          txt.leg=paste("Cluster", 1:13), lab="",
          ncol=1, xlab="individuals")
dapc2.post.prob <- as.data.frame.matrix(dapc.out2$posterior)
dapc2.post.prob$ind <- rownames(dapc2.post.prob)
dapc2.post.prob$ind <- factor(dapc2.post.prob$ind, levels=unique(dapc2.post.prob$ind))
dapc2.post.prob$K <- 13
dapc2.post.prob$geo_sample <- factor(pop(my_o), levels = unique(pop(my_o))) 
dapc2.post.prob <- pivot_longer(dapc2.post.prob, -c(K, ind, geo_sample), names_to = "cluster", values_to = "prob")
dapc2.post.prob$cluster <- factor(dapc2.post.prob$cluster, levels = unique(dapc2.post.prob$cluster)) 
dapc2.post.prob <- arrange(dapc2.post.prob, cluster)  
p2.5 <- postprob_plot(dapc2.post.prob, dapc2.col) +
        theme(legend.position = "bottom") + guides(col=guide_legend(nrow = 2))


p2.5 <- postprob_plot(dapc2.post.prob, my.col, my.lab) +
        theme(plot.background = element_rect(color='white'))

ggsave(filename = "dapc2-posterior-prob.pdf", plot = p2.5, path = "./my_data/report/", dpi = 600, width =21, height= 15, units = "cm", device = "pdf")




