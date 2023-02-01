library(adegenet)
library(BiodiversityR)
library(ggplot2)
library(reshape2)
library(vegan)

#f fourier
#s shape
#b both


#plot param----
strong <- c("#11A4C8","#63C2C5","#1D4F9F", "#2A2771",
            "#396D35","#80C342", "#808080",
            "#ED2224", "#ED3995","#7E277C",
            "#F8941E","#F7EC16","#8C2A1C")
my.col <- strong[c(3,5,8,11,12)]
my.lab <- c("GSA7_2009","GSA9_2009",
            "GSA17_2009",
            "GSA22_2009", "GSA24_2009")




##BASE----
## fourier----
f.f <- read.csv("~/Documents/LAVORO_BETTA/my_data/otolith-analysis/fourier-momocs-13.csv", header = T, na.strings = "")

## Prepare the response and explanatory variables
# Calculate Euclidean distances to estimate the response variable
# extract matrix of allele frequency
Y.f <- f.f[,4:ncol(f.f)]
Y.f <- dist(Y.f, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
# explanatory
X.f <- f.f[,2:3]


## CAP
# plot change in classification success against m 
# plot(seq(1:100), rep(-1000, 100), xlim=c(1, 100), ylim=c(0, 100), xlab="m", 
#      ylab="Classification success (percent)", type="n")
# perc.max <- vector()
# for (mseq in 1:100) {
#   res.tmp <- CAPdiscrim(Y~population.code, data=X, axes=2, m=mseq)
#   perc.max <- c(perc.max, res.tmp$percent)
#   points(mseq, res.tmp$percent)
# }
# which.max(perc.max)

#13
cap.f <- CAPdiscrim(Y.f~population.code, X.f, m=30)
# Overall classification success (m=30) : 69.3121693121693 percent
# A (n=47) correct: 95.7446808510638 percent
# GL (n=30) correct: 60 percent
# NA (n=44) correct: 65.9090909090909 percent
# NT (n=39) correct: 48.7179487179487 percent
# TU (n=29) correct: 68.9655172413793 percent
cap.f$varm


#plot oridination
X.f$population.code <- ordered(X.f$population.code,
                             levels = c("GL", "NT", "NA", "A", "TU"))


# res <- capscale(Y ~ population.code, data=X, na.action = na.omit)
# s <- summary(res)
# s$concont$importance[2,1]
# s$concont$importance[2,2]




#plot 13
pdf("./CAP-fourier13.pdf")
#par(mar = c(5, 5, 1, 7) + 0.1)
p.f <- ordiplot(cap.f,type = "none",
              xlab='CAP1', ylab='CAP2', cex.lab=1.2, cex.axis=1.2)
#orditorp(p,lwd=2.5, display="sites",col=my.col[as.numeric(X$population.code)], pch = c(15,16,17,18,19)[as.numeric(X$population.code)], labels = F, cex = 1.3)
orditorp(p.f,lwd=1.2, display="sites", col="black", bg=alpha(my.col[as.numeric(X.f$population.code)], .8), pch = 21, labels = F, cex = 1.3)
#Create 95% ellipse
pl.f <- ordiellipse(cap.f, groups = X.f$population.code, kind="se", conf=0.95, lwd=0.1, lty=0, draw = "polygon", label = F, col=my.col, alpha = 0.4)
#legend(5.5,-0.7, legend=levels(X$population.code), bty="n", col=my.col, pch=c(15,16,17,18,19), cex=1.3,y.intersp = 0.5, x.intersp = 0.5, xpd=T)
#legend(5.5,-0.7, legend=my.lab, bty="n", pch=21, col="black", pt.bg=unique(alpha(my.col[as.numeric(X.f$population.code)], .8)), cex=1.3,y.intersp = 0.5, x.intersp = 0.5, xpd=T)


#plot 13
# cap$group <- ordered(cap$group,
#                      levels = c("GL", "NT", "NA", "A", "TU"))
# cap$CV <- ordered(cap$CV,
#                   levels = c("GL", "NT", "NA", "A", "TU"))
# t <- table(cap$group, cap$CV)
# GL NT NA  A TU
# GL 18  4  6  1  1
# NT  4 19  6  7  3
# NA  4 11 29  0  0
# A   1  1  0 45  0
# TU  2  4  3  0 20

# grp <- as.data.frame.matrix(t)
# grp$actual <- rownames(t)
# grp$actual <- factor(grp$actual, levels = c("GL", "NT", "NA", "A", "TU"))
# grpL <- melt(grp, id.vars = "actual", variable.name = "inferred", value.name = "size")
# 
# 
# ggplot(grpL, aes(inferred, actual, size=size, col=actual)) +
#   geom_point() +
#   xlab("Inferred group") + ylab("") +
#   scale_color_manual(values = my.col) +
#   theme_base(base_size = 15) +
#   guides(size=guide_legend())



## shape----
f.s <- read.csv("~/Documents/LAVORO_BETTA/my_data/otolith-analysis/shape-index-corrected.csv", header = T, na.strings = "")

## Prepare the response and explanatory variables
# Calculate Euclidean distances to estimate the response variable
# extract matrix of allele frequency
Y.s <- f.s[,4:6]
Y.s <- dist(Y.s, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
# explanatory
X.s <- f.s[,1:2]


## CAP
# # plot change in classification success against m 
# plot(seq(1:100), rep(-1000, 100), xlim=c(1, 100), ylim=c(0, 100), xlab="m", 
#      ylab="Classification success (percent)", type="n")
# perc.max <- vector()
# for (mseq in 1:100) {
#   res.tmp <- CAPdiscrim(Y~population.code, data=X, axes=2, m=mseq)
#   perc.max <- c(perc.max, res.tmp$percent)
#   points(mseq, res.tmp$percent)
# }
# which.max(perc.max)

#cap
cap.s <- CAPdiscrim(Y.s~population.code, X.s, m=3)
# Overall classification success (m=3) : 37.9120879120879 percent
# A (n=45) correct: 40 percent
# GL (n=28) correct: 75 percent
# NA (n=42) correct: 9.52380952380952 percent
# NT (n=38) correct: 34.2105263157895 percent
# TU (n=29) correct: 44.8275862068966 percent
cap.s$varm

#plot oridination
X.s$population.code <- ordered(X.s$population.code,
                             levels = c("GL", "NT", "NA", "A", "TU"))
# res <- capscale(Y ~ population.code, data=X, na.action = na.omit)
# s <- summary(res)
# s$concont$importance[2,1]
# s$concont$importance[2,2]


pdf("./CAP-shape.pdf")
#par(mar = c(5, 5, 1, 7) + 0.1)
p.s <- ordiplot(cap.s,type = "none",
              xlab='CAP1', ylab='CAP2', cex.lab=1.2, cex.axis=1.2)
#orditorp(p,lwd=2.5, display="sites",col=my.col[as.numeric(X$population.code)], pch = c(15,16,17,18,19)[as.numeric(X$population.code)], labels = F, cex = 1.3)
orditorp(p.s,lwd=1.2, display="sites", col="black", bg=alpha(my.col[as.numeric(X.s$population.code)], .8), pch = 21, labels = F, cex = 1.3)
#Create 95% ellipse
pl.s <- ordiellipse(cap.s, groups = X.s$population.code, kind="se", conf=0.95, lwd=0.1, lty=0, draw = "lines", label = F, col=my.col, alpha = 0.4)
#legend(5.4,0.8, legend=levels(X$population.code), bty="n", col=my.col, pch=c(15,16,17,18,19), cex=1.3,y.intersp = 0.4, x.intersp = 0.5, xpd=T)
#legend(5.4,0.8, legend=my.lab, bty="n", col="black", pt.bg=unique(alpha(my.col[as.numeric(X.s$population.code)], .8)), pch=21, cex=1.3,y.intersp = 0.4, x.intersp = 0.4, xpd=T)
# Get Ordiellipse output
summary(pl)

#plot group assignments
# cap$group <- ordered(cap$group,
#                      levels = c("GL", "NT", "NA", "A", "TU"))
# cap$CV <- ordered(cap$CV,
#                   levels = c("GL", "NT", "NA", "A", "TU"))
# t <- table(cap$group, cap$CV)
# GL NT NA  A TU
# GL 21  0  0  7  0
# NT  2 13 11 11  1
# NA  1 12  4 19  6
# A   6 10  4 18  7
# TU  0  6  4  6 13

# grp <- as.data.frame.matrix(t)
# grp$actual <- rownames(t)
# grp$actual <- factor(grp$actual, levels = c("GL", "NT", "NA", "A", "TU"))
# grpL <- melt(grp, id.vars = "actual", variable.name = "inferred", value.name = "size")
# 
# 
# ggplot(grpL, aes(inferred, actual, size=size, col=actual)) +
#   geom_point() +
#   xlab("Inferred group") + ylab("") +
#   scale_color_manual(values = my.col) +
#   theme_base(base_size = 15) +
#   guides(size=guide_legend()) 


#check which individuals are misslassified



##both----
f.b <- read.csv("~/Documents/LAVORO_BETTA/my_data/otolith-analysis/shape-index-fourier-13-181ind-cap.csv", header = T, na.strings = "")
ncol(f.b)

## Prepare the response and explanatory variables
# Calculate Euclidean distances to estimate the response variable
# extract matrix of allele frequency
Y.b <- f.b[,3:ncol(f.b)]
Y.b <- dist(Y.b, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
# explanatory
X.b <- f.b[,1:2]


# plot change in classification success against m 
plot(seq(1:100), rep(-1000, 100), xlim=c(1, 100), ylim=c(0, 100), xlab="m",
     ylab="Classification success (percent)", type="n")
perc.max <- vector()
for (mseq in 1:100) {
  res.tmp <- CAPdiscrim(Y.b~population.code, data=X.b, axes=2, m=mseq)
  perc.max <- c(perc.max, res.tmp$percent)
  points(mseq, res.tmp$percent)
}
which.max(perc.max)

## CAP
cap.b <- CAPdiscrim(Y.b~population.code, X.b, axes=4, m=17)
# Overall classification success (m=17) : 72.9281767955801 percent
# A (n=45) correct: 93.3333333333333 percent
# GL (n=28) correct: 82.1428571428571 percent
# NA (n=42) correct: 73.8095238095238 percent
# NT (n=38) correct: 47.3684210526316 percent
# TU (n=28) correct: 64.2857142857143 percent
cap.b$varm

#plot oridination
X.b$population.code <- ordered(X.b$population.code,
                             levels = c("GL", "NT", "NA", "A", "TU"))
#plot param
# res <- capscale(Y ~ population.code, data=X, na.action = na.omit)
# s <- summary(res)
# s$concont$importance[2,1]
# s$concont$importance[2,2]

#plot 30
pdf("./CAP-shape-fourier30.pdf")
#par(mar = c(5, 5, 1, 7) + 0.1)
p.b <- ordiplot(cap.b,type = "none",
              xlab='CAP1', ylab='CAP2', cex.lab=1.2, cex.axis=1.2)
#orditorp(p,lwd=2.5, display="sites",col=my.col[as.numeric(X$population.code)], pch = c(15,16,17,18,19)[as.numeric(X$population.code)], labels = F, cex = 1.3)
orditorp(p.b,lwd=1.2, display="sites", col="black", bg=alpha(my.col[as.numeric(X.b$population.code)], .8), pch = 21, labels = F, cex = 1.3)
#Create 95% ellipse
pl.b <- ordiellipse(cap.b, groups = X.b$population.code, kind="se", conf=0.95, lwd=0.1, lty=0, draw = "polygon", label = F, col=my.col, alpha = 0.4)
#legend(6,0.3, legend=levels(X$population.code), bty="n", col=my.col, pch=c(15,16,17,18,19), cex=1.3,y.intersp = 0.5, x.intersp = 0.5, xpd=T)
legend(6,0.3, legend=my.lab, bty="n", col="black", pt.bg=unique(alpha(my.col[as.numeric(X.b$population.code)], .8)), pch=21, cex=1.3,y.intersp = 0.4, x.intersp = 0.4, xpd=T)
# Get Ordiellipse output
summary(pl)




#plot group assignment
#plot 30
cap.b$group <- ordered(cap.b$group,
                     levels = c("GL", "NT", "NA", "A", "TU"))
cap.b$CV <- ordered(cap.b$CV,
                  levels = c("GL", "NT", "NA", "A", "TU"))
t <- table(cap.b$group, cap.b$CV)
# GL NT NA  A TU
# GL 23  0  2  3  0
# NT  2 18 11  4  3
# NA  2  9 31  0  0
# A   0  2  0 42  1
# TU  0  3  5  2 18


# grp <- as.data.frame.matrix(t)
# grp$actual <- rownames(t)
# grp$actual <- factor(grp$actual, levels = c("GL", "NT", "NA", "A", "TU"))
# grpL <- melt(grp, id.vars = "actual", variable.name = "inferred", value.name = "size")
# 
# 
# ggplot(grpL, aes(inferred, actual, size=size, col=actual)) +
#   geom_point() +
#   xlab("Inferred group") + ylab("") +
#   scale_color_manual(values = my.col) +
#   theme_base(base_size = 15) +
#   guides(size=guide_legend()) 





#GGPLOT----
library(ggthemes)
library(ggforce)
source("./ordiellipse-long-function.R")
#' CAP scatterplot
#'
#' @param dat dataframe of coordinates of the individuals with columns LD1, LD2, ..., cluster
#' @param clrs color vector for clusters
#' @param labs labels vector for clusters
#' @param x_ld LD x axis 
#' @param y_ld LD y axis
#'
#' @return cap plot
#' @export
#'
#' @examples scat_plot(dapc.ld, dapc.col)
cap_plot <- function(dat, clrs, labs, x_ld, y_ld) {
  x_ld <- enquo(x_ld)
  y_ld <- enquo(y_ld)
  p <- ggplot(dat, aes(!!x_ld,!!y_ld, fill = pop))
  p <- p + geom_point(size = 3.5,
                      shape = 21,
                      alpha = .8)
  p <- p + scale_fill_manual(values = clrs, labels = labs)
  p <- p + xlab("CAP1") + ylab("CAP2")
  p <- p + theme_base()
}

d.f <- data.frame(pop= f.f$population.code, p.f$sites)
d.f$pop <- ordered(d.f$pop, levels = c("GL", "NT", "NA", "A", "TU"))
ellips.f <- ordiellipse.long(pl.f, grouping.name = "population.code")
ellips.f$population.code <- ordered(ellips.f$population.code, levels = c("GL", "NT", "NA", "A", "TU"))
fourier.p <- cap_plot(d.f, my.col, my.lab, LD1, LD2) +
  stat_ellipse(type = "t") +
  theme(plot.background = element_rect(color = 'white'),
        legend.title = element_blank(),
        text = element_text(size = 16)
        # legend.text = element_text(size = 16),
        # axis.title = element_text(size = 16),
        # axis.text = element_text(size = 14)
        )

d.s <- data.frame(pop= f.s$population.code, p.s$sites)
d.s$pop <- ordered(d.s$pop, levels = c("GL", "NT", "NA", "A", "TU"))
ellips.s <- ordiellipse.long(pl.s, grouping.name = "population.code")
shape.p <- cap_plot(d.s, my.col, LD1, LD2) +
  theme(plot.background = element_rect(color = 'white'),
        legend.title = element_blank())

d.b <- data.frame(pop= f.b$population.code, p.b$sites)
d.b$pop <- ordered(d.b$pop, levels = c("GL", "NT", "NA", "A", "TU"))
ellips.b <- ordiellipse.long(pl.b, grouping.name = "population.code")
both.p <- cap_plot(d.b, my.col, LD1, LD2) +
  theme(plot.background = element_rect(color = 'white'),
        legend.title = element_blank())


library(patchwork)
panel <- ((shape.p / fourier.p) |both.p) + plot_layout(guides = 'collect')


##GGPLOT-BIS----
library(ggordiplots)
#fourier----
fourier <- gg_ordiplot(
  cap.f,
  d.f$pop,
  kind = "se",
  conf = 0.95,
  show.groups = "all",
  ellipse = TRUE,
  label = FALSE,
  hull = FALSE,
  spiders = FALSE,
  pt.size = 5,
  plot = TRUE, 
) 

d.f2 <- fourier$df_ord
ellips.f2 <- fourier$df_ellipse

fourier.p2 <- ggplot(d.f2, aes(x,y, fill = Group)) + 
  geom_point(size = 5,
                    shape = 21,
                    alpha = .8) +
  geom_polygon(data = ellips.f2, aes(color=Group, fill= Group), alpha = .4, show.legend = FALSE) +
  scale_fill_manual(values = my.col, labels = my.lab) +
  scale_color_manual(values = my.col) +
  xlab("CAP1") + ylab("CAP2") +
  theme_base() +
  theme(plot.background = element_rect(color = 'white'),
        legend.title = element_blank(),
        text = element_text(size = 16)
        # legend.text = element_text(size = 16),
        # axis.title = element_text(size = 16),
        # axis.text = element_text(size = 14)
  )


#shape----
shape <- gg_ordiplot(
  cap.s,
  d.s$pop,
  kind = "se",
  conf = 0.95,
  show.groups = "all",
  ellipse = TRUE,
  label = FALSE,
  hull = FALSE,
  spiders = FALSE,
  pt.size = 5,
  plot = TRUE, 
) 

d.s2 <- shape$df_ord
ellips.s2 <- shape$df_ellipse

shape.p2 <- ggplot(d.s2, aes(x,y, fill = Group)) + 
  geom_point(size = 5,
             shape = 21,
             alpha = .8) +
  geom_polygon(data = ellips.s2, aes(color=Group, fill= Group), alpha = .4, show.legend = FALSE) +
  scale_fill_manual(values = my.col, labels = my.lab) +
  scale_color_manual(values = my.col) +
  xlab("CAP1") + ylab("CAP2") +
  theme_base() +
  theme(plot.background = element_rect(color = 'white'),
        legend.title = element_blank(),
        text = element_text(size = 16)
        # legend.text = element_text(size = 16),
        # axis.title = element_text(size = 16),
        # axis.text = element_text(size = 14)
  )


#both----
both <- gg_ordiplot(
  cap.b,
  d.b$pop,
  kind = "se",
  conf = 0.95,
  show.groups = "all",
  ellipse = TRUE,
  label = FALSE,
  hull = FALSE,
  spiders = FALSE,
  pt.size = 5,
  plot = TRUE, 
) 

d.b2 <- both$df_ord
ellips.b2 <- both$df_ellipse

both.p2 <- ggplot(d.b2, aes(x,y, fill = Group)) + 
  geom_point(size = 5,
             shape = 21,
             alpha = .8) +
  geom_polygon(data = ellips.b2, aes(color=Group, fill= Group), alpha = .4, show.legend = FALSE) +
  scale_fill_manual(values = my.col, labels = my.lab) +
  scale_color_manual(values = my.col) +
  xlab("CAP1") + ylab("CAP2") +
  theme_base() +
  theme(plot.background = element_rect(color = 'white'),
        legend.title = element_blank(),
        text = element_text(size = 16)
        # legend.text = element_text(size = 16),
        # axis.title = element_text(size = 16),
        # axis.text = element_text(size = 14)
  )

#panel----
library(patchwork)
# panel2 <- ((shape.p2 / fourier.p2) | both.p2) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A')
# ggsave(filename = "cap-shape-panels.pdf", plot = panel2, dpi = 600, width = 379, height = 189, units = "mm", device = "pdf")


#par(mfrow=c(1,3), mar = c(2.5, 2.5, 1, 7) + 0.1)

panel3 <- (shape.p2 | fourier.p2 | both.p2) + 
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A') 

panel3[[1]] = panel3[[1]] +
  theme(text = element_text(size = 28),
        legend.key.size = unit(1.1, "cm"),
        plot.margin = unit(c(0,60,0,0), "pt"))
panel3[[2]] = panel3[[2]] +
  theme(text = element_text(size = 28),
        legend.key.size = unit(1.1, "cm"),
        plot.margin = unit(c(0,60,0,0), "pt"))
panel3[[3]] = panel3[[3]] +
  theme(text = element_text(size = 28),
        legend.key.size = unit(1.1, "cm"),
        plot.margin = unit(c(0,60,0,0), "pt"))
            
ggsave(filename = "cap-shape-panels-se.pdf", plot = panel3, dpi = 600, width = 859, height = 259, units = "mm", device = "pdf")
