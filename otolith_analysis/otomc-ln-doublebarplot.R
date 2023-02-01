library(purrr)
library(foreach)
library(glue)
library(PMCMR)
library(multcompView)
library(dplyr)
library(tidyverse)
library(ggthemes)


setwd("./chem-april2021-ln")

#core-----
core <- read.csv("./otomc-corr-ln-core-10%na.csv", header = T, na.strings = "na")
core$population.code <- ordered(core$population.code,
                                     levels = c("GL", "NT", "NA", "A", "TU"))

library(rstatix)
kw.formula <- foreach(n=5:12) %do% as.formula(glue("{colnames(core)[n]} ~ population.code"))
res.kruskal.core <- foreach(n=1:8) %do% kruskal_test(core, kw.formula[[n]])
res.kruskal.core

post.hoc.formula <- foreach(n=5:12) %do% as.formula(glue("{colnames(core)[n]} ~ population.code"))
post.hoc.kw.core <- foreach(n=1:8) %do% posthoc.kruskal.dunn.test(post.hoc.formula[[n]], data = core, p.adjust = "bonf")
names(post.hoc.kw.core) <- colnames(core)[c(5:12)]

out.p.core <- map(post.hoc.kw.core, get.pvalues)
out.mcV.core <- map(out.p.core, ~multcompLetters(.x, threshold=0.05))



#edge----
edge <- read.csv("./otomc-corr-ln-edge-10%na.csv", header = T, na.strings = "na")
edge$population.code <- ordered(edge$population.code,
                                     levels = c("GL", "NT", "NA", "A", "TU"))


kw.formula <- foreach(n=5:12) %do% as.formula(glue("{colnames(edge)[n]} ~ population.code"))
res.kruskal.edge <- foreach(n=1:8) %do% kruskal_test(edge, kw.formula[[n]])
res.kruskal.edge

post.hoc.formula <- foreach(n=5:12) %do% as.formula(glue("{colnames(edge)[n]} ~ population.code"))
post.hoc.kw.edge <- foreach(n=1:8) %do% posthoc.kruskal.dunn.test(post.hoc.formula[[n]], data = edge, p.adjust = "bonf")
names(post.hoc.kw.edge) <- colnames(edge)[5:12]

out.p.edge <- map(post.hoc.kw.edge, get.pvalues)
out.mcV.edge <- map(out.p.edge, ~multcompLetters(.x, threshold=0.05))


#plot param----
my.col <- c("#1D4F9F","#7E277C","#F8941E","#80C342", "#ED1943")
my.col2 <- c("chocolate2", "cadetblue3")
#my.lab <- c("GL", "NT", "NA", "A", "TU")
my.lab <- c("GSA7_2009", "GSA9_2009",
            "GSA17_2009",
            "GSA22_2009", "GSA24_2009")
identical(colnames(core)[5:12],colnames(edge)[5:12])
my.vec <- foreach(n=1:8) %do% strsplit(colnames(core)[5:12], "_")[[n]][1]
my.vec <-  unlist(foreach(n=1:8) %do% strsplit(unlist(my.vec), "[1-9]")[[n]][1])


#plot----
chem.plot <- rbind(core, edge)
#chem.plot$population.code <- ordered(chem.plot$population.code, levels = my.lab)
chem.plot[,c(5:ncol(chem.plot))] <- exp(chem.plot[,c(5:ncol(chem.plot))])#exp fo log to have data as original for visualisation

my.ylab <- expression(paste(
  "Mean (",
  mu, g, " ", g^-1,
  ")", sep=""))



#create the datasets per element
k <- 1
data.list <- list()
for (i in 5:12){
  pop <- chem.plot[,2]
  ar <- chem.plot[,3]
  el <- chem.plot[,i]
  data <-  data.frame(population.code = pop,
                      area = ar,
                      conc = el)
  #data$element <- rank(data$element)
  #data <- melt(data, value.name = "conc")
  tim.core <- as.vector(table(core$population.code))
  tim.edge <- as.vector(table(edge$population.code))
  lab.core <- as.vector(out.mcV.core[[k]][[1]])
  lab.edge <- as.vector(out.mcV.edge[[k]][[1]])
  data$groups <- c(rep(lab.core, tim.core), rep(lab.edge, tim.edge))
  
  #add sd for barplot
  data <- data %>% group_by(area, population.code) %>% 
    dplyr::mutate(mean = mean(conc, na.rm = T),
                  sd = sd(conc, na.rm = T),
                  n = n(),
                  se = sd(conc, na.rm = T)/sqrt(n),
                  meanSE = mean_se(conc),
                  ylsd = mean + sd,
                  ylse = mean + se)
  data.list[[k]] <- data
  k <- k+1
}  


textpos.core <- c(480,4.5,0.085, 0.05,0.088,300, 1.10, 0.0035)
textpos.edge <- c(435,3.5,0.06, 0.05,0.015,450, 1.02, 0.02)
textpos <- c(textpos.core, textpos.edge)
ylv <- c(4300, 41, 0.82, 0.55, 1, 4200, 15, 0.145)
#plot per element
#OLD
k <- 1
plot.list <- list()
for (i in 1:length(data.list)){
  datat <- data.list[[i]]
  datac <- datat[datat$area == "core",]
  datae <- datat[datat$area == "edge",]
  #u <- max(datat$ylsd)+max(datat$sd)
  p <- ggplot(data= datat, aes(x= population.code,
                               y=mean,
                               fill=area))
  p <- p+ geom_bar(stat = "identity", position = position_dodge(0.85))
  p <- p +  geom_errorbar(data = datat,
                          aes(ymin=mean-se, ymax=mean+se),
                          width = 0.2,position = position_dodge(0.85))
  # p <- p +  geom_errorbar(data = datae,
  #                         aes(ymin=mean-se, ymax=mean+se),
  #                         width = 0.2,position = position_dodge(0.7))
  p <- p + geom_text(data = datat, 
                     aes(x=population.code, y=(mean+se), label=groups),
                     size = 3.8,position = position_dodge(0.85), vjust =-.7)
  # p <- p + geom_text(data = datae, 
  #                     aes(x=population.code, y=(mean+se), label=groups),
  #                     size = 5, position = position_nudge(y = textpos.edge[i]))
  p <- p + scale_fill_manual(values = my.col2)
  p <- p + xlab("") + ylab(my.ylab) + ggtitle(my.vec[k])
  p <- p + ylim(NA, ylv[i])
  p <- p + scale_x_discrete(breaks=c("GL", "NT", "NA", "A", "TU"),
                          labels=my.lab)
  #p <-  p + coord_flip()
  p <- p + theme_base()
  p <- p + theme(plot.background = element_rect(color='white'),
                 legend.position = "right", legend.title = element_blank(),
                 axis.text.x = element_text(size=9, angle=16),
                 axis.text.y = element_text(size=9),
                 text = element_text(size = 16),
                 plot.title = element_text(face = "plain"),
                 legend.key.size = unit(0.85, "cm")
                 # legend.text = element_text(size = 16),
                 # axis.title = element_text(size = 16),
                 # axis.text = element_text(size = 14)
  )
  
  plot.list[[k]] <- p
  k <- k+1
} #ylab concentration units!! error bar and significance of ANOVA

#all plots together
library(patchwork)
panel <- wrap_plots(plot.list) + plot_layout(guides = 'collect') + theme(plot.margin = unit(c(0,60,0,0), "pt"))

ggsave(filename = "barplot-ln-core-edge.pdf", plot = panel, path = "~/Desktop/manuscript-solea/graphs/", dpi = 600, width = 299, height = 210.748, units = "mm", device = "pdf")
