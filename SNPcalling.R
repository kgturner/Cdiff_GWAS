#Quality checks, pre-SNP calling in dDocent
#2/21/2017

#John Purtiz sez:  The best way to find the similarity cutoff is to plot # of 
# contigs vs. % similarity, making a line for each combination of K1 and K2.  
# You should see that each curve should have a similar inflection point, 
# even though the absolute values differ.

kopt <- read.delim("kopt.data", header = FALSE, sep=" ")
head(kopt)
names(kopt)[1] <- "k1"
names(kopt)[2] <- "k2"
names(kopt)[3] <- "percSim"
names(kopt)[4] <- "NoContigs"
kopt$ID <- paste0(kopt$k1, kopt$k2)
kopt$ID <- as.factor(kopt$ID)
subkopt <- subset(kopt, NoContigs>3000)

library(ggplot2)
p <- qplot(data = kopt, x = percSim, y=NoContigs)
p2 <- ggplot(kopt, aes(x=percSim, y=NoContigs, group=ID, color=ID))+
  geom_line()+
  geom_point()
p2

p3 <- ggplot(subkopt, aes(x=percSim, y=NoContigs, group=ID, color=ID))+
  geom_line()+
  geom_point()
p3
ggsave("Cdiff_GBS_qaplot.png", height = 10)
