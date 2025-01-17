library("ggplot2")

d1<-read.csv("delme_depth1.tsv",sep="\t",header=FALSE)
colnames(d1)<-c("chr","pos","depth")

d100<-read.csv("delme_depth100.tsv",sep="\t",header=FALSE)
colnames(d100)<-c("chr","pos","depth")


ggplot(d1,aes(x=pos,y=depth)) + geom_point() + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+ labs(x="Position",y="Depth")
ggsave("simv1_depth1_subsampled.png")

ggplot(d100,aes(x=pos,y=depth)) + geom_point() + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+ labs(x="Position",y="Depth")
ggsave("simv1_depth100_subsampled.png")

d1<-read.csv("delme_depth1_v2.tsv",sep="\t",header=FALSE)
colnames(d1)<-c("chr","pos","depth")



ggplot(d1,aes(x=pos,y=depth)) + geom_point() + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+ labs(x="Position",y="Depth")
ggsave("simv1_depth1_subsampledv2.png")



library(dplyr)


p<-d100%>% 
ggplot(aes(x = pos, y = ifelse(depth==0,0,log10(depth)))) +
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) + 
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Position", y = "Depth (log10)")

p2<-d100%>% 
ggplot(aes(x = pos, y = log10(depth))) + 
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) + 
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Position", y = "Depth (log10)

ggsave(plot=p,filename="simv1_depth100_realdata_log.png")


ggsave(plot=p2,filename="simv1_depth100_realdata_log.png")
