
lut2<-c("0"="isConcordant","1"="isDiscordant")

#define HOM_TO_HOM 1
#define HET_TO_HET 2
#define HOM_TO_HET 3
#define HET_TO_HOM 4
lut<-c("1"="HOM_TO_HOM","2"="HET_TO_HET","3"="HOM_TO_HET","4"="HET_TO_HOM")



d1<-read.csv("delme1_persite",header=FALSE,sep="\t")
d2<-read.csv("delme2_persite",header=FALSE,sep="\t")


hdr<-c("Pos","Ind","isDiscordant","refToCallType")

colnames(d1)<-hdr
colnames(d2)<-hdr

library(dplyr)
# left_join(d1,d2,by=c("Pos","Ind")) -> x3

# x3$disagree<-(x3$isDiscordant.x != x3$isDiscordant.y)



# find the sites where "beta0" and "beta5" disagree
d3<-merge(d1,d2,by=c("Pos","Ind"),all=TRUE)
d3$disagree<-(d3$isDiscordant.x != d3$isDiscordant.y)

# find the sites where "beta0" and "beta5" disagree and "beta0" is discordant
(d3$disagree & d3$isDiscordant.x) -> d3$disagree_beta0discordant

# find the sites where "beta0" and "beta5" disagree and "beta5" is discordant
(d3$disagree & d3$isDiscordant.y) -> d3$disagree_beta5discordant


sum(d3$disagree_beta0discordant)
sum(d3$disagree_beta5discordant)

# > sum(d3$disagree_beta0discordant)
# [1] 24486
# > sum(d3$disagree_beta5discordant)
# [1] 9379


(!d3$disagree & d3$isDiscordant.x) -> d3$agree_discordant
(!d3$disagree & !d3$isDiscordant.x) -> d3$agree_concordant

sum(d3$agree_discordant)
sum(d3$agree_concordant)


# x = beta var 0
# y = beta var 1e-5

# > sum(d3$disagree_beta0discordant)
# [1] 24486
# > sum(d3$disagree_beta5discordant)
# [1] 9379


# number of sites where beta0 is discordant and beta5 is concordant
sum(d3$disagree_beta0discordant)

# number of sites where beta5 is discordant and beta0 is concordant
sum(d3$disagree_beta5discordant)

# print the rows where beta0 is discordant and beta5 is concordant
d3[d3$disagree_beta0discordant,]

# it is a bad site if beta0 is discordant and beta5 is concordant
d3$badSite<-(d3$isDiscordant.x & !d3$isDiscordant.y)
# d3[d3$badSite,]->d4
# d3$goodSite<-(!d3$isDiscordant.x & d3$isDiscordant.y)
# d3$goodSite<-NULL


# # [d4$Ind==0,][1,]

#           Pos Ind isDiscordant.x refToCallType.x Beta.x isDiscordant.y
# 10022162   0              1               4  beta0              0
#      refToCallType.y Beta.y disagree disagree_beta0discordant
# 2  beta5     TRUE                     TRUE
#      disagree_beta5discordant
# FALSE

# what is the frequency of different paste0(retToCallType.x,"_",refToCallType.y) values for sites where beta0 is discordant and beta5 is concordant



# get X from X_TO_Y
# sum(gsub("_.*","",lut[as.character(d3$refToCallType.x)]) != gsub("_.*","",lut[as.character(d3$refToCallType.y)]))
# == 0! OK!! it is all ok :)

d3$trueType<-gsub("_.*","",lut[as.character(d3$refToCallType.x)])
# get Y from X_TO_Y
d3$beta0call<-gsub(".*_","",lut[as.character(d3$refToCallType.x)])
d3$beta5call<-gsub(".*_","",lut[as.character(d3$refToCallType.y)])
d3$beta0status<-lut2[as.character(d3$isDiscordant.x)]
d3$beta5status<-lut2[as.character(d3$isDiscordant.y)]

d3 %>% group_by(beta0call,beta5call,trueType,beta0status,beta5status) %>% summarise(n=n()) -> dx

# add percentage to dx
dx %>% mutate(perc=n/sum(dx$n)) -> dx


# sum regardless trueType
dx %>% group_by(beta0status,beta5status) %>% summarise(n=sum(n), perc=sum(perc)) -> dx2

library(ggplot2)

dx2%>%
ggplot(aes(x=beta0status,y=beta5status))+
geom_tile()+
geom_text(aes(label=paste0("n=",n,"\n","perc=",round(perc*100,2),"%")),color="white")->plt

ggsave("wd/beta0_beta5_summary_perc.png",plt,width=10,height=10)

dx %>% group_by(trueType,beta0status,beta5status) %>% summarise(n=sum(n), perc=sum(perc)) -> dx3


dx3%>%
ggplot(aes(x=beta0status,y=beta5status,color="white"))+
facet_grid(.~trueType,labeller=label_both)+
geom_tile()+
coord_equal()+
labs(color="",fill="")+
#remove legend
theme(legend.position="none")+
geom_text(aes(label=paste0("n=",n,"\n","perc=",round(perc*100,2),"%")))->plt

ggsave("wd/beta0_beta5_summary_trueType_perc.png",plt,width=10,height=10)


dx %>% 
ggplot(aes(x=interaction(beta0status,beta0call),y=interaction(beta5status,beta5call),fill=trueType))+
geom_tile()+
#move legend to top
theme(legend.position="top")+
#add numbers as text
geom_text(aes(label=paste0("n=",n,"\n","perc=",round(perc*100,5),"%")),color="white")+
xlab("status.call for beta0")+
ylab("status.call for beta5")+
labs(fill="Type of the true genotype")+
facet_grid(.~trueType,scales="free",labeller=label_both)->plt

ggsave("wd/beta0_beta5_withPerc.png",plt,width=10,height=10)




#TSV: type, ind, contig, pos, read_index, error_prob
tsv1fn<-"wd/sim_vcfgl_2312-OutOfAfrica_3G09-chr22-rep0-d1-e0.002_qs2_5_qScores_QSERRPROB.tsv"
qep<-read.csv(tsv1fn,header=FALSE,sep="\t")
colnames(qep)<-c("type","Ind","contig","Pos","read_index","error_prob")
qep$contig<-NULL


# find the "error_prob" of all reads at the sites where beta0 is discordant and beta5 is concordant
# qep[qep$pos %in% d4$Pos ,]$error_prob -> bads_error_prob
# qep[!qep$pos %in% d4$Pos ,]$error_prob -> goods_error_prob
library(dplyr)


qep[qep$Pos %in% d3[d3$badSite,]$Pos ,]$error_prob -> bads_error_prob

qep$error_prob -> all_error_prob
qep[!qep$Pos %in% d3[d3$badSite,]$Pos ,]$error_prob -> notbads_error_prob

mean(bads_error_prob)
mean(notbads_error_prob)
mean(all_error_prob)

# > mean(bads_error_prob)
# [1] 0.001986238
# > mean(notbads_error_prob)
# [1] 0.002000724
# > mean(all_error_prob)
# [1] 0.001999568

qep$IndName<-qep$Ind
qep$Ind <- as.numeric(factor(qep$IndName,levels=unique(qep$IndName)))-1


# get the error probabilities for different combinations
d3 %>% left_join(qep,by=c("Pos","Ind")) -> qepd3

qepd3 %>%
# group_by(beta0call,beta5call,trueType,beta0status,beta5status) %>% 
group_by(beta0status,beta5status) %>% 
summarise(mean_error_prob=mean(error_prob,na.rm=TRUE)) 

qepd3 %>%
group_by(beta0status,beta5status) %>% 
summarise(mean_error_prob=mean(error_prob,na.rm=TRUE))

ggplot(qepd3,aes(y=error_prob,x=interaction(beta0status,beta5status)))+
geom_boxplot()+
geom_hline(yintercept=0.002,color="red")->plt

ggsave("wd/beta0_beta5_error_prob.png",plt,width=10,height=10)

qepd3 %>%
# group_by(beta0call,beta5call,trueType,beta0status,beta5status) %>% 
mutate(err_prob_diff=abs(0.002-error_prob)) %>%
group_by(beta0status,beta5status) %>% 
summarise(mean_abs_error_prob_diff=mean(err_prob_diff,na.rm=TRUE))


qepd3 %>%
group_by(beta0call,beta5call,trueType,beta0status,beta5status) %>% 
summarise(mean_error_prob=mean(error_prob,na.rm=TRUE)) -> qepd3b

qepd3b %>%
ggplot(aes(x=interaction(beta0status,beta0call),y=interaction(beta5status,beta5call),fill=trueType))+
geom_tile()+
facet_grid(.~trueType,labeller=label_both)+
theme(legend.position="top")+
labs(fill="Type of the true genotype",
        x="status.call for beta0",
        y="status.call for beta5",
        title="Mean raw error probability sampled from beta distribution")+
        #rotate x axis labels 90 degrees
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
geom_text(aes(label=paste0("mean\nerror prob\n",round(mean_error_prob,7))),color="white")->plt

ggsave("wd/beta0_beta5_error_prob_detailed.png",plt,width=10,height=10)



# library(ggplot2)
# # compare bads and goods prob with boxplots
# ggplot()+
#     geom_boxplot(aes(x="bad",y=bads_error_prob),data=data.frame(bads_error_prob))+
#     geom_boxplot(aes(x="all",y=all_error_prob),data=data.frame(all_error_prob))+
#     geom_boxplot(aes(x="notbad",y=notbads_error_prob),data=data.frame(notbads_error_prob))->plt

# ggsave("wd/bads_error_probs.png",plt,width=10,height=10)





#TSV: type, ind, contig, pos, read_index, qScore
tsv2fn<-"wd/sim_vcfgl_2312-OutOfAfrica_3G09-chr22-rep0-d1-e0.002_qs2_5_qScores_QS.tsv"
qs<-read.csv(tsv2fn,header=FALSE,sep="\t")
colnames(qs)<-c("type","ind","contig","pos","read_index","qScore")
qs$contig<-NULL

qs$IndName<-qs$ind
qs$Ind <- as.numeric(factor(qs$IndName,levels=unique(qs$IndName)))-1

# find the "qScore" of all reads at the sites where beta0 is discordant and beta5 is concordant
qs[qs$pos %in% d3[d3$badSite,]$Pos ,]$qScore -> bads_qScore
qs[!qs$pos %in% d3[d3$badSite,]$Pos ,]$qScore -> rest_qScore
qs$qScore -> all_qScore

lut3<-c("FALSE_FALSE_TRUE"="bt27","FALSE_TRUE_FALSE"="eq27","TRUE_FALSE_FALSE"="st27")
sum(qs$qScore>27)
sum(qs$qScore<27)
sum(qs$qScore==27)
# summarise the nsites <27, >27, =27, join those new cols, then get the ratios
qs %>% group_by(st27=qScore<27,eq27=qScore==27,bt27=qScore>27) %>% summarise(n=n()) %>% ungroup() %>% mutate(cmp27=paste0(st27,"_",eq27,"_",bt27)) %>% mutate(cmp27=lut3[cmp27]) %>% group_by(cmp27) %>% summarise(n=sum(n)) %>% mutate(perc=n/sum(n))


mean(bads_qScore)
mean(rest_qScore)
mean(all_qScore)


ggplot()+
    geom_boxplot(aes(x="bad",y=bads_qScore),data=data.frame(bads_qScore))+
    geom_boxplot(aes(x="all",y=all_qScore),data=data.frame(all_qScore))+
    geom_hline(yintercept=27,color="black")+
    geom_boxplot(aes(x="rest",y=rest_qScore),data=data.frame(rest_qScore))->plt


ggsave("wd/bads_qScores.png",plt,width=10,height=10)

ggplot()+geom_histogram(aes(x=bads_qScore-27),data=data.frame(bads_qScore-27),binwidth=1,fill="red")

ggplot()+geom_histogram(aes(x=goods_qScore-27),data=data.frame(goods_qScore-27),binwidth=1,fill="green")


badsqs<-data.frame(qScore=bads_qScore)
restqs<-data.frame(qScore=rest_qScore)
allqs<-data.frame(qScore=all_qScore)

badsqs$group<-"bad"
restqs$group<-"rest"
allqs$group<-"all"
rbind(badsqs,restqs,allqs)->allqs

# allqs%>% 
# ggplot(aes(x=qScore-27,fill=group))+
# facet_grid(.~group)+
# geom_histogram(binwidth=1,alpha=0.3,position="identity")->plt
# ggsave("wd/bads_qScores_hist.png",plt,width=10,height=10)
# ggsave("wd/bads_qScores_hist.png",plt,width=10,height=10)

qs$Pos<-qs$pos

# qsd3<-merge(qs,d3,by=c("Pos","Ind"),all=TRUE)
# left_join(qs,d3,by=c("Pos","Ind")) -> qsd3
# _join(qs,d3,by=c("Pos","Ind")) -> qsd3



right_join(qs,d3,by=c("Pos","Ind")) -> qsd3

qsd3 %>% group_by(beta0status,beta5status) %>% summarise(mean_qScore=mean(qScore,na.rm=TRUE))



mean(qsd3$qScore,na.rm=TRUE)



mean(qs$qScore)


sum(qs$qScore>27)
sum(qs$qScore<27)
sum(qs$qScore==27)



sum(qsd3$qScore>27,na.rm=TRUE)

qsd3$qScoreDiff<-qsd3$qScore-27

qsd3 %>% group_by(beta0status,beta5status) %>% summarise(mean_qScoreDiff=mean(qScoreDiff,na.rm=TRUE))

qsd3 %>% group_by(beta0status,beta5status,trueType) %>% summarise(mean_qScoreDiff=mean(qScoreDiff,na.rm=TRUE))

qsd3 %>% group_by(beta0status,beta5status,trueType,beta0call,beta5call) %>% summarise(mean_qScoreDiff=mean(qScoreDiff,na.rm=TRUE))  %>% arrange(desc(mean_qScoreDiff))


-10*log10(0.002)+0.499

qsd3 %>%  mutate(abs_qScoreDiff=abs(qScoreDiff)) %>% group_by(beta0status,beta5status,trueType,beta0call,beta5call) %>% summarise(mean_abs_qScoreDiff=mean(abs_qScoreDiff,na.rm=TRUE))  %>% arrange(desc(mean_abs_qScoreDiff))

qsd3 %>%  mutate(abs_qScoreDiff=abs(qScoreDiff)) %>% 
group_by(beta0status,beta5status,trueType,beta0call,beta5call) %>% 
summarise(mean_abs_qScoreDiff=mean(abs_qScoreDiff,na.rm=TRUE))  %>% arrange(desc(mean_abs_qScoreDiff))

qsd3 %>%  mutate(abs_qScoreDiff=abs(qScoreDiff)) %>% 
group_by(beta0status,beta5status,trueType) %>%
summarise(mean_abs_qScoreDiff=mean(abs_qScoreDiff,na.rm=TRUE))  %>% arrange(desc(mean_abs_qScoreDiff))

qsd3 %>%  mutate(abs_qScoreDiff=abs(qScoreDiff)) %>% 
group_by(beta0status,beta5status) %>%
summarise(mean_abs_qScoreDiff=mean(abs_qScoreDiff,na.rm=TRUE))  %>% arrange(desc(mean_abs_qScoreDiff))


mean(bads_qScore-27)
mean(rest_qScore-27)
mean(all_qScore-27)
mean(27-bads_qScore)
mean(27-rest_qScore)
mean(27-all_qScore)

# I also checked the quality scores of the sites where beta0 is discordant and beta5 is concordant. This is important to be able to see the actual effect the beta-sampled error probabilities may have. The mean quality score is 




mis1fn="wd/in1_pos_misRatio.tsv"
mis2fn="wd/in2_pos_misRatio.tsv"
mis1<-read.csv(mis1fn,header=FALSE,sep=" ")
mis2<-read.csv(mis2fn,header=FALSE,sep=" ")
colnames(mis1)<-c("Pos","misRatio")
colnames(mis2)<-c("Pos","misRatio")
mis1$Beta<-"beta0"
mis2$Beta<-"beta5"

merge(mis1,mis2,by="Pos",all=TRUE)->mis


d3 %>% left_join(mis,by=c("Pos","Beta.x","Beta.y")) -> d3mis

d3mis %>% group_by(beta0status,beta5status) %>% 
summarise(mean_misRatio0=mean(misRatio.x,na.rm=TRUE),mean_misRatio5=mean(misRatio.y,na.rm=TRUE)) %>% arrange(desc(mean_misRatio0))

d3mis %>% group_by(beta0status,beta5status,trueType) %>%
summarise(mean_misRatio0=mean(misRatio.x,na.rm=TRUE),mean_misRatio5=mean(misRatio.y,na.rm=TRUE)) %>% arrange(desc(mean_misRatio0))

d3mis %>% group_by(beta0status,beta5status,trueType,beta0call,beta5call) %>%
summarise(mean_misRatio0=mean(misRatio.x,na.rm=TRUE),mean_misRatio5=mean(misRatio.y,na.rm=TRUE)) %>% arrange(desc(mean_misRatio0))



mean_misRatio0 = mean missingness ratio for beta0 at these sites
mean_misRatio5 = mean missingness ratio for beta5 at these sites
beta0status = isDiscordant/isConcordant status for beta0
beta5status = isDiscordant/isConcordant status for beta5
trueType = true genotype type (HOM/HET)
beta0call = call type for beta0 (HOM/HET)
beta5call = call type for beta5 (HOM/HET)

# find the "misRatio" of all reads at the sites where beta0 is discordant and beta5 is concordant
# mis1[mis1$Pos %in% d4$Pos ,]$misRatio -> bads_misRatio
# mis1[mis1$Pos %in% d3[d3$badSite,]$Pos ,]$misRatio -> bads_misRatio
# mis1$misRatio -> all_misRatio 
# mis1[!mis1$Pos %in% d3[d3$badSite,]$Pos ,]$misRatio -> rest_misRatio

mean_error_prob = mean error probability sampled from beta distribution (raw probability value) for beta5 runs at these sites
mean_qScore = mean qScore for beta5 runs at these sites (qScore is calculated from the error probability sampled from beta distribution with the phred transformation using the formula -10*log10(p)+0.499

mean(bads_misRatio)
mean(all_misRatio)
mean(rest_misRatio)



# compare bads and goods misRatio with boxplots
ggplot()+
    geom_boxplot(aes(x="bad",y=bads_misRatio),data=data.frame(bads_misRatio))+
    geom_boxplot(aes(x="good",y=goods_misRatio),data=data.frame(goods_misRatio))+
    geom_hline(yintercept=0.002,color="black")+
    geom_hline(yintercept=mean(bads_misRatio),color="red")+
    geom_hline(yintercept=mean(goods_misRatio),color="green")->plt

ggsave("wd/bads_misRatio.png",plt,width=10,height=10)

ggplot()+
geom_histogram(aes(x=bads_misRatio),data=data.frame(bads_misRatio),binwidth=0.0001,fill="red")+

# Hi,

# I checked examples of sites where beta variance=0 parameter runs genotype is discordant while beta variance=1e-5 is concordant. 
Below, when I use the term "number of sites", I mean the number of sites * number of individuals as each individual has one genotype call for each site.

beta0_isRight_beta5_isRight beta0_isRight_beta5_isWrong 
                   10731039                        9379 
beta0_isWrong_beta5_isRight beta0_isWrong_beta5_isWrong 
                      24486                      771448 

beta0_isRight_with_HET_TO_HET_beta5_isRight_with_HET_TO_HET 
                                                     680146 
beta0_isRight_with_HET_TO_HET_beta5_isWrong_with_HET_TO_HET 
                                                         90 
beta0_isRight_with_HET_TO_HET_beta5_isWrong_with_HET_TO_HOM 
                                                       3411 
beta0_isRight_with_HOM_TO_HOM_beta5_isRight_with_HOM_TO_HOM 
                                                   10050893 
beta0_isRight_with_HOM_TO_HOM_beta5_isWrong_with_HOM_TO_HET 
                                                       5878 
beta0_isWrong_with_HET_TO_HOM_beta5_isRight_with_HET_TO_HET 
                                                      22706 
beta0_isWrong_with_HET_TO_HOM_beta5_isWrong_with_HET_TO_HET 
                                                        222 
beta0_isWrong_with_HET_TO_HOM_beta5_isWrong_with_HET_TO_HOM 
                                                     558213 
beta0_isWrong_with_HOM_TO_HET_beta5_isRight_with_HOM_TO_HOM 
                                                       1780 
beta0_isWrong_with_HOM_TO_HET_beta5_isWrong_with_HOM_TO_HET 
                                                     211766 
beta0_isWrong_with_HOM_TO_HET_beta5_isWrong_with_HOM_TO_HOM 
                                                         72 
beta0_isWrong_with_HOM_TO_HOM_beta5_isWrong_with_HOM_TO_HET 
                                                        457 
beta0_isWrong_with_HOM_TO_HOM_beta5_isWrong_with_HOM_TO_HOM 
                                                        718 



For replicate=0 depth=1, I compared beta variance=0 and beta variance=1e-5 genotype calls for all sites.
Number of sites where beta0 is discordant and beta1e-5 is concordant: 24486
Number of sites where beta1e-5 is discordant and beta0 is concordant: 9379


I compared the missingness ratio (number of individuals with no data / total number of individuals) of the sites where (1) beta0 is discordant and beta1e-5 is concordant and (2) rest of the sites.




Also, in the meeting I said there is a big difference in the number of SNP sites between beta variance=0 and beta variance=1e-5 runs. Since it being big is a bit subjective, I am providing the numbers below:
Number of SNPs in beta variance=0 run: 92074
Number of SNPs in beta variance=1e-5 run: 107298
16.5 % increase in number of SNPs at beta variance=1e-5 run compared to beta variance=0 run




In terms of discordance type,
Number of sites where beta0=discordant beta5=concordant (in all of these beta0 misidentifies the genotype as HOM while the true genotype is HET): 22706
Number of sites where beta0=concordant beta5=discordant (in all of these beta5 misidentifies the 



In the attached plots, you can find the count of sites for different combinations of beta0 and beta5 status (isDiscordant,isConcordant) and beta0 and beta5 genotype calls (HOM,HET) and the true genotype (HOM,HET). These are the numbers where I compare the individual genotype calls where the only difference is in one scenario quality score errors were not added (beta0) and in the other scenario quality score errors were added (beta5). All the other parameters, such as read depth and base counts etc are the same (I used separate RNG for beta distribution so that we can keep all the other parameters fixed).


I also paste the table below, sorted by percentage:

  beta0status  beta5status         n     perc
1 isConcordant isConcordant 10731039 0.930   
2 isDiscordant isDiscordant   771448 0.0669  
3 isDiscordant isConcordant    24486 0.00212 
4 isConcordant isDiscordant     9379 0.000813


 With regards to the true type and call type:

   beta0call beta5call trueType beta0status  beta5status         n       perc
 1 HOM       HOM       HOM      isConcordant isConcordant 10050893 0.871     
 2 HET       HET       HET      isConcordant isConcordant   680146 0.0590    
 3 HOM       HOM       HET      isDiscordant isDiscordant   558213 0.0484    
 4 HET       HET       HOM      isDiscordant isDiscordant   211766 0.0184    
 5 HOM       HET       HET      isDiscordant isConcordant    22706 0.00197   
 6 HOM       HET       HOM      isConcordant isDiscordant     5878 0.000510  
 7 HET       HOM       HET      isConcordant isDiscordant     3411 0.000296  
 8 HET       HOM       HOM      isDiscordant isConcordant     1780 0.000154  
 9 HOM       HOM       HOM      isDiscordant isDiscordant      718 0.0000622 
10 HOM       HET       HOM      isDiscordant isDiscordant      457 0.0000396 
11 HOM       HET       HET      isDiscordant isDiscordant      222 0.0000192 
12 HET       HET       HET      isConcordant isDiscordant       90 0.00000780
13 HET       HOM       HOM      isDiscordant isDiscordant       72 0.00000624

For the sites where beta0 is discordant and beta5 is concordant, the mean error probability sampled from beta distribution is 0.001986238. For the rest of the sites, the mean error probability sampled from beta distribution is 0.002000724. For all sites, the mean error probability sampled from beta distribution is 0.001999568. This is the raw values sampled from beta distribution. 

if we compare it for each site for beta0 and beta5:
  beta0status  beta5status  mean_error_prob
1 isConcordant isConcordant         0.00200
2 isConcordant isDiscordant         0.00201
3 isDiscordant isConcordant         0.00198
4 isDiscordant isDiscordant         0.00200

  beta0status  beta5status  mean_qScore
  <chr>        <chr>              <dbl>
1 isConcordant isConcordant        33.9
2 isConcordant isDiscordant        33.9
3 isDiscordant isConcordant        34.0
4 isDiscordant isDiscordant        33.9

  beta0status  beta5status  mean_qScoreDiff
  <chr>        <chr>                  <dbl>
1 isConcordant isConcordant            6.87
2 isConcordant isDiscordant            6.86
3 isDiscordant isConcordant            7.01
4 isDiscordant isDiscordant            6.87

  beta0status  beta5status  trueType mean_qScoreDiff
  <chr>        <chr>        <chr>              <dbl>
1 isConcordant isConcordant HET                 6.86
2 isConcordant isConcordant HOM                 6.87
3 isConcordant isDiscordant HET                 7.00
4 isConcordant isDiscordant HOM                 6.78
5 isDiscordant isConcordant HET                 7.05
6 isDiscordant isConcordant HOM                 6.53
7 isDiscordant isDiscordant HET                 6.87
8 isDiscordant isDiscordant HOM                 6.84


Please also see the attached q-score and error probability (raw values sampled from beta distribution) distributions. The majority of the q-scores have an over-estimated quality score. This is due to the fact that we sample probabilities from the beta distribution and then convert them to q-scores. The q-scores themselves are not sampled from the beta distribution. Since the beta distribution is skewed towards 0, the majority of the q-scores are over-estimated. This is the reason why the mean q-score is 33.9 for all sites, while the mean error probability sampled from beta distribution is 0.001999568.



qs_error_prob	popYRI_ind1	chr22	10022162	0	0.000521
qs	popYRI_ind1	chr22	10022162	0	33

An example site where beta0 is discordant and beta5 is concordant:


beta=0, gl=1
chr22	10022162	.	A	C,<*>	.	PASS	QS=59.3333,0.666667,0;AD=90,2,0	DP:GL:PL:AD	3:-1.79777,0,-4.0523,-2.39984,-4.35333,-6.27605:18,0,41,24,44,63:2,1,0

beta=5, gl=1
chr22	10022162	.	A	C,<*>	.	PASS	QS=59.2716,0.728428,0;AD=90,2,0	DP:GL:PL:AD	3:-2.39712,0,-6.31703,-2.99919,-6.61806,-9.14013:24,0,63,30,66,91:2,1,0

beta=0, gl=2
chr22	10022162	.	A	C,<*>	.	PASS	QS=59.3333,0.666667,0;AD=90,2,0	DP:GL:PL:AD	3:-2.27403,0,-5.45028,-2.87551,-5.75103,-8.62654:23,0,55,29,58,86:2,1,0

beta=5, gl=2
chr22	10022162	.	A	C,<*>	.	PASS	QS=59.2716,0.728428,0;AD=90,2,0	DP:GL:PL:AD	3:-2.87462,0,-7.64977,-3.47595,-7.95073,-11.4267:29,0,76,35,80,114:2,1,0

(gl1=error-dependent gl model, gl2=direct gl model)
I put the gl 2 values too to show that we see a similar pattern among gl models in terms of the calculated gl values with different beta variances

beta=5, gl=1
chr22	10022162	.	A	C	20.7068	PASS	AD=90,2;AC=1;AN=2	GT:DP:GL:PL:AD:GP:GQ	0/1:3:-2.39712,0,-6.31703,-2.99919,-6.61806,-9.14013:24,0,63:2,1:0.139391,0.860609,2.65043e-09:8


beta=0, gl=1
(site is invariable in beta=0 run)
chr22	10022162	.	A	.	17.8655	PASS	AD=90;AN=2	GT:DP:GL:AD	0/0:3:-1.79777,0,-4.0523,-2.39984,-4.35333,-6.27605:2


beta 



pileup for beta0:
chr22	10022162	A	3	CAA	<<<

pileup for beta5:
chr22	10022162	A	3	CAA	B8V





# 240118

library(ggplot2)
library(dplyr)

alpha <- 0.397200
beta <- 198.202800

# alpha=399.198000 and beta=199199.802000 
# alpha <- 399.198000
# beta <- 199199.802000

num_trials <- 1e6
bin_size <- 1e-5
bins <- seq(0, 1, by = bin_size)
#remove first val (0->later inf)
bins <- bins[-1]
density <- dbeta(bins, alpha, beta)

expected_hits <- density * bin_size * num_trials

d <- data.frame(Bin = bins, ExpectedHits = expected_hits)

summary(expected_hits)

#remove inf values
d[is.infinite(d$ExpectedHits),]

sum(d[d$Bin<0.0005,]$ExpectedHits)/sum(d$ExpectedHits)



btRatio<-sum(d[d$Bin>0.002,]$ExpectedHits)/sum(d$ExpectedHits)
stRatio<-sum(d[d$Bin<0.002,]$ExpectedHits)/sum(d$ExpectedHits)
eqRatio<-sum(d[d$Bin==0.002,]$ExpectedHits)/sum(d$ExpectedHits)
max1em5<-max(d[d$ExpectedHits>1e-5,]$Bin)

# d%>% filter( Bin<max(d[d$ExpectedHits > 0,]$Bin))%>% 
d%>%
filter(ExpectedHits>0)%>%
ggplot(aes(x = Bin, y = ExpectedHits)) +
xlim(0,max1em5)+
geom_point()+
geom_vline(xintercept = 0.002,color="red")+
annotate("text",x=0.002,y=0.0000001,label="0.002",color="red")




ggplot(d, aes(x = Bin, y = ExpectedHits)) +
geom_point()
  # geom_line() +
  labs(title = "Expected Number of Hits in Each Bin",
       x = "Bin",
       y = "Expected Number of Hits")






rm(list=ls());gc()


library(dplyr)
library(ggplot2)
library(tidyr)

#TSV: type, ind, contig, pos, read_index, error_prob
tsv1fn<-"wd/sim_vcfgl_2312-OutOfAfrica_3G09-chr22-rep0-d1-e0.002_qs2_5_qScores_QSERRPROB.tsv"
qep<-read.csv(tsv1fn,header=FALSE,sep="\t")
colnames(qep)<-c("Type","Ind","Contig","Pos","Read","e")
qep$contig<-NULL

qep$IndName<-qep$Ind
qep$Ind <- as.numeric(factor(qep$IndName,levels=unique(qep$IndName)))-1




lut2<-c("0"="isConcordant","1"="isDiscordant")

#define HOM_TO_HOM 1
#define HET_TO_HET 2
#define HOM_TO_HET 3
#define HET_TO_HOM 4
lut<-c("1"="HOM_TO_HOM","2"="HET_TO_HET","3"="HOM_TO_HET","4"="HET_TO_HOM")






library(dplyr)
d1<-read.csv("~/Projects/VCFGL/vcfgl_test_mcall_240118/discordance_qs0_persite",header=FALSE,sep="\t")
d2<-read.csv("~/Projects/VCFGL/vcfgl_test_mcall_240118/discordance_qs5_persite",header=FALSE,sep="\t")
hdr<-c("Pos","Ind","isDiscordant","refToCallType")
colnames(d1)<-hdr
colnames(d2)<-hdr
d1%>%mutate(isDiscordant=1==isDiscordant)%>%summarise(DiscordanceRate=sum(isDiscordant)/nrow(d1))
d2%>%mutate(isDiscordant=1==isDiscordant)%>%summarise(DiscordanceRate=sum(isDiscordant)/nrow(d2))



d<-full_join(d1,d2,by=c("Pos","Ind"),suffix=c("_qs0","_qs5"))

d$disagree<-(d$isDiscordant_qs0 != d$isDiscordant_qs5)

d$isDiscordant_qs0<-1==d$isDiscordant_qs0
d$isDiscordant_qs5<-1==d$isDiscordant_qs5


# save(d,file="d.RData")

d%>%group_by(lut[refToCallType_qs0],lut[refToCallType_qs5],disagree)%>%summarise(n=n())%>%ungroup()%>%mutate(perc=(n/sum(n))*100)%>%arrange(desc(perc))
d%>%group_by(lut[refToCallType_qs0],lut[refToCallType_qs5],disagree)%>%summarise(n=n())%>%ungroup()%>%mutate(perc=(n/sum(n))*100)%>%arrange(desc(perc))


d%>%group_by(isDiscordant_qs0,isDiscordant_qs5)%>%summarise(n=n())%>%ungroup()%>%mutate(perc=(n/sum(n))*100)%>%arrange(desc(perc))

#get discordance rate for qs0
d%>%group_by(isDiscordant_qs0)%>%summarise(n=n())%>%ungroup()%>%mutate(perc=(n/sum(n))*100)%>%arrange(desc(perc))

















d %>% group_by(disagree,isDiscordant_qs0,isDiscordant_qs5) %>% summarise(mean_e=mean(e,na.rm=TRUE)) %>% arrange(desc(mean_e))


#TSV: type, ind, contig, pos, read_index, qScore
tsv2fn<-"wd/sim_vcfgl_2312-OutOfAfrica_3G09-chr22-rep0-d1-e0.002_qs2_5_qScores_QS.tsv"
qs<-read.csv(tsv2fn,header=FALSE,sep="\t")
colnames(qs)<-c("type","ind","contig","pos","read_index","qScore")
qs$contig<-NULL

qs$IndName<-qs$ind
qs$Ind <- as.numeric(factor(qs$IndName,levels=unique(qs$IndName)))-1

