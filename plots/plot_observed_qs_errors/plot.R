# isinaltinkaya

###############################################################################
# PRODUCE DATA

# # vcfgl_exec="/home/isin/Projects/VCFGL/vcfgl/vcfgl"
# vcfgl_exec="/maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_paper_analyses/tools/vcfgl/vcfgl"
# infile="/maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_paper_analyses/simulation/sim/sim_vcfgl_2312/model_OutOfAfrica_3G09/contig_chr22/vcf/sim_vcfgl_2312-OutOfAfrica_3G09-chr22.vcf"

# betavar=c(0,5,6,7)
# preciseopt=c(0,1)
# alld<-NULL
# for(i in seq_along(betavar)){
#   for (j in seq_along(preciseopt)){
#   if(betavar[i]==0){
#     system(paste0(vcfgl_exec," -i ",infile," -d 0.1 --error-qs 0 -e 0.002 --seed 42 -printGlError 1 -printQsError 1 -printQScores 1 -printBasePickError 1  --precise-gl ",preciseopt[j]," > qscores_",betavar[i],"_",preciseopt[j]," 2> qscores_",betavar[i],"_",preciseopt[j],".log"))
#   }else{
#     system(paste0(vcfgl_exec," -i ",infile," -d 0.1 --error-qs 2 -e 0.002 --beta-variance 1e-",betavar[i]," --seed 42 -printGlError 1 -printQsError 1 -printQScores 1 -printBasePickError 1 --precise-gl ",preciseopt[j]," > qscores_",betavar[i],"_",preciseopt[j]," 2> qscores_",betavar[i],"_",preciseopt[j],".log"))
#   }
#   dd<-read.csv(paste0("qscores_",betavar[i],"_",preciseopt[j]),sep="\t",header=FALSE)
#   colnames(dd)<-c("type","sample_id","contig","site","read_index","value")
#   dd$betavar=betavar[i]
#   dd$precise=preciseopt[j]
#   # dd<-data.frame(probs=(10 ^ ((-1*d$V1)/10)), qs=d$V1)
#   alld<-rbind(alld,dd)
# }
# }



# alld$type<-as.factor(alld$type)

# alld$betavar<-as.factor(alld$betavar)

# save(alld,file="data.RData")


##########################i#####################################################
# LOAD LIBRARIES

library(tidyr)
library(ggplot2)
library(dplyr)

###############################################################################
# FUNCTIONS

pToQ<-function(p){floor(-10*log10(p)+0.499)}
qToP<-function(q){10^(-q/10)}

parse_betavar <- function(x){
  for(i in seq_along(levels(x))){
    if(levels(x)[i] != 0){
      levels(x)[i]<-gsub("X","10^-",paste0("X",levels(x)[i]))
    }
  }
  x
}

precise_lut<-c("0"="Discretized GL error","1"="Precise GL error")
type_lut<-c("gl_error_prob"="e used in GL calculation","qs_error_prob"="e used in Q-score calculation","base_pick_error_prob"="e used in base picking","qs"="Q-score")


labeller1fn <- function(variable,value){
  if (variable=='betavar') {
      # return(label_parsed(gsub("X","10^-",paste0("X",value))))
      return(label_parsed(parse_betavar(value)))
    # return(as.character(value))
  }else if(variable=='precise'){
    precise_lut[as.character(value)]
  }else if (variable=='type'){
    return(type_lut[as.character(value)])
  }else{
    return(as.character(value))
  }
}

plot_outdir<-getwd()

save_plt <- function(plt,plot_filename,plot_outdir,overwrite=FALSE,...){
  if(!grepl("/$",plot_outdir)){
    plot_outdir<-paste0(plot_outdir,"/")
  }
  dir.create(plot_outdir, showWarnings = FALSE, recursive = TRUE)
  plot_fullpath<-paste0(plot_outdir,plot_filename)
  if(file.exists(plot_fullpath)){
    if(overwrite){
      print(paste0("Overwriting ",plot_fullpath))
      ggsave(plot_fullpath, plt,...)
    }else{
      stop(paste0("\n\n\n\n\nFile ",plot_fullpath," already exists.\n\n\n\n\n"))
    }
  }else{
    print(paste0("Saving plot to ",plot_fullpath))
    ggsave(plot_fullpath, plt,...)
  }
}

###############################################################################
# LOAD DATA

load("data.RData")


basepick_error_prob<-unique(alld[alld$type=="base_pick_error_prob",]$value)


color_betavars<-scale_color_manual(labels=c("0"="0","7"=expression(10^-7), "6"=expression(10^-6),"5"=expression(10^-5)), values=c("0"="black","7"="green","6"="blue","5"="red"))
fill_betavars<-scale_fill_manual(labels=c("0"="0","7"=expression(10^-7), "6"=expression(10^-6),"5"=expression(10^-5)), values=c("0"="black","7"="green","6"="blue","5"="red"))


# ###############################################################################



alld %>% 
filter(betavar!=0)%>%
filter(type=="qs_error_prob")%>%
filter(precise==0)%>%
#calculate observed mean for each betavar
group_by(betavar)%>%
summarise(mean=mean(value))->observed_mean_qs_error_prob
#OK: all equal to 0.002



alld %>% 
filter(betavar!=0)%>%
filter(type=="qs_error_prob")%>%
filter(precise==0)%>%
ggplot(aes(x=value,fill=betavar,group=betavar))+
stat_bin(alpha=0.7,binwidth=0.0001,position="identity")+
# geom_bar(aes(y=after_stat(prop),group=betavar))+
# facet_wrap(.~betavar,scale="free",labeller=labeller1fn)+
theme_bw()+
# xlim(0,0.02)+
theme(strip.background = element_rect(colour="black", fill="white"))+
geom_vline(xintercept=basepick_error_prob,color="black",size=0.5,linetype="dashed")+
labs(
  x="Error probability used in Q-score calculation",
  y="Count",
  fill="Beta variance"
)+
theme(legend.position="top")+
fill_betavars+
guides(color=FALSE)



save_plt(plt=last_plot(),plot_filename="plot_histogram_qs_error_prob.png",plot_outdir=plot_outdir,overwrite=TRUE,width=10,height=5,scale=1.5,dpi=300,units="cm")

true_qScore<-pToQ(basepick_error_prob)

alld %>% 
filter(betavar!=0)%>%
filter(type=="qs")%>%
filter(precise==0)%>%
mutate(across(c('value'),as.factor))%>%
ggplot(aes(x=value,fill=betavar,group=betavar))+
# stat_count(alpha=0.7, width=1)+
geom_bar(position="identity", alpha=0.5)+
geom_vline(xintercept=as.factor(true_qScore),color="black",size=0.5, linetype="dashed")+
theme_bw()+
theme(strip.background = element_rect(colour="black", fill="white"))+
labs(
  x="Q-score",
  y="Count",
  fill="Beta variance"
)+
theme(legend.position="top")+
fill_betavars+
guides(color=FALSE)+
scale_x_discrete(breaks=seq(0,max(alld[alld$type=="qs",]$value),2))



save_plt(plt=last_plot(),plot_filename="plot_histogram_qs.png",plot_outdir=plot_outdir,overwrite=TRUE,width=12,height=10,scale=1.5,dpi=300,units="cm")


gldata0<-alld[alld$betavar==0 & alld$type=="gl_error_prob",]
gldata5<-alld[alld$betavar==5 & alld$type=="gl_error_prob",]
gldata6<-alld[alld$betavar==6 & alld$type=="gl_error_prob",]
gldata7<-alld[alld$betavar==7 & alld$type=="gl_error_prob",]
gldata0$gldiff<-gldata0$value-basepick_error_prob
gldata5$gldiff<-gldata5$value-basepick_error_prob
gldata6$gldiff<-gldata6$value-basepick_error_prob
gldata7$gldiff<-gldata7$value-basepick_error_prob

qsdata0<-alld[alld$betavar==0 & alld$type=="qs_error_prob",]
qsdata5<-alld[alld$betavar==5 & alld$type=="qs_error_prob",]
qsdata6<-alld[alld$betavar==6 & alld$type=="qs_error_prob",]
qsdata7<-alld[alld$betavar==7 & alld$type=="qs_error_prob",]
qsdata0$qsdiff<-qsdata0$value-basepick_error_prob
qsdata5$qsdiff<-qsdata5$value-basepick_error_prob
qsdata6$qsdiff<-qsdata6$value-basepick_error_prob
qsdata7$qsdiff<-qsdata7$value-basepick_error_prob

rbind(qsdata0,qsdata5,qsdata6,qsdata7)->qsdatas
rbind(gldata0,gldata5,gldata6,gldata7)->gldatas




# plot the difference between the error probability used in GL calculation and the error probability used in base picking per beta variance
gldatas%>%
filter(precise==0)%>%
mutate(betavar=factor(betavar,levels=c(0,7,6,5)))%>%
# filter(betavar!=0)%>%
ggplot(aes(y=gldiff,color=betavar))+
# facet_wrap(.~precise,scale="fixed",labeller=labeller1fn)+
# geom_violin()
geom_boxplot()+
theme_bw()+
theme(strip.background = element_rect(colour="black", fill="white"))+
labs(
  x="Beta variance",
  y="Change in the value of error probability used in GL calculation",
  color="Beta variance"
)+
color_betavars

save_plt(plt=last_plot(),plot_filename="plot_boxplot_gl_error_prob.png",plot_outdir=plot_outdir,overwrite=TRUE,width=10,height=10,scale=1.5,dpi=300,units="cm")


gldatas%>%
filter(precise==0)%>%
mutate(betavar=factor(betavar,levels=c(0,7,6,5)))%>%
# filter(betavar!=0)%>%
ggplot(aes(x=gldiff,color=betavar))+
facet_wrap(.~betavar,scale="fixed",labeller=labeller1fn)+
# geom_violin()
# geom_boxplot()+
stat_bin(alpha=0.7,position="identity")+
# geom_col(position="identity", alpha=0.5)+
theme_bw()+
theme(strip.background = element_rect(colour="black", fill="white"))+
labs(
  x="Beta variance",
  y="Change in the value of error probability used in GL calculation",
  color="Beta variance"
)+
color_betavars

#mean change per beta variance
gldatas%>%
filter(precise==0)%>%
mutate(betavar=factor(betavar,levels=c(0,7,6,5)))%>%
group_by(betavar)%>%
summarise(mean=mean(abs(gldiff)))->observed_mean_gl_error_prob


gldatas%>%
filter(precise==0)%>%
# filter(betavar==7)%>%
# mutate(betavar=factor(betavar,levels=c(0,7,6,5)))%>%
filter(betavar!=0)%>%
ggplot(aes(x=gldiff,color=betavar,fill=betavar))+
facet_wrap(.~betavar,scale="free",labeller=labeller1fn)+
# geom_histogram(aes(y=after_stat(count/sum(count))), position="identity",binwidth=0.0001)+
geom_density(aes(y=after_stat(count/sum(count))), position="identity",alpha=0.7)+
# geom_violin()
# geom_boxplot()+
# stat_bin(alpha=0.7,position="identity")+
# geom_col(position="identity", alpha=0.5)+
theme_bw()+
theme(strip.background = element_rect(colour="black", fill="white"))+
labs(
  x="Change in the value of error probability used in GL calculation",
  y="Density",
  fill="Beta variance",
  color=""
)+
color_betavars+
fill_betavars+
theme(legend.position="top")+
guides(color=FALSE)


save_plt(plt=last_plot(),plot_filename="plot_density_gl_error_prob.png",plot_outdir=plot_outdir,overwrite=TRUE,width=10,height=10,scale=1.5,dpi=300,units="cm")



alld %>% 
filter(betavar!=0)%>%
filter(type=="qs")%>%
filter(precise==0)%>%
ggplot(aes(x=value,fill=betavar))+
stat_count(alpha=0.7,binwidth=1)+
facet_wrap(.~betavar,scale="free",labeller=labeller1fn)+
geom_vline(xintercept=true_qScore,color="black",size=1)+
theme_bw()+
theme(strip.background = element_rect(colour="black", fill="white"))+
labs(
  x="Q-score",
  y="Count",
  subtitle=paste0("Base picking error probability = ",round(basepick_error_prob,5),"\nTrue Q-score = ",true_qScore)
)+
fill_betavars


save_plt(plt=last_plot(),plot_filename="plot_histogram_qs.png",plot_outdir=plot_outdir,overwrite=TRUE,width=10,height=10,scale=1.5,dpi=300,units="cm")



