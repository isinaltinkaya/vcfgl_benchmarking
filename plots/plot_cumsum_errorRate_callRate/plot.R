###############################################################################
# isinaltinkaya
#

rm(list=ls())
gc()

################################################################################
# LOAD LIBRARIES

# require(data.table) 
# require(readr)

# library(gridExtra)
library(ggplot2)
# library(ggpubr)
library(dplyr)
# library(tidyr) # pivot_longer

today<-format(Sys.Date(), "%y%m%d")


# setwd(paste0("plots_",today))
plot_outdir<-getwd()



################################################################################
# FUNCTIONS

assert <- function(x) {
  if (!x) stop("Assertion failed.")
}

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


parse_betavar <- function(x){
  for(i in seq_along(levels(x))){
    if(levels(x)[i] != 0){
      levels(x)[i]<-gsub("X","10^-",paste0("X",levels(x)[i]))
    }
  }
  x
}

gcmethod_lut<-c("Genotype calling across populations","Genotype calling within populations")
names(gcmethod_lut)<-c("genotype_calling","genotype_calling_perpop")

glmethod_lut<-c("GL method 1","GL method 2","GL method 2 \nwith precise error")
names(glmethod_lut)<-c("gl1","gl2","precise1_gl2")

betavar_lut<-c("0","10^-7","10^-6","10^-5")
names(betavar_lut)<-c("0","7","6","5")


labeller1fn <- function(variable,value){
  if (variable=='Betavar') {
      return(label_parsed(ifelse(value==0,"0",paste0("10^-",value))))
  }else if (variable=='Depth'){
    return(as.character(paste0("Depth: ",as.character(value))))
  }else if (variable=='GcMethod'){
    return(as.character(gcmethod_lut[as.character(value)]))
  }else if (variable=='Gl'){
    return(as.character(glmethod_lut[as.character(value)]))
  }else{
    return(as.character(value))
  }
}


glmethod_lut2<-c("GL1","GL2","P_GL2")
names(glmethod_lut2)<-c("gl1","gl2","precise1_gl2")

labeller2fn <- function(variable,value){
  if (variable=='Betavar') {
      return(label_parsed(ifelse(value==0,"0",paste0("10^-",value))))
  }else if (variable=='Depth'){
    return(as.character(paste0("Depth: ",as.character(value))))
  }else if (variable=='GcMethod'){
    return(as.character(gcmethod_lut[as.character(value)]))
  }else if (variable=='Gl'){
    return(as.character(glmethod_lut2[as.character(value)]))
  }else{
    return(as.character(value))
  }
}



fill_trueType<-scale_fill_manual(values=c(trueHetDiscordanceRateOverAll="blue",trueHomDiscordanceRateOverAll="red"),labels=c(trueHetDiscordanceRateOverAll="True heterozygous",trueHomDiscordanceRateOverAll="True homozygous"))
color_trueType<-scale_color_manual(values=c(trueHetDiscordanceRateOverAll="blue",trueHomDiscordanceRateOverAll="red"),labels=c(trueHetDiscordanceRateOverAll="True heterozygous",trueHomDiscordanceRateOverAll="True homozygous"))
color_brewer_betavars<-scale_color_brewer(palette="Set1",labels=c("0"="0","5"=expression(10^-5),"6"=expression(10^-6),"7"=expression(10^-7)))
fill_brewer_betavars<-scale_fill_brewer(palette="Set1",labels=c("0"="0","5"=expression(10^-5),"6"=expression(10^-6),"7"=expression(10^-7)))
linetype_gls<-scale_linetype_manual(values=c("gl1"="solid","gl2"="dashed","precise1_gl2"="dotted"),labels=c("gl1"="GL method 1","gl2"="GL method 2","precise1_gl2"="GL method 2 with precise error"))
fill_brewer_typeToType<-scale_fill_brewer(palette="Set1",labels=c(HetToHetDiscordanceRateOverAll="Het to het",HomToHomDiscordanceRateOverAll="Hom to hom",HomToHetDiscordanceRateOverAll="Hom to het",HetToHomDiscordanceRateOverAll="Het to hom"))

################################################################################
# READ DATA



simdir="../../simulation/sim/"

Method_i=c("genotype_calling", "genotype_calling_perpop")
GlMethod_i=c("gl1","gl2","precise1_gl2")
hdr<-c("Sample","GQ","nDiscordant","nHomToHomDiscordant","nHomToHetDiscordant","nHetToHomDiscordant","nHetToHetDiscordant","nConcordant","nHomToHomConcordant","nHetToHetConcordant","nSitesComparedForSample","Rep","Depth","Betavar")

# sd<-read.csv(paste0(simdir,"sim_vcfgl_2312/model_OutOfAfrica_3G09/stats/nSites_non0dp/sim_vcfgl_2312-OutOfAfrica_3G09_nSitesNon0DP.csv"),sep=",",header=TRUE)
# # convert sample names to sample indices (0-indexed)
# sd$Sample<-as.numeric(factor(sd$Sample,levels=unique(sd$Sample)))-1

# d<-NULL
# for(mi in seq_along(Method_i)){
#   for(gli in seq_along(GlMethod_i)){

#   di<-read.csv(paste0(simdir, "sim_vcfgl_2312/model_OutOfAfrica_3G09/gc_evaluation/genotype_discordance_",GlMethod_i[gli],"/sim_vcfgl_2312-OutOfAfrica_3G09-",Method_i[mi],"-doGQ7.tsv" ),sep="\t",header=FALSE)

#   colnames(di)<-hdr
#   di$GcMethod<-Method_i[mi]
#   di$Gl<-GlMethod_i[gli]
#   di$nSitesGQ<-di$nDiscordant+di$nConcordant
#   di<-merge(di,sd,by=c("Sample","Rep","Depth"))
#   d<-rbind(d,di)
# }
# }




# save(d,file=paste0("./data_",today,".RData"))


load(file=paste0("./data_",today,".RData"))


################################################################################



rm(d1,d2,d3,plt);gc();

d1<-d[d$GcMethod=="genotype_calling_perpop",]



(d1%>% 
group_by(Depth, Betavar, Gl,GcMethod,GQ)%>%
summarize(nDiscordant= sum(nDiscordant),
nConcordant= sum(nConcordant),
nSitesGQ= sum(nSitesGQ),
nSitesNon0DP= sum(nSitesNon0DP))%>%
ungroup()%>% 
arrange(desc(GQ)) %>%
group_by(Depth, Betavar, Gl,GcMethod)%>%
mutate(DiscordanceRate= nDiscordant/(nDiscordant+nConcordant),
      ratioCallsTotal=cumsum(as.numeric(nSitesGQ))/nSitesNon0DP,
      min_GQ_threshold= GQ,
      ratioCallsIncludedSoFar=(cumsum(as.numeric(nDiscordant))+cumsum(as.numeric(nConcordant)))/(sum(nDiscordant)+sum(nConcordant)),
      DiscordanceRateOverCumsums= cumsum(as.numeric(nDiscordant))/(cumsum(as.numeric(nDiscordant))+cumsum(as.numeric(nConcordant))),
      DiscordanceRateOverCumsums2= cumsum(as.numeric(nDiscordant))/(1+(cumsum(as.numeric(nDiscordant))+cumsum(as.numeric(nConcordant)))),
      )%>% 
      ungroup()%>%
mutate(across(c('Gl','GcMethod','Depth','Betavar'),as.factor))%>%
mutate(Betavar=factor(Betavar,levels=c("0","7","6","5")))%>%
ungroup())->d2



rm(plt)
(plt<-ggplot(d2,aes(x=ratioCallsTotal,y=DiscordanceRateOverCumsums, color=Depth ,fill=Depth,linetype=Gl))+
facet_wrap(.~Betavar,scale="fixed",labeller=labeller1fn)+
geom_point(data=d2[d2$min_GQ_threshold==20,],aes(x = ratioCallsIncludedSoFar),size=3)+
geom_step()+
# geom_line()+
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
theme_bw()+
theme(legend.position = "top")+
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
  scale_color_brewer(palette="RdYlGn")+
  scale_fill_brewer(palette="RdYlGn")+
labs(
    x="Call rate",
    y="Error rate",
    color="Depth",
    linetype=""
    )+
theme(aspect.ratio=1)+
guides(fill=FALSE)+
guides(shape=FALSE)+
theme(strip.background = element_rect(colour="black", fill="white"))+
linetype_gls+
guides(linetype=guide_legend(nrow=3)))

plt


save_plt(plt,plot_filename="plot-step_x-ratioCallsTotal_y-DiscordanceRateOverCumsums_grid-Betavar_color-Depth_linetype-Gl_data-GenotypeCallingPerPop.png",plot_outdir=plot_outdir,overwrite=TRUE,width=10,height=10,units="in",scale=0.6,dpi=300)


rm(plt)
(plt<-ggplot(d2,aes(x=ratioCallsTotal,y=DiscordanceRateOverCumsums, color=Depth ,fill=Depth,linetype=Gl))+
facet_wrap(.~Betavar,scale="fixed",labeller=labeller1fn)+
geom_point(data=d2[d2$min_GQ_threshold==20,],aes(x = ratioCallsIncludedSoFar),size=3)+
# geom_step()+
geom_line()+
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
theme_bw()+
theme(legend.position = "top")+
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
  scale_color_brewer(palette="RdYlGn")+
  scale_fill_brewer(palette="RdYlGn")+
labs(
    x="Call rate",
    y="Error rate",
    color="Depth",
    linetype=""
    )+
theme(aspect.ratio=1)+
guides(fill=FALSE)+
guides(shape=FALSE)+
theme(strip.background = element_rect(colour="black", fill="white"))+
linetype_gls+
guides(linetype=guide_legend(nrow=3)))

plt


save_plt(plt,plot_filename="plot-line_x-ratioCallsTotal_y-DiscordanceRateOverCumsums_grid-Betavar_color-Depth_linetype-Gl_data-GenotypeCallingPerPop.png",plot_outdir=plot_outdir,overwrite=TRUE,width=10,height=10,units="in",scale=0.6,dpi=300)




################################################################################


d1<-d[d$GcMethod=="genotype_calling_perpop",]

d1$DiscordanceRate= d1$nDiscordant/(d1$nDiscordant+d1$nConcordant)


(d1%>% 
group_by(Sample,Rep,Depth, Betavar, Gl,GcMethod)%>%
summarize(
  nDiscordant= sum(nDiscordant),
  nConcordant= sum(nConcordant),
  nHomToHomDiscordant= sum(nHomToHomDiscordant),
  nHomToHetDiscordant= sum(nHomToHetDiscordant),
  nHetToHomDiscordant= sum(nHetToHomDiscordant),
  nHetToHetDiscordant= sum(nHetToHetDiscordant),
  nHomToHomConcordant= sum(nHomToHomConcordant),
  nHetToHetConcordant= sum(nHetToHetConcordant),
  nSitesNon0DP= mean(nSitesNon0DP),
  nSitesGQ= sum(nSitesGQ),
  nSitesComparedForSample= sum(nSitesComparedForSample),
)%>% 
ungroup()%>%
mutate(
  DiscordanceRate= nDiscordant/(nDiscordant+nConcordant),
)%>% 
ungroup()%>%
group_by(Depth, Betavar, Gl,GcMethod)%>%
summarize(
  DiscordanceRate= mean(DiscordanceRate,na.rm=TRUE),
)%>%
ungroup()%>%
mutate(
  across(c('Gl','GcMethod','Depth','Betavar'),as.factor)
)) -> d2


(d1%>% 
group_by(Depth, Betavar, Gl,GcMethod,GQ)%>%
summarize(
  nDiscordant= sum(nDiscordant),
  nConcordant= sum(nConcordant),
  nHomToHomDiscordant= sum(nHomToHomDiscordant),
  nHomToHetDiscordant= sum(nHomToHetDiscordant),
  nHetToHomDiscordant= sum(nHetToHomDiscordant),
  nHetToHetDiscordant= sum(nHetToHetDiscordant),
  nHomToHomConcordant= sum(nHomToHomConcordant),
  nHetToHetConcordant= sum(nHetToHetConcordant),
  nSitesNon0DP= mean(nSitesNon0DP),
  nSitesGQ= sum(nSitesGQ),
  nSitesComparedForSample= sum(nSitesComparedForSample),
)%>% 
ungroup()%>%
group_by(Depth, Betavar, Gl,GcMethod,GQ)%>%
mutate(
  DiscordanceRate= nDiscordant/(nDiscordant+nConcordant),
)%>% 
ungroup()%>%
mutate(
  across(c('Gl','GcMethod','Depth','Betavar'),as.factor)
)) -> d2



d2%>%
filter("genotype_calling_perpop"==GcMethod)%>%
group_by(Depth, Betavar, Gl,GcMethod)%>%
summarize(sumnDiscordant=sum(nDiscordant),sumnConcordant=sum(nConcordant),meanDiscordanceRate=mean(sum(nDiscordant)/(sum(nDiscordant)+sum(nConcordant))))%>% 
ungroup()%>%
mutate(
  across(c('Gl','GcMethod','Depth','Betavar'),as.factor)
)->d3

d3%>%
ggplot(aes(x=Depth,y=meanDiscordanceRate, color=Betavar,fill=Betavar))+
facet_grid(.~Gl,scale="free",labeller=labeller1fn)+
geom_col(position="dodge2")+
theme_bw()+
theme(strip.background = element_rect(colour="black", fill="white"))


d2%>%
filter("genotype_calling_perpop"==GcMethod)%>%
group_by(Depth, Betavar, Gl,GcMethod)%>%
summarize(
  nDiscordant= sum(nDiscordant),
  nConcordant= sum(nConcordant),
  nHomToHomDiscordant= sum(nHomToHomDiscordant),
  nHomToHetDiscordant= sum(nHomToHetDiscordant),
  nHetToHomDiscordant= sum(nHetToHomDiscordant),
  nHetToHetDiscordant= sum(nHetToHetDiscordant),
  nHomToHomConcordant= sum(nHomToHomConcordant),
  nHetToHetConcordant= sum(nHetToHetConcordant),
  # nSitesNon0DP= mean(nSitesNon0DP),
  nSitesGQ= sum(nSitesGQ),
  nSitesComparedForSample= sum(nSitesComparedForSample),
)%>%
ungroup()%>%
mutate(DiscordanceRate= nDiscordant/(nDiscordant+nConcordant))%>%
ungroup()%>%
group_by(Depth, Betavar, Gl,GcMethod)%>%
summarize(
  meanDiscordanceRate=mean(DiscordanceRate,na.rm=TRUE),
)%>%
ungroup()%>%
mutate(across(c('Gl','GcMethod','Depth','Betavar'),as.factor)
)->d4


d4%>%
ggplot(aes(x=Depth,y=meanDiscordanceRate, color=Betavar,fill=Betavar))+
facet_grid(GcMethod~Gl,scale="free",labeller=labeller1fn)+
geom_point()+
# geom_boxplot()+
theme_bw()+
theme(strip.background = element_rect(colour="black", fill="white"))


d2%>% 
# filter(Depth=="0.1")%>%
mutate(across(c('Gl','GcMethod','Depth','Betavar'),as.factor))%>%
ggplot(aes( x=GQ,y=DiscordanceRate, color=Betavar,fill = Betavar))+
  facet_grid(Gl+Betavar~Depth,scale="free",labeller=labeller1fn)+
  geom_col(position="identity",alpha=0.5)+
  # geom_density()
  # geom_point()+
  # geom_line()+
  theme_bw()



d1$Sample<-NULL
d1$Rep<-NULL


(d1%>% 
group_by(Depth, Betavar, Gl,GcMethod,GQ)%>%
summarize(nDiscordant= sum(nDiscordant),
nConcordant= sum(nConcordant),
nHomToHomDiscordant= sum(nHomToHomDiscordant),
nHomToHetDiscordant= sum(nHomToHetDiscordant),
nHetToHomDiscordant= sum(nHetToHomDiscordant),
nHetToHetDiscordant= sum(nHetToHetDiscordant),
nHomToHomConcordant= sum(nHomToHomConcordant),
nHetToHetConcordant= sum(nHetToHetConcordant),
# nSitesNon0DP= mean(nSitesNon0DP),
)%>%
ungroup()%>% 
mutate(nSitesGQ= nDiscordant+nConcordant)%>%
arrange(desc(GQ)) %>%
group_by(Depth, Betavar, Gl,GcMethod)%>%
mutate(
      # ratioCallsTotal=cumsum(as.numeric(nSitesGQ))/nSitesNon0DP,
      min_GQ_threshold= GQ,
      ratioCallsIncludedSoFar=(cumsum(as.numeric(nDiscordant))+cumsum(as.numeric(nConcordant)))/(sum(nDiscordant)+sum(nConcordant)),
      DiscordanceRateOverCumsums= cumsum(as.numeric(nDiscordant))/(cumsum(as.numeric(nDiscordant))+cumsum(as.numeric(nConcordant))),
      trueHomDiscordanceRateOverCumsums= cumsum(as.numeric(nHomToHomDiscordant+nHomToHetDiscordant))/(cumsum(as.numeric(nHomToHomDiscordant+nHomToHetDiscordant+nHomToHomConcordant))),
      trueHetDiscordanceRateOverCumsums= cumsum(as.numeric(nHetToHomDiscordant+nHetToHetDiscordant))/(cumsum(as.numeric(nHetToHomDiscordant+nHetToHetDiscordant+nHetToHetConcordant))),
      )%>% 
    ungroup()%>%
    mutate(across(c('Gl','GcMethod','Depth','Betavar'),as.factor)) %>%
ungroup())->d2



# # ################################################################################


# rm(d1,d2,d3,plt);gc();

# d1<-d[d$GcMethod=="genotype_calling_perpop",]


# d1%>%
# group_by(Depth, Betavar, Gl,GcMethod,GQ)%>%
# summarize(nDiscordant= sum(nDiscordant),
# nConcordant= sum(nConcordant))%>%
# ungroup()%>%
# mutate(DiscordanceRate= nDiscordant/(nDiscordant+nConcordant))%>%
# mutate(across(c('Gl','GcMethod','Depth','Betavar'),as.factor))->d2



# d2%>%
# ggplot(aes(x=GQ,y=DiscordanceRate))+
# facet_grid(Gl+Betavar~Depth, labeller=labeller1fn)+
# ylim(0,1)+
# geom_col()+
# theme_bw()+
# theme(strip.background = element_rect(colour="black", fill="white"))+
# #set theme text size
# theme(text = element_text(size=8))

# save_plt(plt= last_plot(),plot_filename="plot_col_x-GQ_y-DiscordanceRate_grid-GlBetavarDepth_gcPerPop.png",plot_outdir=plot_outdir,overwrite=TRUE,width=10,height=10,units="in",scale=0.6,dpi=300)

# # # # ####

# beta0<-d2%>%filter(Betavar==0)
# beta5<-d2%>%filter(Betavar==5)


# # # # beta0$nDiscordant0Minus5<-beta0$nDiscordant-beta5$nDiscordant
# # # # beta0$Diff0Minus5<-beta0$DiscordanceRate-beta5$DiscordanceRate


# # # # beta0$Betavar<-NULL
# # # # beta0%>%
# # # # ggplot(aes(x=GQ,y=Diff0Minus5))+
# # # # facet_grid(Gl~Depth, labeller=labeller1fn)+
# # # # ylim(min(beta0$Diff0Minus5),max(beta0$Diff0Minus5))+
# # # # xlim(1,129)+
# # # # # ylim(0,1)+
# # # # ylim(
# # # #   -max(abs(beta0$Diff0Minus5),na.rm=TRUE),
# # # #   max(abs(beta0$Diff0Minus5),na.rm=TRUE))+
# # # # geom_line()+
# # # # theme_bw()+
# # # # #add horizontal 0 line
# # # # geom_hline(yintercept=0,color="red")+
# # # # geom_text(data=data.frame(x=70,y=-0.10),aes(x=x,y=y,label="D(0) < D(10^-5)"),parse=T,color="red")+
# # # # geom_text(data=data.frame(x=70,y=+0.10),aes(x=x,y=y,label="D(10^-5) < D(0)"),parse=T,color="red")+
# # # # theme(strip.background = element_rect(colour="black", fill="white"))+
# # # # labs(
# # # #     x="GQ",
# # # #     y="Difference between discordance rates (1e-5 - 0)",
# # # #     color="Beta variance",
# # # #     linetype=""
# # # #     )
  

# # # # save_plt(plt= last_plot(), plot_filename="plot_x-GQ_y-Diff0Minus5_grid-BetavarDepth_gcPerPop.png",plot_outdir=plot_outdir,overwrite=TRUE,width=12,height=8,units="in",scale=0.9,dpi=300)


# # # # # 



# # # # ####


# # # # rm(d1,d2,d3,plt);gc();

# # # # d1<-d[d$GcMethod=="genotype_calling_perpop",]
# # # # d1<-d1[d1$Gl %in% c("gl2","precise1_gl2"),]




# # # # ggplot(d2,aes(x=min_GQ_threshold,y=DiscordanceRateOverCumsums, color=Betavar ,linetype=Gl))+
# # # # facet_wrap(.~Depth,scale="free",labeller=labeller1fn)+
# # # # geom_line()+
# # # # geom_vline(xintercept=20,color="black")+
# # # # theme_bw()+
# # # # theme(strip.background = element_rect(colour="black", fill="white"))

# # # # ggsave("plot_x-min_GQ_threshold_y-DiscordanceRateOverCumsums_grid-Depth_color-Betavar_linetype-Gl.png",width=10,height=10,units="in",scale=0.6,dpi=300)



# # # # #   Betavar overallDiscordanceRate
# # # # #   <fct>                    <dbl>
# # # # # 1 0                       0.0246
# # # # # 2 5                       0.0263
# # # # # 3 6                       0.0239


# # # # d%>%
# # # # filter(GcMethod=="genotype_calling_perpop")%>%
# # # # filter(Gl=="precise1_gl2")%>% 
# # # # group_by(Depth, Betavar, Gl,GcMethod,GQ)%>%
# # # # summarize(nDiscordant= sum(nDiscordant),
# # # # nConcordant= sum(nConcordant))%>%
# # # # ungroup()%>%
# # # # mutate(DiscordanceRate= nDiscordant/(nDiscordant+nConcordant))%>%
# # # # group_by(Betavar,Gl,GcMethod,Depth)%>% 
# # # # summarize(overallDiscordanceRate=mean(DiscordanceRate,na.rm=TRUE))%>%
# # # # ungroup()%>%
# # # # mutate(across(c('Gl','GcMethod','Depth','Betavar'),as.factor)) %>%
# # # # ggplot(aes(x=Depth,y=overallDiscordanceRate,color=Betavar,fill=Betavar))+
# # # # facet_grid(.~Gl,scale="fixed")+
# # # # geom_col(position="dodge2")+
# # # # theme_bw()+
# # # # theme(strip.background = element_rect(colour="black", fill="white"))+
# # # # theme(legend.position = "top")





# # # # #    Betavar Gl           GcMethod                Depth overallDiscordanceRate
# # # # #    <fct>   <fct>        <fct>                   <fct>                  <dbl>
# # # # #  1 0       gl2          genotype_calling_perpop 0.1                   0.0770
# # # # #  2 0       gl2          genotype_calling_perpop 0.5                   0.0342
# # # # #  3 0       gl2          genotype_calling_perpop 1                     0.0254
# # # # #  4 0       gl2          genotype_calling_perpop 2                     0.0213
# # # # #  5 0       gl2          genotype_calling_perpop 10                    0.0193
# # # # #  6 0       gl2          genotype_calling_perpop 20                    0.0177
# # # # #  7 0       precise1_gl2 genotype_calling_perpop 0.1                   0.0769
# # # # #  8 0       precise1_gl2 genotype_calling_perpop 0.5                   0.0329
# # # # #  9 0       precise1_gl2 genotype_calling_perpop 1                     0.0246
# # # # # 10 0       precise1_gl2 genotype_calling_perpop 2                     0.0214

# # # # # colnames(d)

# # # # d2%>%
# # # # group_by(Betavar,Gl,GcMethod,Depth,GQ)%>%
# # # # summarize(DiscordanceRate=sum(nDiscordant)/(sum(nDiscordant)+sum(nConcordant))) %>% 
# # # # ungroup()%>%
# # # # ggplot(aes(x=GQ,y=DiscordanceRate,color=Betavar,fill=Betavar))+
# # # # facet_wrap(Gl~Depth,scale="free",labeller=labeller1fn)+
# # # # geom_col(position="fill")+
# # # # theme_bw()+
# # # # theme(strip.background = element_rect(colour="black", fill="white"))+
# # # # labs(
# # # #   x="GQ",
# # # #   y="Discordance rate"
# # # # )+
# # # # #change the beta variance labels, 5->1e-5, 6->1e-6, 0->0
# # # # scale_color_brewer(palette="Set1",labels=c("0"="0","5"=expression(10^-5),"6"=expression(10^-6)))+
# # # # scale_fill_brewer(palette="Set1",labels=c("0"="0","5"=expression(10^-5),"6"=expression(10^-6)))

# # # # ggsave("plot_col_x-GQ_y-DiscordanceRate_grid-GlDepth_gcPerPop.png",width=10,height=10,units="in",scale=0.6,dpi=300)




# # # # (plt<-ggplot(d2,aes(x=ratioCallsIncludedSoFar,y=DiscordanceRateOverCumsums, color=Betavar ,linetype=Gl))+
# # # # facet_wrap(.~Depth,scale="free",labeller=labeller1fn)+
# # # # # geom_point()+
# # # # geom_point(data=d2[d2$min_GQ_threshold==20,],aes(x = ratioCallsIncludedSoFar),size=3)+
# # # # geom_line()+
# # # # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
# # # # theme_bw()+
# # # # theme(legend.position = "top")+
# # # #   theme(text = element_text(size=10))+
# # # #   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
# # # #   # scale_color_brewer(palette="RdYlGn")+
# # # #   # scale_fill_brewer(palette="RdYlGn")+
# # # # labs(
# # # #     x="Call rate (Over total genotype calls)",
# # # #     y="Error rate",
# # # #     color="Beta variance",
# # # #     linetype=""
# # # #     )+
# # # # theme(aspect.ratio=1)+
# # # # guides(fill=FALSE)+
# # # # guides(shape=FALSE)+
# # # # theme(strip.background = element_rect(colour="black", fill="white"))+
# # # # scale_color_brewer(palette="Set1",labels=c("0"="0","5"=expression(10^-5),"6"=expression(10^-6)))+
# # # # scale_fill_brewer(palette="Set1",labels=c("0"="0","5"=expression(10^-5),"6"=expression(10^-6)))+
# # # # guides(linetype=guide_legend(nrow=3)))


# # # # save_plt(plt,plot_filename="plot_x-ratioCallsIncludedSoFar_y-DiscordanceRateOverCumsums_grid-Depth_color-Betavar_linetype-Gl.png",plot_outdir=plot_outdir,overwrite=TRUE,width=10,height=10,units="in",scale=0.6,dpi=300)


# # # # # (plt<-ggplot(d2,aes(x=ratioCallsTotal,y=DiscordanceRateOverCumsums, color=Betavar ,linetype=Gl))+
# # # # # facet_wrap(.~Depth,scale="free",labeller=labeller1fn)+
# # # # # # geom_point()+
# # # # # geom_point(data=d2[d2$min_GQ_threshold==20,],aes(x = ratioCallsTotal),size=3)+
# # # # # geom_line()+
# # # # # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
# # # # # theme_bw()+
# # # # # theme(legend.position = "top")+
# # # # #   theme(text = element_text(size=10))+
# # # # #   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
# # # # #   # scale_color_brewer(palette="RdYlGn")+
# # # # #   # scale_fill_brewer(palette="RdYlGn")+
# # # # # labs(
# # # # #     x="Call rate (Over total non-zero depth sites)",
# # # # #     y="Error rate",
# # # # #     color="Beta variance",
# # # # #     linetype=""
# # # # #     )+
# # # # # theme(aspect.ratio=1)+
# # # # # guides(fill=FALSE)+
# # # # # guides(shape=FALSE)+
# # # # # theme(strip.background = element_rect(colour="black", fill="white"))+
# # # # # scale_linetype_manual(values=c("solid","dotted"),labels=c("GL method 2","GL method 2 with precise error"))+
# # # # # scale_color_manual(values=c("black","blue","red"),labels=c("0","10^-6","10^-5"))+
# # # # # guides(linetype=guide_legend(nrow=3)))


# # # # # save_plt(plt,plot_filename="plot_x-ratioCallsTotal_y-DiscordanceRateOverCumsums_grid-Depth_color-Betavar_linetype-Gl.png",plot_outdir=plot_outdir,overwrite=TRUE,width=10,height=10,units="in",scale=0.6,dpi=300)







# # # # rm(d1,d2,d3,plt);gc();

# # # # d1<-d[d$GcMethod=="genotype_calling_perpop",]
# # # # d1<-d1[d1$Gl %in% c("gl2","precise1_gl2"),]

# d1$Sample<-NULL
# d1$Rep<-NULL
# (d1%>% 
# group_by(Depth, Betavar, Gl,GcMethod,GQ)%>%
# summarize(nDiscordant= sum(nDiscordant),
# nConcordant= sum(nConcordant),
# nHomToHomDiscordant= sum(nHomToHomDiscordant),
# nHomToHetDiscordant= sum(nHomToHetDiscordant),
# nHetToHomDiscordant= sum(nHetToHomDiscordant),
# nHetToHetDiscordant= sum(nHetToHetDiscordant),
# nHomToHomConcordant= sum(nHomToHomConcordant),
# nHetToHetConcordant= sum(nHetToHetConcordant),
# # nSitesNon0DP= mean(nSitesNon0DP),
# )%>%
# ungroup()%>% 
# mutate(nSitesGQ= nDiscordant+nConcordant)%>%
# arrange(desc(GQ)) %>%
# group_by(Depth, Betavar, Gl,GcMethod)%>%
# mutate(
#       # ratioCallsTotal=cumsum(as.numeric(nSitesGQ))/nSitesNon0DP,
#       min_GQ_threshold= GQ,
#       ratioCallsIncludedSoFar=(cumsum(as.numeric(nDiscordant))+cumsum(as.numeric(nConcordant)))/(sum(nDiscordant)+sum(nConcordant)),
#       DiscordanceRateOverCumsums= cumsum(as.numeric(nDiscordant))/(cumsum(as.numeric(nDiscordant))+cumsum(as.numeric(nConcordant))),
#       trueHomDiscordanceRateOverCumsums= cumsum(as.numeric(nHomToHomDiscordant+nHomToHetDiscordant))/(cumsum(as.numeric(nHomToHomDiscordant+nHomToHetDiscordant+nHomToHomConcordant))),
#       trueHetDiscordanceRateOverCumsums= cumsum(as.numeric(nHetToHomDiscordant+nHetToHetDiscordant))/(cumsum(as.numeric(nHetToHomDiscordant+nHetToHetDiscordant+nHetToHetConcordant))),
#       )%>% 
#     ungroup()%>%
#     mutate(across(c('Gl','GcMethod','Depth','Betavar'),as.factor)) %>%
# ungroup())->d2

# # # # (plt<-ggplot(d2,aes(x=ratioCallsIncludedSoFar,y=trueHomDiscordanceRateOverCumsums, color=Betavar ,linetype=Gl))+
# # # # facet_wrap(.~Depth,scale="free",labeller=labeller1fn)+
# # # # # geom_point()+
# # # # geom_point(data=d2[d2$min_GQ_threshold==20,],aes(x = ratioCallsIncludedSoFar),size=3)+
# # # # geom_line()+
# # # # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
# # # # theme_bw()+
# # # # theme(legend.position = "top")+
# # # #   theme(text = element_text(size=10))+
# # # #   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
# # # # labs(
# # # #     x="Call rate (Over total genotype calls)",
# # # #     y="Error rate",
# # # #     color="Beta variance",
# # # #     linetype=""
# # # #     )+
# # # # theme(aspect.ratio=1)+
# # # # guides(fill=FALSE)+
# # # # guides(shape=FALSE)+
# # # # theme(strip.background = element_rect(colour="black", fill="white"))+
# # # # scale_linetype_manual(values=c("gl2"="solid","precise1_gl2"="dotted"),labels=c("GL method 2","GL method 2 with precise error")) +
# # # # scale_color_manual(values=c("0"="black","5"="blue","6"="red"),labels=c("0"="0","5"=expression(10^-5),"6"=expression(10^-6)))+
# # # # guides(linetype=guide_legend(nrow=3)))


# # # # save_plt(plt,plot_filename="plot-line_x-ratioCallsIncludedSoFar_y-trueHomDiscordanceRateOverCumsums_grid-Depth_color-Betavar_linetype-Gl.png",plot_outdir=plot_outdir,overwrite=TRUE,width=10,height=10,units="in",scale=0.6,dpi=300)

# # # # (plt<-ggplot(d2,aes(x=ratioCallsIncludedSoFar,y=trueHetDiscordanceRateOverCumsums, color=Betavar ,linetype=Gl))+
# # # # facet_wrap(.~Depth,scale="free",labeller=labeller1fn)+
# # # # # geom_point()+
# # # # geom_point(data=d2[d2$min_GQ_threshold==20,],aes(x = ratioCallsIncludedSoFar),size=3)+
# # # # geom_line()+
# # # # ggtitle("True genotype: Heterozygous")+
# # # # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
# # # # theme_bw()+
# # # # theme(legend.position = "top")+
# # # #   theme(text = element_text(size=10))+
# # # #   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
# # # # labs(
# # # #     x="Call rate (Over total genotype calls)",
# # # #     y="Error rate",
# # # #     color="Beta variance",
# # # #     linetype=""
# # # #     )+
# # # # theme(aspect.ratio=1)+
# # # # guides(fill=FALSE)+
# # # # guides(shape=FALSE)+
# # # # theme(strip.background = element_rect(colour="black", fill="white"))+
# # # # scale_linetype_manual(values=c("gl2"="solid","precise1_gl2"="dotted"),labels=c("GL method 2","GL method 2 with precise error")) +
# # # # scale_color_manual(values=c("0"="black","5"="blue","6"="red"),labels=c("0"="0","5"=expression(10^-5),"6"=expression(10^-6)))+
# # # # guides(linetype=guide_legend(nrow=3)))


# # # # save_plt(plt,plot_filename="plot-line_x-ratioCallsIncludedSoFar_y-trueHetDiscordanceRateOverCumsums_grid-Depth_color-Betavar_linetype-Gl.png",plot_outdir=plot_outdir,overwrite=TRUE,width=10,height=10,units="in",scale=0.6,dpi=300)






# # # # d1<-d[d$GcMethod=="genotype_calling_perpop",]
# # # # d1<-d1[d1$Gl %in% c("gl2","precise1_gl2"),]

d1%>%
group_by(Depth, Betavar, Gl,GcMethod,GQ)%>%
summarize(nDiscordant= sum(nDiscordant),
nConcordant= sum(nConcordant),
nHomToHomDiscordant= sum(nHomToHomDiscordant),
nHomToHetDiscordant= sum(nHomToHetDiscordant),
nHetToHomDiscordant= sum(nHetToHomDiscordant),
nHetToHetDiscordant= sum(nHetToHetDiscordant),
nHomToHomConcordant= sum(nHomToHomConcordant),
nHetToHetConcordant= sum(nHetToHetConcordant),
# nSitesNon0DP= mean(nSitesNon0DP)
)%>%
ungroup()%>%
mutate(DiscordanceRate= nDiscordant/(nDiscordant+nConcordant),
trueHomDiscordanceRate= (nHomToHomDiscordant+nHomToHetDiscordant)/(nHomToHomDiscordant+nHomToHetDiscordant+nHomToHomConcordant),
trueHetDiscordanceRate= (nHetToHomDiscordant+nHetToHetDiscordant)/(nHetToHomDiscordant+nHetToHetDiscordant+nHetToHetConcordant),
HomToHomDiscordanceRate= nHomToHomDiscordant/(nHomToHomDiscordant+nHomToHomConcordant),
HetToHetDiscordanceRate= nHetToHetDiscordant/(nHetToHetDiscordant+nHetToHetConcordant))%>%
mutate(
Gl=as.factor(Gl),
GcMethod=as.factor(GcMethod),
Betavar=as.factor(Betavar),
Depth=as.factor(Depth))->d2



# melt the data so that new column is the type of discordance rate 
# (true, hom2hom, het2het):
d3<-d1%>%
group_by(Depth, Betavar, Gl,GcMethod,GQ)%>%
summarize(nDiscordant= sum(nDiscordant),
nConcordant= sum(nConcordant),
nHomToHomDiscordant= sum(nHomToHomDiscordant),
nHomToHetDiscordant= sum(nHomToHetDiscordant),
nHetToHomDiscordant= sum(nHetToHomDiscordant),
nHetToHetDiscordant= sum(nHetToHetDiscordant),
nHomToHomConcordant= sum(nHomToHomConcordant),
nHetToHetConcordant= sum(nHetToHetConcordant),
# nSitesNon0DP= mean(nSitesNon0DP)
)%>%
ungroup()%>%
mutate(DiscordanceRate= nDiscordant/(nDiscordant+nConcordant),
trueHomDiscordanceRate= (nHomToHomDiscordant+nHomToHetDiscordant)/(nHomToHomDiscordant+nHomToHetDiscordant+nHomToHomConcordant),
trueHetDiscordanceRate= (nHetToHomDiscordant+nHetToHetDiscordant)/(nHetToHomDiscordant+nHetToHetDiscordant+nHetToHetConcordant),
HomToHomDiscordanceRate= nHomToHomDiscordant/(nHomToHomDiscordant+nHomToHomConcordant),
HetToHetDiscordanceRate= nHetToHetDiscordant/(nHetToHetDiscordant+nHetToHetConcordant),
trueHomDiscordanceRateOverAll = (nHomToHomDiscordant+nHomToHetDiscordant)/(nDiscordant+nConcordant),
trueHetDiscordanceRateOverAll = (nHetToHomDiscordant+nHetToHetDiscordant)/(nDiscordant+nConcordant),
HomToHomDiscordanceRateOverAll = nHomToHomDiscordant/(nDiscordant+nConcordant),
HetToHetDiscordanceRateOverAll = nHetToHetDiscordant/(nDiscordant+nConcordant),
)%>%
mutate(
Gl=as.factor(Gl),
GcMethod=as.factor(GcMethod),
Betavar=as.factor(Betavar),
Depth=as.factor(Depth))%>%
gather(key, value, trueHomDiscordanceRateOverAll, trueHetDiscordanceRateOverAll)%>%
mutate(key=factor(key))

d3 %>%
ggplot(aes(x=GQ,y=value,color=key,fill=key))+
facet_grid(Gl+Betavar~Depth, labeller=labeller2fn)+
geom_col(position="fill")+
ylim(0,1)+
theme_bw()+
theme(strip.background = element_rect(colour="black", fill="white"))+
theme(legend.position = "top")+
theme(text = element_text(size=10))+
xlab("GQ")+
ylab("Discordance rate")+
labs(x="GQ",y="Discordance rate",
color="",
fill="")+
#rename key values with lookup
scale_fill_manual(values=c(trueHetDiscordanceRateOverAll="blue",trueHomDiscordanceRateOverAll="red"),labels=c(trueHetDiscordanceRateOverAll="True heterozygous",trueHomDiscordanceRateOverAll="True homozygous"))+
scale_color_manual(values=c(trueHetDiscordanceRateOverAll="blue",trueHomDiscordanceRateOverAll="red"),labels=c(trueHetDiscordanceRateOverAll="True heterozygous",trueHomDiscordanceRateOverAll="True homozygous"))

save_plt(plt= last_plot(),plot_filename="plot_colFill_x-GQ_y-value_grid-GlBetavarDepth_fill-trueHomtrueHetKey_gcPerPop.png",plot_outdir=plot_outdir,overwrite=TRUE,width=10,height=10,units="in",scale=0.8,dpi=300)






d3 %>%
ggplot(aes(x=GQ,y=value,color=key,fill=key))+
facet_grid(Gl+Betavar~Depth, labeller=labeller2fn)+
geom_col()+
ylim(0,1)+
theme_bw()+
theme(strip.background = element_rect(colour="black", fill="white"))+
theme(legend.position = "top")+
theme(text = element_text(size=10))+
xlab("GQ")+
ylab("Discordance rate")+
labs(x="GQ",y="Discordance rate",
color="",
fill="")+
#rename key values with lookup
scale_fill_manual(values=c(trueHetDiscordanceRateOverAll="blue",trueHomDiscordanceRateOverAll="red"),labels=c(trueHetDiscordanceRateOverAll="True heterozygous",trueHomDiscordanceRateOverAll="True homozygous"))+
scale_color_manual(values=c(trueHetDiscordanceRateOverAll="blue",trueHomDiscordanceRateOverAll="red"),labels=c(trueHetDiscordanceRateOverAll="True heterozygous",trueHomDiscordanceRateOverAll="True homozygous"))


save_plt(plt= last_plot(),plot_filename="plot_col_x-GQ_y-value_grid-GlBetavarDepth_fill-trueHomtrueHetKey_gcPerPop.png",plot_outdir=plot_outdir,overwrite=TRUE,width=10,height=10,units="in",scale=0.8,dpi=300)




d3<-d1%>%
group_by(Depth, Betavar, Gl,GcMethod,GQ)%>%
summarize(nDiscordant= sum(nDiscordant),
nConcordant= sum(nConcordant),
nHomToHomDiscordant= sum(nHomToHomDiscordant),
nHomToHetDiscordant= sum(nHomToHetDiscordant),
nHetToHomDiscordant= sum(nHetToHomDiscordant),
nHetToHetDiscordant= sum(nHetToHetDiscordant),
nHomToHomConcordant= sum(nHomToHomConcordant),
nHetToHetConcordant= sum(nHetToHetConcordant),
# nSitesNon0DP= mean(nSitesNon0DP)
)%>%
ungroup()%>%
mutate(DiscordanceRate= nDiscordant/(nDiscordant+nConcordant),
trueHomDiscordanceRate= (nHomToHomDiscordant+nHomToHetDiscordant)/(nHomToHomDiscordant+nHomToHetDiscordant+nHomToHomConcordant),
trueHetDiscordanceRate= (nHetToHomDiscordant+nHetToHetDiscordant)/(nHetToHomDiscordant+nHetToHetDiscordant+nHetToHetConcordant),
HomToHomDiscordanceRate= nHomToHomDiscordant/(nHomToHomDiscordant+nHomToHomConcordant),
HetToHetDiscordanceRate= nHetToHetDiscordant/(nHetToHetDiscordant+nHetToHetConcordant),
trueHomDiscordanceRateOverAll = (nHomToHomDiscordant+nHomToHetDiscordant)/(nDiscordant+nConcordant),
trueHetDiscordanceRateOverAll = (nHetToHomDiscordant+nHetToHetDiscordant)/(nDiscordant+nConcordant),
HomToHomDiscordanceRateOverAll = nHomToHomDiscordant/(nDiscordant+nConcordant),
HetToHetDiscordanceRateOverAll = nHetToHetDiscordant/(nDiscordant+nConcordant),
HetToHomDiscordanceRateOverAll = nHetToHomDiscordant/(nDiscordant+nConcordant),
HomToHetDiscordanceRateOverAll = nHomToHetDiscordant/(nDiscordant+nConcordant))%>%
mutate(
Gl=as.factor(Gl),
GcMethod=as.factor(GcMethod),
Betavar=as.factor(Betavar),
Depth=as.factor(Depth))%>%
gather(key, value, HomToHomDiscordanceRateOverAll, HetToHetDiscordanceRateOverAll, HomToHetDiscordanceRateOverAll, HetToHomDiscordanceRateOverAll)%>%
mutate(key=factor(key))

d3 %>%
ggplot(aes(x=GQ,y=value,color=key,fill=key))+
facet_grid(Gl+Betavar~Depth, labeller=labeller2fn)+
geom_col(position="fill")+
ylim(0,1)+
theme_bw()+
theme(strip.background = element_rect(colour="black", fill="white"))+
theme(legend.position = "top")+
theme(text = element_text(size=10))+
xlab("GQ")+
ylab("Discordance rate")+
labs(x="GQ",y="Discordance rate",
color="",
fill="")+
#rename key values with lookup
scale_fill_brewer(palette="Set1",labels=c(HetToHetDiscordanceRateOverAll="Het to het",HomToHomDiscordanceRateOverAll="Hom to hom",HomToHetDiscordanceRateOverAll="Hom to het",HetToHomDiscordanceRateOverAll="Het to hom"))+
scale_color_brewer(palette="Set1",labels=c(HetToHetDiscordanceRateOverAll="Het to het",HomToHomDiscordanceRateOverAll="Hom to hom",HomToHetDiscordanceRateOverAll="Hom to het",HetToHomDiscordanceRateOverAll="Het to hom"))

save_plt(plt= last_plot(),plot_filename="plot_col_x-GQ_y-value_grid-GlBetavarDepth_fill-HomHetKey_gcPerPop.png",plot_outdir=plot_outdir,overwrite=TRUE,width=10,height=10,units="in",scale=0.8,dpi=300)





d3<-d1%>%
group_by(Depth, Betavar, Gl,GcMethod,GQ)%>%
summarize(nDiscordant= sum(nDiscordant),
nConcordant= sum(nConcordant),
nHomToHomDiscordant= sum(nHomToHomDiscordant),
nHomToHetDiscordant= sum(nHomToHetDiscordant),
nHetToHomDiscordant= sum(nHetToHomDiscordant),
nHetToHetDiscordant= sum(nHetToHetDiscordant),
nHomToHomConcordant= sum(nHomToHomConcordant),
nHetToHetConcordant= sum(nHetToHetConcordant),
nSitesNon0DP= sum(nSitesNon0DP)
)%>%
ungroup()%>%
mutate(DiscordanceRate= nDiscordant/(nDiscordant+nConcordant),
trueHomDiscordanceRate= (nHomToHomDiscordant+nHomToHetDiscordant)/(nHomToHomDiscordant+nHomToHetDiscordant+nHomToHomConcordant),
trueHetDiscordanceRate= (nHetToHomDiscordant+nHetToHetDiscordant)/(nHetToHomDiscordant+nHetToHetDiscordant+nHetToHetConcordant),
HomToHomDiscordanceRate= nHomToHomDiscordant/(nHomToHomDiscordant+nHomToHomConcordant),
HetToHetDiscordanceRate= nHetToHetDiscordant/(nHetToHetDiscordant+nHetToHetConcordant),
trueHomDiscordanceRateOverAll = (nHomToHomDiscordant+nHomToHetDiscordant)/(nDiscordant+nConcordant),
trueHetDiscordanceRateOverAll = (nHetToHomDiscordant+nHetToHetDiscordant)/(nDiscordant+nConcordant),
HomToHomDiscordanceRateOverAll = nHomToHomDiscordant/(nDiscordant+nConcordant),
HetToHetDiscordanceRateOverAll = nHetToHetDiscordant/(nDiscordant+nConcordant),
HetToHomDiscordanceRateOverAll = nHetToHomDiscordant/(nDiscordant+nConcordant),
HomToHetDiscordanceRateOverAll = nHomToHetDiscordant/(nDiscordant+nConcordant))%>%
mutate(
  across(c('Gl','GcMethod','Depth','Betavar'),as.factor))

d3 %>%
ggplot(aes(x=GQ,y=trueHomDiscordanceRateOverAll,color=Betavar,fill=Betavar))+
facet_wrap(.~Depth+Gl,scale="free",labeller=labeller1fn)+
# geom_point(data=d2[d2$min_GQ_threshold==20,],aes(x = ratioCallsIncludedSoFar),size=3)+
# geom_step()+
geom_line()+
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
theme_bw()+
theme(legend.position = "top")+
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
  scale_color_brewer(palette="RdYlGn")+
  scale_fill_brewer(palette="RdYlGn")+
labs(
    # x="Call rate",
    # y="Error rate",
    # color="Depth",
    linetype=""
    )+
theme(aspect.ratio=1)+
guides(fill=FALSE)+
guides(shape=FALSE)+
theme(strip.background = element_rect(colour="black", fill="white"))+
# linetype_gls+
guides(linetype=guide_legend(nrow=3))+
color_brewer_betavars






# # # # ggsave("plot_x-ratioCallsIncludedSoFar_y-DiscordanceRateOverCumsums_grid-BetavarDepth_gcPerPop.png",width=10,height=10,units="in")

# # # ggplot(d2,aes(x=min_GQ_threshold,y=DiscordanceRateOverCumsums, color=Depth ,fill=Depth))+
# # # facet_grid(Gl~Betavar,scale="fixed")+
# # # # geom_point()+
# # # geom_line()+
# # # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
# # # theme_bw()+
# # # theme(legend.position = "top")+
# # #   theme(text = element_text(size=10))+
# # #   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
# # #   scale_color_brewer(palette="RdYlGn")+
# # #   scale_fill_brewer(palette="RdYlGn")+
# # # # make it square
# # # theme(aspect.ratio=1)

# # # ggsave("plot_x-GQ_y-DiscordanceRate_grid-BetavarDepth_gcPerPop.png",width=10,height=10,units="in")

# # # ggplot(d2,aes(x=min_GQ_threshold,y=DiscordanceRate, color=Depth ,fill=Depth))+
# # # facet_grid(Gl~Betavar,scale="fixed")+
# # # # geom_point()+
# # # geom_line()+
# # # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
# # # theme_bw()+
# # # theme(legend.position = "top")+
# # #   theme(text = element_text(size=10))+
# # #   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
# # #   scale_color_brewer(palette="RdYlGn")+
# # #   scale_fill_brewer(palette="RdYlGn")+
# # # # make it square
# # # theme(aspect.ratio=1)

# # # ggsave("plot_x-GQ_y-DiscordanceRatePerGQ_grid-BetavarDepth_gcPerPop.png",width=10,height=10,units="in")

# # # ggplot(d2,aes(x=GQ,y=DiscordanceRateOverCumsums, color=Depth ,fill=Depth))+
# # # facet_grid(Gl~Betavar+Depth,scale="fixed")+
# # # # geom_point()+
# # # geom_line()+
# # # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
# # # theme_bw()+
# # # theme(legend.position = "top")+
# # #   theme(text = element_text(size=10))+
# # #   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
# # #   scale_color_brewer(palette="RdYlGn")+
# # #   scale_fill_brewer(palette="RdYlGn")+
# # # # make it square
# # # theme(aspect.ratio=1)





# # # # x axis only 


# # # ggplot(d2,aes(x=GQ,y=DiscordanceRate, color=Depth ,fill=Depth,linetype=Gl))+
# # # facet_grid(Gl~Betavar)+
# # # # geom_col()+
# # # geom_point()+
# # # geom_line()+
# # # #turn x axis 90 degrees
# # # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
# # # theme_bw()+
# # # theme(legend.position = "top")+
# # #   theme(text = element_text(size=10))+
# # #   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
# # #   scale_color_brewer(palette="RdYlGn")+
# # #   scale_fill_brewer(palette="RdYlGn")+
# # # # make it square
# # # theme(aspect.ratio=1)
# # # # x axis only 




































# # # d2%>%
# # # ggplot(aes(x=min_GQ_threshold, y = DiscordanceRateOverCumsums, color = Depth, fill=Gl))+
# # # geom_line()


# # # rm(plt)
# # # # if DiscordanceRate is NA, remove the row
# # # d3<-d2[d2$Gl==2,]
# # # (plt<-(
# # # ggplot(data=d3,aes(x = ratioCallsTotal, y = DiscordanceRateOverCumsums, color = GcMethod, fill=GcMethod))+
# # #   geom_line(linewidth=1)+
# # #   geom_point(data=d3[d3$min_GQ_threshold==20,],aes(x = ratioCallsTotal),size=3)+
# # #   facet_wrap(.~Depth,labeller=labeller1fn,scale="free")+
# # #   xlim(0,1)+
# # #   labs(
# # #     x="Call rate",
# # #     y="Error rate",
# # #     color="",
# # #     linetype="Genotype likelihood method",
# # #     )+
# # #   scale_color_manual(labels=c("Genotype calling across populations","Genotype calling within populations"),values=c("red","blue"))+
# # #   guides(color = guide_legend(nrow= 2))))






# # # d2

# # # ggplot(d2,aes(x=GQ,y=DiscordanceRate, color=Depth ,fill=Depth,linetype=Gl))+
# # # facet_grid(Gl~Betavar)+
# # # # geom_col()+
# # # # geom_point()+
# # # geom_line()+
# # # #turn x axis 90 degrees
# # # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
# # # theme_bw()+
# # # theme(legend.position = "top")+
# # #   theme(text = element_text(size=10))+
# # #   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
# # #   scale_color_brewer(palette="RdYlGn")+
# # #   scale_fill_brewer(palette="RdYlGn")+
# # # # make it square
# # # theme(aspect.ratio=1)
# # # # x axis only 

# # # ggsave("plot_x-GQ_y-DiscordanceRate_grid-BetavarDepth_gcAllSamples.png",width=10,height=10,units="in")


# # # d1$ErrorRate=d1$nDiscordant/(d1$nDiscordant+d1$nConcordant)
# # # d1$CallRate=(d1$nDiscordant+d1$nConcordant)/d1$totnComparisons
# # # d1$Gl=as.factor(d1$Gl)
# # # d1$GcMethod=as.factor(d1$GcMethod)
# # # d1$Betavar=as.factor(d1$Betavar)
# # # d1$Depth=as.factor(d1$Depth)
# # # # d1$GQ=as.factor(d1$GQ)


# # # d1$ErrorRate[is.na(d1$ErrorRate)]<-0
# # # d1$CallRate[is.na(d1$CallRate)]<-0
# # # d1%>%
# # # ggplot(aes(x = GQ, y = ErrorRate, color = Gl, fill=Gl))+
# # # facet_grid(Depth~Betavar,scale="free")+
# # # geom_col()


# # # ggplot(d1,aes(x = GQ, y = ErrorRate, color = Gl, fill=Gl))+
# # # geom_point()
# # # geom_line()+
# # # facet_grid(Depth~Betavar,scale="free")

# # # d1<-d
# # # d2<-(d1%>%
# # # group_by(Sample, Rep,Depth, Betavar, Gl,GcMethod)%>%
# # # mutate(
# # #       ratioCallsTotal=cumsum(nSitesGQ)/nSitesNon0DP,
# # #       min_GQ_threshold= GQ,
# # #       DiscordanceRateOverCumsums= cumsum(nDiscordant)/(cumsum(nDiscordant)+cumsum(nConcordant)),
# # #       )%>%
# # #       ungroup()%>%
# # #       group_by(Depth, Betavar, Gl,GcMethod,min_GQ_threshold)%>%
# # #       summarize(DiscordanceRateOverCumsums=mean(DiscordanceRateOverCumsums,na.rm=TRUE),
# # #                 ratioCallsTotal=mean(ratioCallsTotal,na.rm=TRUE),
# # #                 min_GQ_threshold=mean(min_GQ_threshold))%>%ungroup()%>%    
# # # mutate(Gl=as.factor(Gl),
# # # GcMethod=as.factor(GcMethod),
# # # Betavar=factor(Betavar,levels=c(0,6,5)),
# # # Depth=as.factor(Depth))%>%
# # # ungroup()
# # # )


# # # d1<-d
# # # d2<-(d1%>%
# # # arrange(desc(GQ)) %>%
# # # group_by(Sample, Rep,Depth, Betavar, Gl,GcMethod)%>%
# # # mutate(
# # #       ratioCallsTotal=cumsum(nSitesGQ)/nSitesNon0DP,
# # #       min_GQ_threshold= GQ,
# # #       DiscordanceRateOverCumsums= cumsum(nDiscordant)/(cumsum(nDiscordant)+cumsum(nConcordant)),
# # #       )%>%
# # #       ungroup()%>%
# # #       group_by(Depth, Betavar, Gl,GcMethod,min_GQ_threshold)%>%
# # #       summarize(DiscordanceRateOverCumsums=mean(DiscordanceRateOverCumsums,na.rm=TRUE),
# # #                 ratioCallsTotal=mean(ratioCallsTotal,na.rm=TRUE),
# # #                 min_GQ_threshold=mean(min_GQ_threshold))%>%ungroup()%>%    
# # # mutate(Gl=as.factor(Gl),
# # # GcMethod=as.factor(GcMethod),
# # # Betavar=factor(Betavar,levels=c(0,6,5)),
# # # Depth=as.factor(Depth))%>%
# # # ungroup()
# # # )

# # # d2%>%
# # # mutate(min_GQ_threshold=as.factor(min_GQ_threshold))%>%
# # # ggplot(aes(x=min_GQ_threshold, y = DiscordanceRateOverCumsums, color = Gl, fill=Gl))+
# # # geom_col()+
# # # facet_wrap(Depth~Betavar,labeller=labeller1fn,scale="free")


# # # rm(plt)
# # # # if DiscordanceRate is NA, remove the row
# # # d3<-d2
# # # (plt<-(
# # # ggplot(data=d3,aes(x = ratioCallsTotal, y = DiscordanceRateOverCumsums, color = Depth, fill=Depth,linetype=Gl))+
# # # # ggplot(data=d3,aes(x = ratioCallsTotal, y = NonConcordanceRateOverCumsums, color = Depth, fill=Depth,linetype=Gl))+
# # # # ggplot(data=d3,aes(x = ratioCallsTotal, y = NonConcordanceRateOverTotalNon0DP,color=Depth,fill=Depth,linetype=Gl))+
# # # # ggplot(data=d3,aes(x = ratioCallsTotal, y = DiscordanceRateOverTotalNon0DP, color = Depth, fill=Depth,linetype=Gl))+
# # #   geom_line(linewidth=1)+
# # #   # add the ratioCallsIncludedSoFar for GQ=20 line
# # #   # geom_point(data=d3[d3$min_GQ_threshold==20,],aes(x = ratioCallsIncludedSoFar,shape=Gl),size=4)+
# # #   # geom_vline(data=d3[d3$min_GQ_threshold==20,],aes(xintercept = ratioCallsIncludedSoFar,color=Depth,linetype=Gl), size=1,alpha=0.7)+
# # #   geom_point(data=d3[d3$min_GQ_threshold==20,],aes(x = ratioCallsTotal),size=3)+
# # #   # geom_vline(data=d3[d3$min_GQ_threshold==20,],aes(xintercept = ratioCallsTotal,color=Depth,linetype=Gl), size=1,alpha=0.7)+
# # #   # facet_grid(.~GcMethod,labeller=labeller1fn)+
# # #   facet_grid(.~Betavar,labeller=labeller1fn)+
# # #   xlim(0,1)+
# # #   # ylim(0,max(d2$DiscordanceRate,na.rm=TRUE))+
# # #   labs(
# # #     x="Call rate",
# # #     y="Error rate",
# # #     color="Depth",
# # #     linetype="Genotype likelihood method",
# # #     )+
# # #   theme_bw()+
# # #   theme(legend.position = "top")+
# # #   theme(legend.box = "vertical")+
# # #   theme(aspect.ratio=1)+
# # #   guides(fill=FALSE)+
# # #   guides(shape=FALSE)+
# # #   theme(strip.background = element_rect(colour="black", fill="white"))+
# # #   theme(text = element_text(size=10))+
# # #   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
# # #   scale_color_brewer(palette="RdYlGn")+
# # #   scale_fill_brewer(palette="RdYlGn")
# # # ))


# # # d2%>%
# # # ggplot(aes(x = ratioCallsTotal, y = DiscordanceRateOverCumsums, color = Depth, fill=Depth,linetype=Gl))+
# # #   geom_line()+
# # #   # geom_point(data=d2[d2$min_GQ_threshold==20,],aes(x = ratioCallsTotal),size=3)+
# # #   facet_grid(.~Betavar,labeller=labeller1fn)
# # #   xlim(0,1)+
# # #   labs(
# # #     x="Call rate",
# # #     y="Error rate",
# # #     color="Depth",
# # #     linetype="Genotype likelihood method",
# # #     )+
# # #   theme_bw()+
# # #   theme(legend.position = "top")+
# # #   theme(legend.box = "vertical")+
# # #   theme(aspect.ratio=1)+
# # #   guides(fill=FALSE)+
# # #   guides(shape=FALSE)+
# # #   theme(strip.background = element_rect(colour="black", fill="white"))+
# # #   theme(text = element_text(size=10))+
# # #   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))

# # # d2 %>% 
# # # ggplot(aes(x = ratioCallsTotal, y = DiscordanceRateOverCumsums, color = Gl, fill=Gl))+
# # # geom_line()+
# # # geom_point(data=d2[d2$min_GQ_threshold==20,],aes(x = ratioCallsTotal),size=3)+
# # # facet_wrap(Depth~Betavar,labeller=labeller1fn,scale="free")+
# # # theme_bw()



# # # # d# 1<-d[d$GcMethod=="genotype_calling_perpop",]
# # # # d2<-(d1%>%
# # # # ungroup()%>%
# # # # #sort by GQ
# # # # arrange(desc(GQ)) %>%
# # # # group_by(Sample,Rep,Depth, Betavar, Gl,GcMethod)%>%
# # # # mutate(
# # # #       DiscordanceRateOverCumsums= cumsum(nDiscordant)/(cumsum(nDiscordant)+cumsum(nConcordant)),
# # # #       ratioCallsTotal=cumsum(nSitesGQ)/nSitesNon0DP,
# # # #       min_GQ_threshold= GQ)%>%
# # # #   ungroup()%>%
# # # # group_by(Depth, Betavar, Gl,GcMethod,min_GQ_threshold)%>%
# # # # summarize(DiscordanceRateOverCumsums=mean(DiscordanceRateOverCumsums,na.rm=TRUE),
# # # #           ratioCallsTotal=mean(ratioCallsTotal,na.rm=TRUE),
# # # #           min_GQ_threshold=mean(min_GQ_threshold))%>%ungroup()%>%
# # # #   mutate(Gl=as.factor(Gl),
# # # #   GcMethod=as.factor(GcMethod),
# # # #   Betavar=as.factor(Betavar),
# # # #   Depth=as.factor(Depth))%>%
# # # #   ungroup()
# # # # )


# # # # d2%>%
# # # # ggplot(aes(x = ratioCallsTotal, y = DiscordanceRateOverCumsums, color= Betavar,linetype=Gl))+
# # # #   # geom_point()+
# # # #   geom_line()+
# # # #   facet_wrap(.~Depth,scale="free",labeller=labeller1fn)+
# # # #   geom_point(data=d2[d2$min_GQ_threshold==20,],aes(x = ratioCallsTotal,shape=Gl),size=3)+
# # # #   xlim(0,1)+
# # # #   labs(
# # # #     x="Call rate",
# # # #     y="Error rate",
# # # #     # color="Depth",
# # # #     # linetype="Genotype likelihood method",
# # # #     # color="Genotype likelihood method",
# # # #     )+
# # # #   theme_bw()+
# # # #   theme(legend.position = "top")+
# # # #   theme(legend.box = "vertical")+
# # # #   theme(aspect.ratio=1)+
# # # #   guides(fill=FALSE)+
# # # #   guides(shape=FALSE)+
# # # #   theme(strip.background = element_rect(colour="black", fill="white"))+
# # # #   theme(text = element_text(size=10))
# # # #   scale_color_manual(values=c("red","blue"),labels=c("Genotype likelihood method 2 with phred-scaled error rate","Genotype likelihood method 2 with precise error rate"))+
# # # #   #make color legend 2 row
# # # #   guides(color=guide_legend(nrow=2))+
# # # #   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))


# # # d1$Sample<-paste0(d1$Sample,"_",d1$Rep)





# # # d1%>%
# # # group_by(Depth, Betavar, Gl,GcMethod,GQ)%>%
# # # summarize(nDiscordant= sum(nDiscordant),
# # # nConcordant= sum(nConcordant),
# # # nSitesGQ= sum(nSitesGQ))%>%
# # # ungroup()%>% 
# # #   mutate(DiscordanceRate= nDiscordant/(nDiscordant+nConcordant),
# # #     GQ=as.factor(GQ),
# # #     Gl=as.factor(Gl),
# # #   GcMethod=as.factor(GcMethod),
# # #   Betavar=as.factor(Betavar),
# # #   Depth=as.factor(Depth))%>%
# # # ungroup()%>% 
# # # ggplot(aes(x=GQ,y=DiscordanceRate, color=Depth ,fill=Depth))+
# # # # geom_col()+
# # # geom_density()+
# # # facet_grid(Gl+Betavar~Depth,scale="free")


# # # d1%>%
# # # group_by(Sample,Depth, Betavar, Gl,GcMethod,GQ)%>%
# # # summarize(nDiscordant= sum(nDiscordant)) %>%ungroup()->d2


# # # d1%>%filter(Sample == 1  & Rep == 1 & Depth == 1 & Betavar == 0 & Gl == "gl1" & GcMethod == "genotype_calling_perpop" & GQ == 10)
# # # d2%>%filter(Depth == 1 & Betavar == 0 & Gl == "gl1" & GcMethod == "genotype_calling_perpop" & GQ == 10)

# # # mutate(GQ=as.factor(GQ))%>%
# # # ggplot(aes(x=GQ,y=nDiscordant))+
# # # geom_col()


# # # d1%>% 
# # # group_by(Sample,Depth, Betavar, Gl,GcMethod,GQ)%>%

# # # d2<-(d1%>%
# # # #sort by GQ
# # # arrange(desc(GQ)) %>%
# # # group_by(Depth, Betavar, Gl,GcMethod)%>%
# # # mutate(
# # #       DiscordanceRateOverCumsums= cumsum(nDiscordant)/(cumsum(nDiscordant)+cumsum(nConcordant)),
# # #       ratioCallsIncludedSoFar=cumsum(nSitesGQ)/sum(nSitesGQ),
# # #       ratioCallsTotal=cumsum(nSitesGQ)/nSitesNon0DP,
# # #       min_GQ_threshold= GQ)%>%
# # #   ungroup()%>%
# # #   mutate(Gl=as.factor(Gl),
# # #   GcMethod=as.factor(GcMethod),
# # #   Betavar=as.factor(Betavar),
# # #   Depth=as.factor(Depth))
# # # )

# # # d2%>% 
# # # ggplot(aes(x = ratioCallsTotal, y = DiscordanceRateOverCumsums, color = Gl, fill=Gl))+
# # # geom_line()+
# # # facet_wrap(Depth~Betavar,labeller=labeller1fn,scale="free")





# # # d2<-(d1%>%
# # # #sort by GQ
# # # arrange(desc(GQ)) %>%
# # # group_by(Depth, Betavar, Gl,GcMethod)%>%
# # # mutate(
# # #       DiscordanceRateOverCumsums= cumsum(nDiscordant)/(cumsum(nDiscordant)+cumsum(nConcordant)),
# # #       ratioCallsIncludedSoFar=cumsum(nSitesGQ)/sum(nSitesGQ),
# # #       ratioCallsTotal=cumsum(nSitesGQ)/nSitesNon0DP,
# # #       min_GQ_threshold= GQ)%>%
# # #   ungroup()%>%
# # #   mutate(Gl=as.factor(Gl),
# # #   GcMethod=as.factor(GcMethod),
# # #   Betavar=as.factor(Betavar),
# # #   Depth=as.factor(Depth))
# # # )

# # # # d2

# # # # d2%>%
# # # # ggplot(aes(x = ratioCallsTotal, y = DiscordanceRateOverCumsums, color = Depth, fill=Depth,linetype=Gl))+
# # # #   geom_line()+
# # # #   # geom_point(data=d2[d2$min_GQ_threshold==20,],aes(x = ratioCallsTotal),size=3)+
# # # #   facet_grid(.~Betavar,labeller=labeller1fn)
# # # #   xlim(0,1)+
# # # #   labs(
# # # #     x="Call rate",
# # # #     y="Error rate",
# # # #     color="Depth",
# # # #     linetype="Genotype likelihood method",
# # # #     )+
# # # #   theme_bw()+
# # # #   theme(legend.position = "top")+
# # # #   theme(legend.box = "vertical")+
# # # #   theme(aspect.ratio=1)+
# # # #   guides(fill=FALSE)+
# # # #   guides(shape=FALSE)+
# # # #   theme(strip.background = element_rect(colour="black", fill="white"))+
# # # #   theme(text = element_text(size=10))+
# # # #   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))



# # # # d2<-(d1%>%
# # # # #sum across samples for each rep
# # # # group_by(Rep,Sample,Depth, Betavar, Gl, GQ,GcMethod)%>%
# # # # summarize(nDiscordant= sum(nDiscordant),
# # # # nConcordant= sum(nConcordant),
# # # # nSitesGQ= sum(nSitesGQ))%>%
# # # # ungroup()%>%
# # # # #sort by GQ
# # # # arrange(desc(GQ)) %>%
# # # # group_by(Rep,Depth, Betavar, Gl,GcMethod)%>%
# # # # mutate(
# # # #       DiscordanceRateOverCumsums= cumsum(nDiscordant)/(cumsum(nDiscordant)+cumsum(nConcordant)),
# # # #       ConcordanceRateOverCumsums= cumsum(nConcordant)/(cumsum(nDiscordant)+cumsum(nConcordant)),
# # # #       DiscordanceRateOverCumsums2= cumsum(nDiscordant)/(1+(cumsum(nDiscordant)+cumsum(nConcordant))),
# # # #       ConcordanceRateOverCumsums2= cumsum(nConcordant)/(1+(cumsum(nDiscordant)+cumsum(nConcordant))),
# # # #       DiscordanceCumSum=cumsum(nDiscordant),
# # # #       DiscordanceRatio=cumsum(nDiscordant)/sum(nDiscordant),
# # # #       nCallsIncludedSoFar=cumsum(nSitesGQ),
# # # #       ratioCallsIncludedSoFar=cumsum(nSitesGQ)/sum(nSitesGQ),
# # # #       min_GQ_threshold= GQ)%>%
# # # # ungroup()%>%
# # # # #get replicate means
# # # # group_by(Depth, Betavar, Gl, GcMethod, min_GQ_threshold)%>%
# # # # summarize(
# # # #   meanDiscordanceRateOverCumsums=mean(DiscordanceRateOverCumsums),
# # # #   meanConcordanceRateOverCumsums=mean(ConcordanceRateOverCumsums),
# # # #   meanDiscordanceRateOverCumsums2=mean(DiscordanceRateOverCumsums2),
# # # #   meanConcordanceRateOverCumsums2=mean(ConcordanceRateOverCumsums2),
# # # #   meanDiscordanceCumSum=mean(DiscordanceCumSum),
# # # #   meanDiscordanceRatio=mean(DiscordanceRatio),
# # # #   meanNCallsIncludedSoFar=mean(nCallsIncludedSoFar),
# # # #   meanRatioCallsIncludedSoFar=mean(ratioCallsIncludedSoFar),
# # # #   meanMin_GQ_threshold=mean(min_GQ_threshold)
# # # #   )%>%
# # # #   ungroup()%>%
# # # #   mutate(Gl=as.factor(Gl),
# # # #   GcMethod=as.factor(GcMethod),
# # # #   Betavar=as.factor(Betavar),
# # # #   Depth=as.factor(Depth))
# # # # )





























# # # # # method 0) use nsitesnon0dp
# # # # rm(d1,d2,d3,plt);gc()
# # # # d1<-d
# # # # d1$Sample<-NULL
# # # # d2<-(d1%>%
# # # # #sum across samples AND replicates
# # # # group_by(Depth, Betavar, Gl, GcMethod, GQ)%>%
# # # # summarize(nDiscordant= sum(nDiscordant),
# # # # nConcordant= sum(nConcordant),
# # # # nSitesNon0DP= sum(nSitesNon0DP),
# # # # nSitesGQ= sum(nSitesGQ))%>%
# # # # ungroup()%>%
# # # # #sort by GQ
# # # # arrange(desc(GQ)) %>%
# # # # group_by(Depth, Betavar, Gl,GcMethod)%>%
# # # # mutate(
# # # #       DiscordanceRateOverTotalNon0DP=cumsum(nDiscordant)/nSitesNon0DP,
# # # #       DiscordanceRateOverCumsums= cumsum(nDiscordant)/(cumsum(nDiscordant)+cumsum(nConcordant)),
# # # #       ConcordanceRateOverCumsums= cumsum(nConcordant)/(cumsum(nDiscordant)+cumsum(nConcordant)),
# # # #       DiscordanceRateOverCumsums2= cumsum(nDiscordant)/(1+(cumsum(nDiscordant)+cumsum(nConcordant))),
# # # #       ConcordanceRateOverCumsums2= cumsum(nConcordant)/(1+(cumsum(nDiscordant)+cumsum(nConcordant))),
# # # #       DiscordanceCumSum=cumsum(nDiscordant),
# # # #       DiscordanceRatio=cumsum(nDiscordant)/sum(nDiscordant),
# # # #       nCallsIncludedSoFar=cumsum(nSitesGQ),
# # # #       ratioCallsIncludedSoFar=cumsum(nSitesGQ)/sum(nSitesGQ),
# # # #       ratioCallsTotal=cumsum(nSitesGQ)/nSitesNon0DP,
# # # #       min_GQ_threshold= GQ)%>%
# # # # ungroup()%>%
# # # # mutate(Gl=as.factor(Gl),
# # # # GcMethod=as.factor(GcMethod),
# # # # Betavar=factor(Betavar,levels=c(0,6,5)),
# # # # Depth=as.factor(Depth)))




# # # # rm(plt)
# # # # # if DiscordanceRate is NA, remove the row
# # # # d3<-d2[d2$GcMethod=="genotype_calling_perpop",]
# # # # (plt<-(
# # # # ggplot(data=d3,aes(x = ratioCallsTotal, y = DiscordanceRateOverCumsums, color = Gl, fill=Gl,linetype=Gl))+
# # # # # ggplot(data=d3,aes(x = ratioCallsTotal, y = DiscordanceRateOverTotalNon0DP, color = Depth, fill=Depth,linetype=Gl))+
# # # #   geom_line(linewidth=1)+
# # # #   # add the ratioCallsIncludedSoFar for GQ=20 line
# # # #   # geom_point(data=d3[d3$min_GQ_threshold==20,],aes(x = ratioCallsIncludedSoFar,shape=Gl),size=4)+
# # # #   # geom_vline(data=d3[d3$min_GQ_threshold==20,],aes(xintercept = ratioCallsIncludedSoFar,color=Depth,linetype=Gl), size=1,alpha=0.7)+
# # # #   geom_point(data=d3[d3$min_GQ_threshold==20,],aes(x = ratioCallsTotal),size=3)+
# # # #   # geom_vline(data=d3[d3$min_GQ_threshold==20,],aes(xintercept = ratioCallsTotal,color=Depth,linetype=Gl), size=1,alpha=0.7)+
# # # #   facet_wrap(Betavar~Depth,labeller=labeller1fn)+
# # # #   xlim(0,1)+
# # # #   # ylim(0,max(d2$DiscordanceRate,na.rm=TRUE))+
# # # #   labs(
# # # #     x="Call rate",
# # # #     y="Error rate",
# # # #     color="Depth",
# # # #     linetype="Genotype likelihood method",
# # # #     )+
# # # #   theme_bw()+
# # # #   theme(legend.position = "top")+
# # # #   theme(legend.box = "vertical")+
# # # #   theme(aspect.ratio=1)+
# # # #   guides(fill=FALSE)+
# # # #   guides(shape=FALSE)+
# # # #   theme(strip.background = element_rect(colour="black", fill="white"))+
# # # #   theme(text = element_text(size=10))+
# # # #   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))
# # # # ))

# # # # # save_plt(plt=plt,plot_filename="plot_x-ratioCallsIncludedSoFar_y-DiscordanceRateOverCumsums_grid-BetavarDepth_gcAllSamples.png",plot_outdir=plot_outdir,
# # # # save_plt(plt=plt,plot_filename="plot_x-ratioCallsTotal_y-DiscordanceRateOverCumsums_color-Depth_fill-Depth_linetype-Gl_grid-GcMethod_point-GQ20.png",plot_outdir=plot_outdir,
# # # #          overwrite=TRUE,width=10,height=10,units="cm",dpi=300,scale=1.2)




# # # # # # save as RData
# # # # # save(d,file=paste0(plot_outdir,"data_",today,".RData"))


# # # # # simdir="../sim/"
# # # # # sd<-read.csv(paste0(simdir,"sim_vcfgl_2312/model_OutOfAfrica_3G09/stats/nSites_non0dp/sim_vcfgl_2312-OutOfAfrica_3G09_nSitesNon0DP.csv"),sep=",",header=TRUE)

# # # # # # convert unique names to indices in unique names list 
# # # # # sd$Sample<-as.numeric(factor(sd$Sample,levels=unique(sd$Sample)))
# # # # # # colnames(sd)
# # # # # # merge this dataframe with data d by columns Sample Rep Depth
# # # # # d<-merge(d,sd,by=c("Sample","Rep","Depth"))
# # # # # # save locally
# # # # # save(d,file="~/Projects/VCFGL/data_231222.RData")


# # # # # load RData
# # # # # load(paste0(plot_outdir,"data_",today,".RData"))

# # # # # load locally
# # # # # load("~/Projects/VCFGL/data_231222.RData")

# # # # # Columns:
# # # # # Sample = Sample index
# # # # # GQ = Genotype quality scores of the called genotypes
# # # # # nDiscordant = the number of discordant genotypes 
# # # # # nConcordant = the number of concordant genotypes
# # # # # nSitesGQ (= nDiscordant + nConcordant) = total number of sites with data for that sample that has a genotype call with the given GQ
# # # # # totnComparisons = the total number of genotypes compared for that sample in total (= nSitesGQ * nGQs for that sample)
# # # # # Rep = replicate index
# # # # # Depth = sequencing depth
# # # # # Betavar = beta variance
# # # # # GcMethod = genotype calling method
# # # # # Gl = GL method ID (1 or 2)



# # # # ################################################################################



# # # # # Aim:
# # # # # 1) Show the differences between the two GL methods

# # # # # colnames(d1)
# # # # #  [1] "Sample"          "GQ"              "nDiscordant"     "nConcordant"    
# # # # #  [5] "totnComparisons" "Rep"             "Depth"           "Betavar"        
# # # # #  [9] "GcMethod"        "Gl"              "nSitesGQ"     
# # # # # group_by(Sample,GQ,Rep,Depth,Betavar,GcMethod,Gl)






# # # # ################################################################################


# # # # ### METHOD 0


# # # # # method 0) use nsitesnon0dp
# # # # rm(d1,d2,d3,plt);gc()
# # # # d1<-d1[d1$Betavar==0,]
# # # # d1$Sample<-NULL
# # # # d2<-(d1%>%
# # # # #sum across samples AND replicates
# # # # group_by(Depth, Betavar, Gl, GcMethod, GQ)%>%
# # # # summarize(nDiscordant= sum(nDiscordant),
# # # # nConcordant= sum(nConcordant),
# # # # nSitesNon0DP= sum(nSitesNon0DP),
# # # # nSitesGQ= sum(nSitesGQ))%>%
# # # # ungroup()%>%
# # # # #sort by GQ
# # # # arrange(desc(GQ)) %>%
# # # # group_by(Depth, Betavar, Gl,GcMethod)%>%
# # # # mutate(
# # # #       DiscordanceRateOverTotalNon0DP=cumsum(nDiscordant)/nSitesNon0DP,
# # # #       DiscordanceRateOverCumsums= cumsum(nDiscordant)/(cumsum(nDiscordant)+cumsum(nConcordant)),
# # # #       ConcordanceRateOverCumsums= cumsum(nConcordant)/(cumsum(nDiscordant)+cumsum(nConcordant)),
# # # #       DiscordanceRateOverCumsums2= cumsum(nDiscordant)/(1+(cumsum(nDiscordant)+cumsum(nConcordant))),
# # # #       ConcordanceRateOverCumsums2= cumsum(nConcordant)/(1+(cumsum(nDiscordant)+cumsum(nConcordant))),
# # # #       DiscordanceCumSum=cumsum(nDiscordant),
# # # #       DiscordanceRatio=cumsum(nDiscordant)/sum(nDiscordant),
# # # #       nCallsIncludedSoFar=cumsum(nSitesGQ),
# # # #       ratioCallsIncludedSoFar=cumsum(nSitesGQ)/sum(nSitesGQ),
# # # #       ratioCallsTotal=cumsum(nSitesGQ)/nSitesNon0DP,
# # # #       min_GQ_threshold= GQ)%>%
# # # # ungroup()%>%
# # # # mutate(Gl=as.factor(Gl),
# # # # GcMethod=as.factor(GcMethod),
# # # # Betavar=factor(Betavar,levels=c(0,6,5)),
# # # # Depth=as.factor(Depth)))



# # # # rm(plt)
# # # # # if DiscordanceRate is NA, remove the row
# # # # d3<-d2
# # # # (plt<-(
# # # # ggplot(data=d3,aes(x = ratioCallsTotal, y = DiscordanceRateOverCumsums, color = Depth, fill=Depth,linetype=Gl))+
# # # # # ggplot(data=d3,aes(x = ratioCallsTotal, y = DiscordanceRateOverTotalNon0DP, color = Depth, fill=Depth,linetype=Gl))+
# # # #   geom_line(linewidth=1)+
# # # #   # add the ratioCallsIncludedSoFar for GQ=20 line
# # # #   # geom_point(data=d3[d3$min_GQ_threshold==20,],aes(x = ratioCallsIncludedSoFar,shape=Gl),size=4)+
# # # #   # geom_vline(data=d3[d3$min_GQ_threshold==20,],aes(xintercept = ratioCallsIncludedSoFar,color=Depth,linetype=Gl), size=1,alpha=0.7)+
# # # #   geom_point(data=d3[d3$min_GQ_threshold==20,],aes(x = ratioCallsTotal),size=3)+
# # # #   # geom_vline(data=d3[d3$min_GQ_threshold==20,],aes(xintercept = ratioCallsTotal,color=Depth,linetype=Gl), size=1,alpha=0.7)+
# # # #   facet_grid(.~GcMethod,labeller=labeller1fn)+
# # # #   xlim(0,1)+
# # # #   # ylim(0,max(d2$DiscordanceRate,na.rm=TRUE))+
# # # #   labs(
# # # #     x="Call rate",
# # # #     y="Error rate",
# # # #     color="Depth",
# # # #     linetype="Genotype likelihood method",
# # # #     )+
# # # #   theme_bw()+
# # # #   theme(legend.position = "top")+
# # # #   theme(legend.box = "vertical")+
# # # #   theme(aspect.ratio=1)+
# # # #   guides(fill=FALSE)+
# # # #   guides(shape=FALSE)+
# # # #   theme(strip.background = element_rect(colour="black", fill="white"))+
# # # #   theme(text = element_text(size=10))+
# # # #   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
# # # #   scale_color_brewer(palette="RdYlGn")+
# # # #   scale_fill_brewer(palette="RdYlGn")
# # # # ))

# # # # # save_plt(plt=plt,plot_filename="plot_x-ratioCallsIncludedSoFar_y-DiscordanceRateOverCumsums_grid-BetavarDepth_gcAllSamples.png",plot_outdir=plot_outdir,
# # # # save_plt(plt=plt,plot_filename="plot_x-ratioCallsTotal_y-DiscordanceRateOverCumsums_color-Depth_fill-Depth_linetype-Gl_grid-GcMethod_point-GQ20.png",plot_outdir=plot_outdir,
# # # #          overwrite=TRUE,width=10,height=10,units="cm",dpi=300,scale=1.2)


# # # # rm(plt)
# # # # # if DiscordanceRate is NA, remove the row
# # # # d3<-d2
# # # # (plt<-(
# # # # ggplot(data=d3,aes(x = ratioCallsTotal, y = DiscordanceRateOverTotalNon0DP, color = Depth, fill=Depth,linetype=Gl))+
# # # #   geom_line(linewidth=1)+
# # # #   # add the ratioCallsIncludedSoFar for GQ=20 line
# # # #   # geom_point(data=d3[d3$min_GQ_threshold==20,],aes(x = ratioCallsIncludedSoFar,shape=Gl),size=4)+
# # # #   # geom_vline(data=d3[d3$min_GQ_threshold==20,],aes(xintercept = ratioCallsIncludedSoFar,color=Depth,linetype=Gl), size=1,alpha=0.7)+
# # # #   geom_point(data=d3[d3$min_GQ_threshold==20,],aes(x = ratioCallsTotal),size=3)+
# # # #   # geom_vline(data=d3[d3$min_GQ_threshold==20,],aes(xintercept = ratioCallsTotal,color=Depth,linetype=Gl), size=1,alpha=0.7)+
# # # #   facet_grid(.~GcMethod,labeller=labeller1fn)+
# # # #   xlim(0,1)+
# # # #   labs(
# # # #     x="Call rate",
# # # #     y="Error rate",
# # # #     color="Depth",
# # # #     linetype="Genotype likelihood method",
# # # #     )+
# # # #   theme_bw()+
# # # #   theme(legend.position = "top")+
# # # #   theme(legend.box = "vertical")+
# # # #   theme(aspect.ratio=1)+
# # # #   guides(fill=FALSE)+
# # # #   guides(shape=FALSE)+
# # # #   theme(strip.background = element_rect(colour="black", fill="white"))+
# # # #   theme(text = element_text(size=10))+
# # # #   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
# # # #   scale_color_brewer(palette="RdYlGn")+
# # # #   scale_fill_brewer(palette="RdYlGn")
# # # # ))

# # # # save_plt(plt=plt,plot_filename="plot_x-ratioCallsTotal_y-DiscordanceRateOverTotalNon0DP_color-Depth_fill-Depth_linetype-Gl_grid-GcMethod_point-GQ20.png",plot_outdir=plot_outdir,
# # # #          overwrite=TRUE,width=10,height=10,units="cm",dpi=300,scale=1.2)


# # # # rm(plt)
# # # # # if DiscordanceRate is NA, remove the row
# # # # d3<-d2
# # # # (plt<-(
# # # # ggplot(data=d3,aes(x = ratioCallsIncludedSoFar, y = DiscordanceRateOverCumsums, color = Depth, fill=Depth,linetype=Gl))+
# # # #   geom_line(linewidth=1)+
# # # #   # add the ratioCallsIncludedSoFar for GQ=20 line
# # # #   geom_point(data=d3[d3$min_GQ_threshold==20,],aes(x = ratioCallsIncludedSoFar,shape=Gl),size=3)+
# # # #   # geom_vline(data=d3[d3$min_GQ_threshold==20,],aes(xintercept = ratioCallsIncludedSoFar,color=Depth,linetype=Gl), size=1,alpha=0.7)+
# # # #   # geom_point(data=d3[d3$min_GQ_threshold==20,],aes(x = ratioCallsTotal),size=3)+
# # # #   # geom_vline(data=d3[d3$min_GQ_threshold==20,],aes(xintercept = ratioCallsTotal,color=Depth,linetype=Gl), size=1,alpha=0.7)+
# # # #   facet_grid(.~GcMethod,labeller=labeller1fn)+
# # # #   xlim(0,1)+
# # # #   labs(
# # # #     x="Call rate",
# # # #     y="Error rate",
# # # #     color="Depth",
# # # #     linetype="Genotype likelihood method",
# # # #     )+
# # # #   theme_bw()+
# # # #   theme(legend.position = "top")+
# # # #   theme(legend.box = "vertical")+
# # # #   theme(aspect.ratio=1)+
# # # #   guides(fill=FALSE)+
# # # #   guides(shape=FALSE)+
# # # #   theme(strip.background = element_rect(colour="black", fill="white"))+
# # # #   theme(text = element_text(size=10))+
# # # #   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
# # # #   scale_color_brewer(palette="RdYlGn")+
# # # #   scale_fill_brewer(palette="RdYlGn")
# # # # ))

# # # # save_plt(plt=plt,plot_filename="plot_x-ratioCallsIncludedSoFar_y-DiscordanceRateOverCumsums_color-Depth_fill-Depth_linetype-Gl_grid-GcMethod_point-GQ20.png",plot_outdir=plot_outdir,
# # # #          overwrite=TRUE,width=10,height=10,units="cm",dpi=300,scale=1.2)


# # # # rm(plt)
# # # # # if DiscordanceRate is NA, remove the row
# # # # d3<-d2
# # # # (plt<-(
# # # # ggplot(data=d3,aes(x = ratioCallsIncludedSoFar, y = DiscordanceRateOverCumsums2, color = Depth, fill=Depth,linetype=Gl))+
# # # #   geom_line(linewidth=1)+
# # # #   # add the ratioCallsIncludedSoFar for GQ=20 line
# # # #   geom_point(data=d3[d3$min_GQ_threshold==20,],aes(x = ratioCallsIncludedSoFar,shape=Gl),size=3)+
# # # #   # geom_vline(data=d3[d3$min_GQ_threshold==20,],aes(xintercept = ratioCallsIncludedSoFar,color=Depth,linetype=Gl), size=1,alpha=0.7)+
# # # #   # geom_point(data=d3[d3$min_GQ_threshold==20,],aes(x = ratioCallsTotal),size=3)+
# # # #   # geom_vline(data=d3[d3$min_GQ_threshold==20,],aes(xintercept = ratioCallsTotal,color=Depth,linetype=Gl), size=1,alpha=0.7)+
# # # #   facet_grid(.~GcMethod,labeller=labeller1fn)+
# # # #   xlim(0,1)+
# # # #   labs(
# # # #     x="Call rate",
# # # #     y="Error rate",
# # # #     color="Depth",
# # # #     linetype="Genotype likelihood method",
# # # #     )+
# # # #   theme_bw()+
# # # #   theme(legend.position = "top")+
# # # #   theme(legend.box = "vertical")+
# # # #   theme(aspect.ratio=1)+
# # # #   guides(fill=FALSE)+
# # # #   guides(shape=FALSE)+
# # # #   theme(strip.background = element_rect(colour="black", fill="white"))+
# # # #   theme(text = element_text(size=10))+
# # # #   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
# # # #   scale_color_brewer(palette="RdYlGn")+
# # # #   scale_fill_brewer(palette="RdYlGn")
# # # # ))



# # # d2[d2$Gl == "precise1_gl2",]$nSitesGQ-d2[d2$Gl == "gl2",]$nSitesGQ


# # # d1<-d1[d1$Gl %in% c("gl2","precise1_gl2"),]
# # # d2<-d2[d2$Betavar=="0",]
# # # rm(plt)
# # # (plt<-ggplot(d2,aes(x=ratioCallsIncludedSoFar,y=DiscordanceRateOverCumsums, color=Gl))+
# # # facet_wrap(Depth~Betavar,scale="free",labeller=labeller1fn)+
# # # # geom_point(data=d2[d2$min_GQ_threshold==20,],aes(x = ratioCallsIncludedSoFar),size=3)+
# # # geom_line()+
# # # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
# # # theme_bw()+
# # # theme(legend.position = "top")+
# # #   theme(text = element_text(size=10))+
# # #   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
# # #   # scale_color_brewer(palette="RdYlGn")+
# # #   scale_fill_brewer(palette="RdYlGn")+
# # # labs(
# # #     x="Call rate",
# # #     y="Error rate",
# # #     color="Depth",
# # #     linetype=""
# # #     )+
# # # # make it square
# # # theme(aspect.ratio=1)+
# # # guides(fill=FALSE)+
# # # guides(shape=FALSE)+
# # # theme(strip.background = element_rect(colour="black", fill="white"))+
# # # #rename linetypes 
# # # # scale_linetype_manual(values=c("solid","dashed","dotted"),labels=c("GL method 1","GL method 2","GL method 2 with precise error rate"))+
# # # # scale_linetype_manual(values=c("gl1"="solid","gl2"="dashed","precise1_gl2"="dotted"),labels=c("gl1"="GL method 1","gl2"="GL method 2","precise1_gl2"="GL method 2 with precise error"))+
# # # # scale_color_manual(values=c("0"="black","5"="blue","6"="red"),labels=c("0"="0","5"=expression(10^-5),"6"=expression(10^-6)))+
# # # guides(linetype=guide_legend(nrow=3)))









# # # # method 0) use nsitesnon0dp
# # # rm(d1,d2,d3,plt);gc()
# # # d1<-d
# # # d1$Sample<-NULL
# # # d2<-(d1%>%
# # # #sum across samples AND replicates
# # # group_by(Depth, Betavar, Gl, GcMethod, GQ)%>%
# # # summarize(nDiscordant= sum(nDiscordant),
# # # nConcordant= sum(nConcordant),
# # # nSitesNon0DP= sum(nSitesNon0DP),
# # # nSitesGQ= sum(nSitesGQ))%>%
# # # ungroup()%>%
# # # #sort by GQ
# # # group_by(Depth, Betavar, Gl,GcMethod)%>%
# # # arrange(desc(GQ)) %>%
# # # mutate(
# # #       DiscordanceRateOverTotalNon0DP=cumsum(nDiscordant)/nSitesNon0DP,
# # #       DiscordanceRateOverCumsums= cumsum(nDiscordant)/(cumsum(nDiscordant)+cumsum(nConcordant)),

# # #       NonConcordanceRateOverCumsums= 1- (cumsum(nConcordant)/(cumsum(nDiscordant)+cumsum(nConcordant))),
# # #       NonConcordanceRateOverTotalNon0DP=1-(cumsum(nConcordant)/nSitesNon0DP),

# # #       ConcordanceRateOverCumsums= cumsum(nConcordant)/(cumsum(nDiscordant)+cumsum(nConcordant)),
# # #       DiscordanceRateOverCumsums2= cumsum(nDiscordant)/(1+(cumsum(nDiscordant)+cumsum(nConcordant))),
# # #       ConcordanceRateOverCumsums2= cumsum(nConcordant)/(1+(cumsum(nDiscordant)+cumsum(nConcordant))),
# # #       DiscordanceCumSum=cumsum(nDiscordant),
# # #       DiscordanceRatio=cumsum(nDiscordant)/sum(nDiscordant),
# # #       nCallsIncludedSoFar=cumsum(nSitesGQ),
# # #       ratioCallsIncludedSoFar=cumsum(nSitesGQ)/sum(nSitesGQ),
# # #       ratioCallsTotal=cumsum(nSitesGQ)/nSitesNon0DP,
# # #       min_GQ_threshold= GQ)%>%
# # # ungroup()%>%
# # # mutate(Gl=as.factor(Gl),
# # # GcMethod=as.factor(GcMethod),
# # # Betavar=factor(Betavar,levels=c(0,6,5)),
# # # Depth=as.factor(Depth)))


d3<-d1%>%
group_by(Depth, Betavar, Gl,GcMethod,GQ)%>%
summarize(nDiscordant= sum(nDiscordant),
nConcordant= sum(nConcordant),
nHomToHomDiscordant= sum(nHomToHomDiscordant),
nHomToHetDiscordant= sum(nHomToHetDiscordant),
nHetToHomDiscordant= sum(nHetToHomDiscordant),
nHetToHetDiscordant= sum(nHetToHetDiscordant),
nHomToHomConcordant= sum(nHomToHomConcordant),
nHetToHetConcordant= sum(nHetToHetConcordant),
nSitesGQ=sum(nDiscordant)+sum(nConcordant),
)%>%
ungroup()%>%
mutate(DiscordanceRate= nDiscordant/(nDiscordant+nConcordant),
trueHomDiscordanceRate= (nHomToHomDiscordant+nHomToHetDiscordant)/(nHomToHomDiscordant+nHomToHetDiscordant+nHomToHomConcordant),
trueHetDiscordanceRate= (nHetToHomDiscordant+nHetToHetDiscordant)/(nHetToHomDiscordant+nHetToHetDiscordant+nHetToHetConcordant),
HomToHomDiscordanceRate= nHomToHomDiscordant/(nHomToHomDiscordant+nHomToHomConcordant),
HetToHetDiscordanceRate= nHetToHetDiscordant/(nHetToHetDiscordant+nHetToHetConcordant),
trueHomDiscordanceRateOverAll = (nHomToHomDiscordant+nHomToHetDiscordant)/(nDiscordant+nConcordant),
trueHetDiscordanceRateOverAll = (nHetToHomDiscordant+nHetToHetDiscordant)/(nDiscordant+nConcordant),
HomToHomDiscordanceRateOverAll = nHomToHomDiscordant/(nDiscordant+nConcordant),
HetToHetDiscordanceRateOverAll = nHetToHetDiscordant/(nDiscordant+nConcordant),
HetToHomDiscordanceRateOverAll = nHetToHomDiscordant/(nDiscordant+nConcordant),
HomToHetDiscordanceRateOverAll = nHomToHetDiscordant/(nDiscordant+nConcordant))%>%
mutate(
  across(c('Gl','GcMethod','Depth','Betavar'),as.factor))



#remove gq bins with nSitesGQ==0
d3<-d3[d3$nSitesGQ>0,]


d3%>%ungroup()%>% 
arrange(GQ) %>%
group_by(Depth, Betavar, Gl,GcMethod)%>%
mutate(DiscordanceRate= nDiscordant/(nDiscordant+nConcordant),
      ratioCallsTotal=cumsum(as.numeric(nSitesGQ))/sum(nSitesGQ),
      # min_GQ_threshold= GQ,
      ratioCallsIncludedSoFar=(cumsum(as.numeric(nDiscordant))+cumsum(as.numeric(nConcordant)))/(sum(nDiscordant)+sum(nConcordant)),
      # DiscordanceRateOverCumsums= cumsum(as.numeric(nDiscordant))/(cumsum(as.numeric(nDiscordant))+cumsum(as.numeric(nConcordant))),
      DiscordanceRateOverCumsums= cumsum(as.numeric(nDiscordant))/(sum(nSitesGQ))
      )%>% 
      ungroup()%>%
mutate(across(c('Gl','GcMethod','Depth','Betavar'),as.factor))%>%
mutate(Betavar=factor(Betavar,levels=c("0","7","6","5")))%>%
ungroup()%>% 
ggplot(aes(x = ratioCallsIncludedSoFar, y = DiscordanceRateOverCumsums, color = Betavar, fill=Betavar))+
  facet_wrap(.~Gl+GcMethod+Depth,labeller=labeller1fn, scales="free")+
  geom_step()+
  # geom_line(linewidth=1)+
  # geom_point(data=d3[d3$min_GQ_threshold==20,],aes(x = ratioCallsTotal),size=3)+
  xlim(0,1)+
  labs(
    x="Call rate",
    y="Error rate",
    # color="Depth",
    # linetype="Genotype likelihood method",
    )+
  theme_bw()+
  theme(legend.position = "top")+
  theme(legend.box = "vertical")+
  theme(aspect.ratio=1)+
  guides(fill=FALSE)+
  guides(shape=FALSE)+
  theme(strip.background = element_rect(colour="black", fill="white"))+
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
  # scale_color_brewer(palette="RdYlGn")+
  # scale_fill_brewer(palette="RdYlGn")+
  scale_linetype_manual(values=c("solid","dashed","dotted"),labels=c("GL method 1","GL method 2","GL method 2 with precise error rate"))


d3%>%
arrange(desc(GQ))%>%
group_by(Depth, Betavar, Gl,GcMethod)%>%
mutate(
  ErrorRate=cumsum(as.numeric(nDiscordant))/sum(as.numeric(nSitesGQ)),
  # ErrorRate=cumsum(as.numeric(nDiscordant))/nSitesNon0DP,
  CallRate=cumsum(as.numeric(nSitesGQ))/sum(as.numeric(nSitesGQ))
  # CallRate=cumsum(as.numeric(nConcordant))/nSitesNon0DP
  )%>%
  ungroup()%>%
ggplot(aes(x=CallRate,y=ErrorRate,color=Betavar,fill=Betavar))+
facet_grid(Gl~Depth,scale="free",labeller=labeller1fn)+
# geom_point(data=d2[d2$min_GQ_threshold==20,],aes(x = ratioCallsIncludedSoFar),size=3)+
# geom_step()+
geom_line()+
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
theme_bw()
theme(legend.position = "top")+
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
  # scale_color_brewer(palette="RdYlGn")+
  # scale_fill_brewer(palette="RdYlGn")+
labs(
    # x="Call rate",
    # y="Error rate",
    # color="Depth",
    linetype=""
    )+
theme(aspect.ratio=1)+
guides(color=FALSE)+
guides(shape=FALSE)+
theme(strip.background = element_rect(colour="black", fill="white"))+
# linetype_gls+
guides(linetype=guide_legend(nrow=3))
# color_brewer_betavars+
# fill_brewer_betavars



d3%>%
arrange(desc(GQ))%>%
group_by(Depth, Betavar, Gl,GcMethod)%>%
mutate(
  # ErrorRate=cumsum(as.numeric(nDiscordant))/sum(as.numeric(nDiscordant)),
  ErrorRate=cumsum(as.numeric(nDiscordant))/sum(nDiscordant),
  CallRate=(cumsum(as.numeric(nConcordant))+cumsum(as.numeric(nDiscordant))/nSitesNon0DP)
  # CallRate=cumsum(as.numeric(nConcordant))/nSitesNon0DP
  )%>%
  ungroup()%>%
ggplot(aes(x=CallRate,y=ErrorRate,color=Betavar,fill=Betavar))+
facet_grid(Gl~Depth,scale="free",labeller=labeller1fn)+
# geom_point(data=d2[d2$min_GQ_threshold==20,],aes(x = ratioCallsIncludedSoFar),size=3)+
# geom_step()+
geom_line()+
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
theme_bw()
theme(legend.position = "top")+
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
  # scale_color_brewer(palette="RdYlGn")+
  # scale_fill_brewer(palette="RdYlGn")+
labs(
    # x="Call rate",
    # y="Error rate",
    # color="Depth",
    linetype=""
    )+
theme(aspect.ratio=1)+
guides(color=FALSE)+
guides(shape=FALSE)+
theme(strip.background = element_rect(colour="black", fill="white"))+
# linetype_gls+
guides(linetype=guide_legend(nrow=3))
# color_brewer_betavars+
# fill_brewer_betavars


(d1%>% 
group_by(Depth, Betavar, Gl,GcMethod,GQ)%>%
summarize(nDiscordant= sum(nDiscordant),
nConcordant= sum(nConcordant),
nSitesGQ= sum(nSitesGQ),
nSitesNon0DP= sum(nSitesNon0DP))%>%
ungroup()%>% 
arrange(desc(GQ)) %>%
group_by(Depth, Betavar, Gl,GcMethod)%>%
mutate(DiscordanceRate= nDiscordant/(nDiscordant+nConcordant),
      ratioCallsTotal=cumsum(as.numeric(nSitesGQ))/nSitesNon0DP,
      GQ_threshold= GQ,
      ratioCallsIncludedSoFar=(cumsum(as.numeric(nDiscordant))+cumsum(as.numeric(nConcordant)))/(sum(nDiscordant)+sum(nConcordant)),
      DiscordanceRateOverCumsums= cumsum(as.numeric(nDiscordant))/(cumsum(as.numeric(nDiscordant))+cumsum(as.numeric(nConcordant))),
      DiscordanceRateOverCumsums2= cumsum(as.numeric(nDiscordant))/(1+(cumsum(as.numeric(nDiscordant))+cumsum(as.numeric(nConcordant)))),
      )%>% 
      ungroup()%>%
mutate(across(c('Gl','GcMethod','Depth','Betavar'),as.factor))%>%
mutate(Betavar=factor(Betavar,levels=c("0","7","6","5")))%>%
ungroup())->d2


rm(plt)
(plt<-ggplot(d2,aes(x=ratioCallsTotal,y=DiscordanceRateOverCumsums, color=Depth ,fill=Depth,linetype=Gl))+
facet_wrap(.~Betavar,scale="fixed",labeller=labeller1fn)+
# geom_point(data=d2[d2$min_GQ_threshold==20,],aes(x = ratioCallsIncludedSoFar),size=3)+
geom_step()+
# geom_line()+
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
theme_bw()+
theme(legend.position = "top")+
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
  scale_color_brewer(palette="RdYlGn")+
  scale_fill_brewer(palette="RdYlGn")+
labs(
    x="Call rate",
    y="Error rate",
    color="Depth",
    linetype=""
    )+
theme(aspect.ratio=1)+
guides(fill=FALSE)+
guides(shape=FALSE)+
theme(strip.background = element_rect(colour="black", fill="white"))+
linetype_gls+
guides(linetype=guide_legend(nrow=3)))

plt


save_plt(plt,plot_filename="plot-step_x-ratioCallsTotal_y-DiscordanceRateOverCumsums_grid-Betavar_color-Depth_linetype-Gl_data-GenotypeCallingPerPop.png",plot_outdir=plot_outdir,overwrite=TRUE,width=10,height=10,units="in",scale=0.6,dpi=300)


rm(plt)
(plt<-ggplot(d2,aes(x=ratioCallsTotal,y=DiscordanceRateOverCumsums, color=Depth ,fill=Depth,linetype=Gl))+
facet_wrap(.~Betavar,scale="fixed",labeller=labeller1fn)+
geom_point(data=d2[d2$min_GQ_threshold==20,],aes(x = ratioCallsIncludedSoFar),size=3)+
# geom_step()+
geom_line()+
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
theme_bw()+
theme(legend.position = "top")+
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5))+
  scale_color_brewer(palette="RdYlGn")+
  scale_fill_brewer(palette="RdYlGn")+
labs(
    x="Call rate",
    y="Error rate",
    color="Depth",
    linetype=""
    )+
theme(aspect.ratio=1)+
guides(fill=FALSE)+
guides(shape=FALSE)+
theme(strip.background = element_rect(colour="black", fill="white"))+
linetype_gls+
guides(linetype=guide_legend(nrow=3)))

plt


save_plt(plt,plot_filename="plot-line_x-ratioCallsTotal_y-DiscordanceRateOverCumsums_grid-Betavar_color-Depth_linetype-Gl_data-GenotypeCallingPerPop.png",plot_outdir=plot_outdir,overwrite=TRUE,width=10,height=10,units="in",scale=0.6,dpi=300)

