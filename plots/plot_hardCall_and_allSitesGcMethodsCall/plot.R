###############################################################################
# isinaltinkaya
#

rm(list=ls())
gc()

################################################################################
# LOAD LIBRARIES

# require(data.table) 
# require(readr)

library(gridExtra)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
# library(tidyr) # pivot_longer

library(ggnewscale)

today<-format(Sys.Date(), "%y%m%d")


# setwd(paste0("plots_",today))
plot_outdir<-paste0(getwd(),"/")

color_betavars<-scale_color_manual(labels=c("0"="0","7"=expression(10^-7), "6"=expression(10^-6),"5"=expression(10^-5)), values=c("0"="black","7"="green","6"="blue","5"="red"))
fill_betavars<-scale_fill_manual(labels=c("0"="0","7"=expression(10^-7), "6"=expression(10^-6),"5"=expression(10^-5)), values=c("0"="black","7"="green","6"="blue","5"="red"))
color_brewer_betavars<-scale_color_brewer(palette="Set1",labels=c("0"="0","5"=expression(10^-5),"6"=expression(10^-6),"7"=expression(10^-7)))
fill_brewer_betavars<-scale_fill_brewer(palette="Set1",labels=c("0"="0","5"=expression(10^-5),"6"=expression(10^-6),"7"=expression(10^-7)))
linetype_gls<-scale_linetype_manual(values=c("gl1"="solid","gl2"="dashed","precise1_gl2"="dotted"),labels=c("gl1"="GL method 1","gl2"="GL method 2","precise1_gl2"="GL method 2 with precise error"))

mypalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"), space = "Lab")

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

gcmethod_lut<-c("Genotype calling\nacross populations","Genotype calling\nwithin populations","Basic\ngenotype calling")
names(gcmethod_lut)<-c("genotype_calling","genotype_calling_perpop","hard_genotype_call")

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

################################################################################
# READ DATA




simdir="../../simulation/sim/"

GlMethod_i=c("gl1","gl2","precise1_gl2")
hdr<-c("Sample","nSitesTotal","nSitesRetained","nSitesCompared","nSitesCallMis","nDiscordantSites","nSitesInTrueNotCall","nConcordantSites","MissingnessRate","DiscordanceRate","ConcordanceRate","T_Hom","T_Het","F_HomToHom","F_HomToHet","F_HetToHom","F_HetToHet","T_Hom_rate","T_Het_rate","F_HomToHom_rate","F_HomToHet_rate","F_HetToHom_rate","F_HetToHet_rate","Rep","Depth","Betavar")

Method_i=c("genotype_calling", "genotype_calling_perpop", "hard_genotype_call")

# d<-NULL
# for(mi in seq_along(Method_i)){
# for(gli in seq_along(GlMethod_i)){
#   di<-read.csv(paste0(simdir, 
#   "sim_vcfgl_2312/model_OutOfAfrica_3G09/gc_evaluation/genotype_discordance_",GlMethod_i[gli],"/sim_vcfgl_2312-OutOfAfrica_3G09-",Method_i[mi],".tsv"),sep="\t",header=FALSE)
#   colnames(di)<-hdr
#   di$Gl<-GlMethod_i[gli]
#   di$GcMethod<-Method_i[mi]
#   d<-rbind(d,di)
# }
# }

# sd<-read.csv(paste0(simdir,"sim_vcfgl_2312/model_OutOfAfrica_3G09/stats/nSites_non0dp/sim_vcfgl_2312-OutOfAfrica_3G09_nSitesNon0DP.csv"),sep=",",header=TRUE)
# d<-merge(d,sd,by=c("Sample","Rep","Depth"))

# save(d,file=paste0("./data_",today,".RData"))


load(file=paste0("./data_",today,".RData"))

################################################################################

rm(d1,d2,d3,plt);gc();

d1<-d 
d1$Sample<-NULL
d1$Rep<-NULL


(d1%>% 
group_by(Gl,Betavar,Depth,GcMethod)%>%
summarize(
  nDiscordant=sum(nDiscordantSites),
  nConcordant=sum(nConcordantSites),
  nSitesRetained=sum(nSitesRetained),
  nSitesCompared=sum(nSitesCompared),
  nSitesCallMis=sum(nSitesCallMis),
  nSitesNon0DP=sum(nSitesNon0DP),
  nSitesInTrueNotCall=sum(nSitesInTrueNotCall),
  meanMissingnessRate=mean(MissingnessRate),
  meanDiscordanceRate=mean(DiscordanceRate),
  meanConcordanceRate=mean(ConcordanceRate)
)%>% 
ungroup()%>%
mutate(DiscordanceRate=nDiscordant/(nDiscordant+nConcordant),
        ConcordanceRate=nConcordant/(nDiscordant+nConcordant),
        MissingnessRate=nSitesCallMis/nSitesCompared,
        DiscordanceRateOverRetained=nDiscordant/nSitesRetained,
        ConcordanceRateOverRetained=nConcordant/nSitesRetained,
        MissingnessRateOverRetained=nSitesCallMis/nSitesRetained,
        DiscordanceRateOverCompared=nDiscordant/nSitesCompared,
        ConcordanceRateOverCompared=nConcordant/nSitesCompared,
        MissingnessRateOverCompared=nSitesCallMis/nSitesCompared,
        DiscordanceRateOverNon0DP=nDiscordant/nSitesNon0DP,
        MissingnessRateOverNon0DP=(nSitesNon0DP-(nDiscordant+nConcordant))/nSitesNon0DP,
        )%>%
mutate(across(c(Betavar,Depth,Gl,GcMethod),as.factor))%>% 
mutate(Betavar=factor(Betavar,levels=c("0","7","6","5"))))->d2




d2%>%
ggplot(aes(x=Depth,y=DiscordanceRate,color=Betavar,group=Betavar))+
  facet_wrap(Gl~GcMethod,scale="free",labeller=labeller1fn)+
  geom_point()+
  geom_line()+
  theme_bw()+
  theme(legend.position="top")+
  theme(strip.background = element_rect(colour="black", fill="white"))+
  #rotate x axis labels 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  color_brewer_betavars+
  labs(
    x="",
    y="Discordance rate",
    color="Beta variance",
    fill=""
  )+
  scale_color_manual(values=mypalette(5)[-1],labels=c("0"="0","7"=expression(10^-7), "6"=expression(10^-6),"5"=expression(10^-5)))

  


save_plt(plt=last_plot(),plot_filename="plot-pointLine_x-Depth_y-DiscordanceRate_color-Betavar_group-Betavar_facetwrap-GlGcMethod_data-allSitesAllGcMethods.png",plot_outdir=plot_outdir,overwrite=TRUE,width=12,height=10)


# d2%>%
# filter(Gl%in%c("gl1","gl2"))->d3
# d3%>%
# ggplot(aes(x=Gl,y=DiscordanceRate,color=Betavar,group=Betavar,fill=Betavar))+
#   facet_wrap(GcMethod~Depth,scale="free",labeller=labeller1fn,ncol=6)+
#   geom_point(pch=21,size=4,colour="black",position = position_dodge(width = 0.5))+
#   theme_bw()+
#   theme(legend.position="top")+
#   theme(strip.background = element_rect(colour="black", fill="white"))+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#   labs(
#     x="",
#     y="Discordance rate (%)",
#     color="Beta variance",
#     fill="Beta variance"
#   )+
#   scale_fill_manual(values=mypalette(5)[-1],labels=c("0"="0","7"=expression(10^-7), "6"=expression(10^-6),"5"=expression(10^-5)))+
#   scale_x_discrete(labels=c("gl1"="GL 1","gl2"="GL 2","precise1_gl2"="NA"))+
#   scale_y_continuous(labels = function(x) format(x*100, digits=8, nsmall=2))



# save_plt(plt=last_plot(),plot_filename="plot-point_x-Gl_y-DiscordanceRate_color-Betavar_group-Betavar_facetwrap-GcMethodDepth_data-allSitesAllGcMethods.png",plot_outdir=plot_outdir,overwrite=TRUE,width=12,height=10)




d2%>%
filter(Gl%in%c("gl1","gl2"))->d3


#get min max y values per depth 
#and set manual per-depth y axis limits for the plots
d3%>%
group_by(Depth)%>%
summarize(
  minDiscordanceRate=min(DiscordanceRate),
  maxDiscordanceRate=max(DiscordanceRate)
)-> d3_ylims



ylimsblank<-data.frame()
for(DPVAL in seq_along(levels(d3$Depth))){
  DPVALi<-d3_ylims$Depth[DPVAL]
  ylimsblanki<-data.frame(Depth=DPVALi,
  Gl=range(d3[d3$Depth%in%DPVALi,]$DiscordanceRate),
  DiscordanceRate=range(d3[d3$Depth%in%DPVALi,]$DiscordanceRate),
  Betavar=NA,
  stringsAsFactors = FALSE)
  ylimsblank<-rbind(ylimsblank,ylimsblanki)
}


betavar_expr_lut=c("0"="0","7"=expression(10^-7), "6"=expression(10^-6),"5"=expression(10^-5))


plotter1<-function(DT,DPVAL,GCMVAL,RMYLAB=FALSE,RMGRIDROWS=TRUE,RMXAXIS=FALSE){
DT%>%
filter(GcMethod%in%GCMVAL)%>%
filter(Depth%in%DPVAL)%>%
ggplot(aes(x=Gl,y=DiscordanceRate,group=Betavar,fill=Betavar))+
  geom_point(pch=21,size=4,colour="black")+
  # geom_blank(data=ylimsblank)+
  theme_bw()+
  theme(legend.position="top")+
  theme(strip.background = element_rect(colour="black", fill="white"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(
    x="",
    y="Discordance rate (%)",
    color="Beta variance",
    fill="Beta variance"
  )+
  scale_y_continuous(labels = function(x) format(x*100, digits=8, nsmall=2),
                      limits=c(d3_ylims[d3_ylims$Depth==DPVAL,]$minDiscordanceRate,d3_ylims[d3_ylims$Depth==DPVAL,]$maxDiscordanceRate))+
  scale_color_manual(values=mypalette(4),labels=betavar_expr_lut)+
  scale_fill_manual(values=mypalette(4),labels=betavar_expr_lut)+
  scale_x_discrete(labels=c("gl1"="GL 1","gl2"="GL 2","precise1_gl2"="NA"))+
  {if(RMYLAB)labs(y="")}+
  {if(RMGRIDROWS)theme(strip.text.y = element_blank())}+
  {if(RMXAXIS)theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())}
}



plist1<-NULL
plistpergcm<-NULL
for(gcmvali in seq_along(names(gcmethod_lut))){
  plisti<-NULL
for(dpvali in seq_along(levels(d3$Depth))){
  gcmval<-names(gcmethod_lut)[gcmvali]
  dpval<-levels(d3$Depth)[dpvali]
  # if(gcmvali==1){}
  if(gcmvali==2 & dpvali==1)rm_y_axis=FALSE else rm_y_axis=TRUE
  #rm x axis tick labels if not the last gcmval
  if(gcmvali!=3)rm_x_axis=TRUE else rm_x_axis=FALSE

#if first row, add facetting by depth 
  if(gcmvali==1){
    pi<-plotter1(d3,dpval,gcmval,RMYLAB=rm_y_axis,RMGRIDROWS=FALSE,RMXAXIS=rm_x_axis)+
       facet_wrap(.~Depth,scale="free",labeller=labeller1fn,ncol=6)
    pi<-pi+theme(plot.title = element_text(face="bold"))
  }else{
  pi<-plotter1(d3,dpval,gcmval,RMYLAB=rm_y_axis,RMGRIDROWS=FALSE,RMXAXIS=rm_x_axis)
  }
  #increase font size for all texts
  pi<-pi+theme(text = element_text(size=14))
  pi<-pi+theme(plot.margin = unit(c(0,0.1,0,0), "cm"))
  plisti<-c(plisti,list(pi))
  plist1<-c(plist1,list(pi))
}
  if(gcmvali==1){
    paneli<-ggarrange(plotlist=plisti,ncol=6,common.legend=TRUE,legend="top",align="hv")
  }else{
    # remove  legend from each plot in list
    plisti<-lapply(plisti, function(x) x + theme(legend.position="none")) 
    paneli<-ggarrange(plotlist=plisti,ncol=6)
  }
  paneli<-paneli+theme(plot.margin = unit(c(0,0,0,0), "cm"))
  plistpergcm<-c(plistpergcm,list(paneli))
}



ggarrange(plotlist=plistpergcm,ncol=1,nrow=3,common.legend=TRUE,legend="top",align="hv")

(ggarrange(plotlist=plistpergcm,
ncol=1,
nrow=3,
common.legend=TRUE,
align="hv",
legend="top",
labels = c("(A)","(B)","(C)")))->plt

annotate_figure(plt,
bottom="Genotype likelihood method")+bgcolor("white")


save_plt(plt=last_plot(),plot_filename="plot-panelABC-GcMethods-sharedyForDepths_geom-point_x-Gl_y-DiscordanceRate_fill-Betavar_group-Betavar_facetwrap-GcMethodDepth_data-allSitesAllGcMethods.png",plot_outdir=plot_outdir,overwrite=TRUE,width=12,height=10)


plotter2<-function(DT,GCMVAL,RMYLAB=FALSE,RMGRIDROWS=TRUE,RMGRIDCOLS=TRUE,RMXAXIS=FALSE){
DT%>%
filter(GcMethod==GCMVAL)%>%
ggplot(aes(x=Gl,y=DiscordanceRate,group=Betavar,fill=Betavar))+
  facet_wrap(.~Depth,scale="free",labeller=labeller1fn,ncol=6)+
  geom_point(pch=21,size=4,colour="black")+
  theme_bw()+
  theme(legend.position="top")+
  theme(strip.background = element_rect(colour="black", fill="white"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(
    x="",
    y="Discordance rate (%)",
    color="Beta variance",
    fill="Beta variance"
  )+
  scale_color_manual(values=mypalette(4),labels=betavar_expr_lut)+
  scale_fill_manual(values=mypalette(4),labels=betavar_expr_lut)+
  scale_x_discrete(labels=c("gl1"="GL 1","gl2"="GL 2","precise1_gl2"="NA"))+
  {if(RMYLAB)labs(y="")}+
  {if(RMGRIDROWS)theme(strip.text.y = element_blank())}+
  {if(RMGRIDCOLS)theme(strip.text.x = element_blank())}+
  {if(RMXAXIS)theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())}+
  scale_y_continuous(labels = function(x) format(x*100, digits=8, nsmall=2))
  # scale_y_continuous(labels = scales::percent_format(accuracy = 0.0001))
}


plist2<-NULL
for(gcmvali in seq_along(gcmethod_lut)){
  gcmval<-names(gcmethod_lut)[gcmvali]
  if(gcmvali==1)rm_grid_cols=FALSE else rm_grid_cols=TRUE
  #rm y axis label for 1 and 3
  if(gcmvali==1 | gcmvali==3)rm_y_axis=TRUE else rm_y_axis=FALSE
  #rm x axis tick labels if not the last gcmval
  if(gcmvali!=3)rm_x_axis=TRUE else rm_x_axis=FALSE

  pi<-plotter2(d3,gcmval,RMYLAB=rm_y_axis,RMGRIDCOLS=rm_grid_cols,RMGRIDROWS=FALSE, RMXAXIS=rm_x_axis  )+theme(text = element_text(size=14))+
  theme(plot.margin = unit(c(0,0.1,0,0), "cm"))
  plist2<-c(plist2,list(pi))
}

ggarrange(plotlist=plist2,ncol=1,nrow=3,common.legend=TRUE,legend="top",widths=c(100,0),align="hv", labels = c("(A)","(B)","(C)"))+bgcolor("white")
#add row labels

save_plt(plt=last_plot(),plot_filename="plot-panelABC-GcMethods_geom-point_x-Gl_y-DiscordanceRate_fill-Betavar_group-Betavar_facetwrap-GcMethodDepth_data-allSitesAllGcMethods.png",plot_outdir=plot_outdir,overwrite=TRUE,width=12,height=10)




# NOW COMPARE GLs


rm(d1,d2,d3,plt);gc();

d1<-d%>%filter(Betavar=="0")
# d1$Sample<-NULL
# d1$Rep<-NULL


(d1%>% 
group_by(Gl,Depth,GcMethod)%>%
summarize(
  nDiscordant=sum(nDiscordantSites),
  nConcordant=sum(nConcordantSites),
  nSitesRetained=sum(nSitesRetained),
  nSitesCompared=sum(nSitesCompared),
  nSitesCallMis=sum(nSitesCallMis),
  nSitesNon0DP=sum(nSitesNon0DP),
  nSitesInTrueNotCall=sum(nSitesInTrueNotCall),
  meanMissingnessRate=mean(MissingnessRate),
  meanDiscordanceRate=mean(DiscordanceRate),
  meanConcordanceRate=mean(ConcordanceRate)
)%>% 
ungroup()%>%
mutate(DiscordanceRate=nDiscordant/(nDiscordant+nConcordant),
        ConcordanceRate=nConcordant/(nDiscordant+nConcordant),
        MissingnessRate=nSitesCallMis/nSitesCompared,
        DiscordanceRateOverRetained=nDiscordant/nSitesRetained,
        ConcordanceRateOverRetained=nConcordant/nSitesRetained,
        MissingnessRateOverRetained=nSitesCallMis/nSitesRetained,
        DiscordanceRateOverCompared=nDiscordant/nSitesCompared,
        ConcordanceRateOverCompared=nConcordant/nSitesCompared,
        MissingnessRateOverCompared=nSitesCallMis/nSitesCompared,
        DiscordanceRateOverNon0DP=nDiscordant/nSitesNon0DP,
        MissingnessRateOverNon0DP=(nSitesNon0DP-(nDiscordant+nConcordant))/nSitesNon0DP,
        )%>%
mutate(across(c(Depth,Gl,GcMethod),as.factor)))->d2



d2%>%
filter(Gl%in%c("gl1","gl2"))->d3


cleanTheme1<-(theme_bw()+
  theme(legend.position="top")+
  theme(strip.background = element_rect(colour="black", fill="white"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)))




d3 %>%
mutate(across(c(Depth,Gl,GcMethod),as.factor))%>%
ggplot(aes(x=GcMethod,y=DiscordanceRate,fill=Gl))+
  geom_point(pch=21,size=4,colour="black",position = position_dodge(width = 0.5))+
  # geom_boxplot()+
  facet_wrap(
    .~Depth,
    scales="free",
    ncol=6,
    labeller=labeller2fn
  )+
  cleanTheme1





# plotter3<-function(DT,DPVAL,RMGRIDROWS=TRUE,RMGRIDCOLS=FALSE){
# DT %>%
# filter(Depth==DPVAL)%>%
# mutate(across(c(Depth,Gl,GcMethod),as.factor))%>%
# ggplot(aes(x="",y=DiscordanceRate,fill=Gl))+
#   # facet_wrap(
#   facet_grid(
#     rows=vars(GcMethod),
#     cols=vars(Depth),
#     # Depth~GcMethod,
#     scales="free",
#     labeller=labeller2fn
#     # ncol=1,
#     # nrow=3,
#   )+
#   cleanTheme1+
#   scale_y_continuous(labels = function(x) format(x*100, digits=8, nsmall=2))+
#   geom_bar(stat="identity",position = position_dodge(width = 0.5),width=0.1)+
#   geom_label(aes(label=glmethod_lut2[Gl],  hjust=1), position = position_dodge(width = 0.5),size=3,show.legend = FALSE)+
#   labs(
#     y="",
#     x="",
#     fill="Genotype likelihood model"
#   )+
#   guides(fill=guide_legend(title="Genotype likelihood model",title.hjust = 0.6))+
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
#   scale_fill_manual(values=c("chocolate1","darkturquoise"),
#   labels=c("gl1"="GL 1","gl2"="GL 2","precise1_gl2"="GL 2 with precise error"))+
#   {if(RMGRIDROWS)theme(strip.text.y = element_blank())}+
#   {if(RMGRIDCOLS)theme(strip.text.x = element_blank())}
# }


# plist3<-NULL
# for(DPVAL in seq_along(levels(d3$Depth))){
#   DPVALi<-levels(d3$Depth)[DPVAL]

#   pi<-plotter3(d3,DPVAL=DPVALi)
#   pi<-pi+theme(plot.title = element_text(face="bold"))
#   pi<-pi+theme(plot.margin = unit(c(0,0.1,0,0), "cm"))
#   pi<-pi+theme(text = element_text(size=14))
#   plist3<-c(plist3,list(pi))
# }

# plt<-(ggarrange(plotlist=plist3,common.legend=TRUE,legend="top",
# ncol=6,
# nrow=1))


# yposA<-0.97
# yposB<-0.65
# yposC<-(yposA-yposB)

# plt<-(plt+
# theme(plot.margin = unit(c(0.5,0,0,0.5), "cm"))+
# annotate("text", x = 0, y = yposA, label = "(A)",size=6,family="sans",fontface="bold")+
# annotate("text", x = 0, y = yposB, label = "(B)",size=6,family="sans",fontface="bold")+
# annotate("text", x = 0, y = yposC, label = "(C)",size=6,family="sans",fontface="bold"))


# annotate_figure(plt,
# left=text_grob("Discordance rate (%)", rot = 90, size=16))+bgcolor("white")


# save_plt(plt=last_plot(),plot_filename="plot-panelABC-GcMethods_geom-bar_x-Gl_y-PctDiscordanceRate_fill-Gl_facetgrid-GcMethodDepth_data-allSitesAllGcMethods.png",plot_outdir=plot_outdir,overwrite=TRUE,width=12,height=10)



plotter4<-function(DT,DPVAL,RMGRIDROWS=TRUE,RMGRIDCOLS=FALSE){
DT %>%
filter(Depth==DPVAL)%>%
mutate(across(c(Depth,Gl,GcMethod),as.factor))%>%
ggplot(aes(x="",y=DiscordanceRate,fill=Gl))+
  # facet_wrap(
  facet_grid(
    # GcMethod~Depth,
    rows=vars(GcMethod),
    cols=vars(Depth),
    # Depth~GcMethod,
    scale="fixed",
    labeller=labeller2fn
    # ncol=1,
    # nrow=3,
  )+
  cleanTheme1+
  scale_y_continuous(labels = function(x) format(x*100, digits=8, nsmall=2))+
  geom_bar(stat="identity",position = position_dodge(width = 0.5),width=0.1)+
  geom_label(aes(label=glmethod_lut2[Gl],  hjust=1), position = position_dodge(width = 0.5),size=3,show.legend = FALSE)+
  labs(
    y="",
    x="",
    fill="Genotype likelihood model"
  )+
  guides(fill=guide_legend(title="Genotype likelihood model",title.hjust = 0.6))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  scale_fill_manual(values=c("chocolate1","darkturquoise"),
  labels=c("gl1"="GL 1","gl2"="GL 2","precise1_gl2"="GL 2 with precise error"))+
  {if(RMGRIDROWS)theme(strip.text.y = element_blank())}+
  {if(RMGRIDCOLS)theme(strip.text.x = element_blank())}
}


plist4<-NULL
for(DPVAL in seq_along(levels(d3$Depth))){
  DPVALi<-levels(d3$Depth)[DPVAL]

  pi<-plotter4(d3,DPVAL=DPVALi)
  pi<-pi+theme(plot.title = element_text(face="bold"))
  pi<-pi+theme(plot.margin = unit(c(0,0.1,0,0), "cm"))
  pi<-pi+theme(text = element_text(size=14))
  plist4<-c(plist4,list(pi))
}

plt<-(ggarrange(plotlist=plist4,common.legend=TRUE,legend="top",
ncol=6,
nrow=1))


yposA<-0.97
yposB<-0.65
yposC<-(yposA-yposB)

plt<-(plt+
theme(plot.margin = unit(c(0.5,0,0,0.5), "cm"))+
annotate("text", x = 0, y = yposA, label = "(A)",size=6,family="sans",fontface="bold")+
annotate("text", x = 0, y = yposB, label = "(B)",size=6,family="sans",fontface="bold")+
annotate("text", x = 0, y = yposC, label = "(C)",size=6,family="sans",fontface="bold"))


annotate_figure(plt,
left=text_grob("Discordance rate (%)", rot = 90, size=16))+bgcolor("white")


save_plt(plt=last_plot(),plot_filename="plot-panelABC-GcMethods_geom-bar_x-Gl_y-PctDiscordanceRate_fill-Gl_facetgrid-GcMethodDepth_data-allSitesAllGcMethods.png",plot_outdir=plot_outdir,overwrite=TRUE,width=12,height=10)















################################################################################

