#rm(list=ls())
library(dplyr)
library(tidyr)




# ###############################
# # ~/.conda/envs/snakemake/bin/Rscript get_joint_common_tgt_discordance.R sim_v1/results/HG00096_chr21_depth100_gl1_subsamplev1_gt_tgt.tsv sim_v1/results/HG00096_chr21_depth100_gl1_subsamplev1_gt_tgt.tsv sim_v1/results/HG00096_chr21_depth100_rep16_gl1_vcfglv3_gt_tgt.tsv "subsamplev1" "vcfglv3" 16 100 1 FALSE sim_v1/results/HG00096_chr21_depth100_rep16_gl1_methods-subsamplev1-vcfglv3_usecommonv1_gt_tgt_discordance_persite.csv sim_v1/results/HG00096_chr21_depth100_rep16_gl1_methods-subsamplev1-vcfglv3_usecommonv1_gt_tgt_discordance_summary.csv v1 sim_v1/depth/HG00096_chr21_subsample_depth100_gl1_bcftools.bcf_avg_persite_depth
# args<-commandArgs(trailingOnly=TRUE)
# args[1]<-"sim_v1/results/HG00096_chr21_depth100_gl1_subsamplev1_gt_tgt.tsv"
# args[2]<-"sim_v1/results/HG00096_chr21_depth100_gl1_subsamplev1_gt_tgt.tsv"
# args[3]<-"sim_v1/results/HG00096_chr21_depth100_rep16_gl1_vcfglv3_gt_tgt.tsv"
# args[4]<-"subsamplev1"
# args[5]<-"vcfglv3"
# args[6]<-16
# args[7]<-100
# args[8]<-1
# args[9]<-FALSE
# args[10]<-"sim_v1/results/HG00096_chr21_depth100_rep16_gl1_methods-subsamplev1-vcfglv3_usecommonv1_gt_tgt_discordance_persite.csv"
# args[11]<-"sim_v1/results/HG00096_chr21_depth100_rep16_gl1_methods-subsamplev1-vcfglv3_usecommonv1_gt_tgt_discordance_summary.csv"
# args[12]<-"v1"
# args[13]<-"sim_v1/depth/HG00096_chr21_subsample_depth100_gl1_bcftools.bcf_avg_persite_depth"
# ###############################


##$                 ~/.conda/envs/snakemake/bin/Rscript get_joint_common_tgt_discordance.R sim_v4/results/HG00096_chr21_depth100_gt_tgt.tsv sim_v4/results/HG00096_chr21_depth100_gl1_subsamplev1_gt_tgt.tsv sim_v4/results/HG00096_chr21_depth100_rep19_gl1_vcfglv4_gt_tgt.tsv "subsamplev1" "vcfglv4" 19 100 1 FALSE sim_v4/results/HG00096_chr21_depth100_rep19_gl1_methods-subsamplev1-vcfglv4_usecommonv2_gt_tgt_discordance_persite.csv sim_v4/results/HG00096_chr21_depth100_rep19_gl1_methods-subsamplev1-vcfglv4_usecommonv2_gt_tgt_discordance_summary.csv v2 sim_v4/depth/HG00096_chr21_subsample_depth100_gl1_bcftools.bcf_avg_persite_depth 
#args[1]="sim_v4/results/HG00096_chr21_depth100_gt_tgt.tsv"
#args[2]="sim_v4/results/HG00096_chr21_depth100_gl1_subsamplev1_gt_tgt.tsv"
#args[3]="sim_v4/results/HG00096_chr21_depth100_rep19_gl1_vcfglv4_gt_tgt.tsv"
#args[4]="subsamplev1"
#args[5]="vcfglv4"
#args[6]=19
#args[7]=100
#args[8]=1
#args[9]=FALSE
#args[10]="sim_v4/results/HG00096_chr21_depth100_rep19_gl1_methods-subsamplev1-vcfglv4_usecommonv2_gt_tgt_discordance_persite.csv"
#args[11]="sim_v4/results/HG00096_chr21_depth100_rep19_gl1_methods-subsamplev1-vcfglv4_usecommonv2_gt_tgt_discordance_summary.csv"
#args[12]="v2"
#args[13]="sim_v4/depth/HG00096_chr21_subsample_depth100_gl1_bcftools.bcf_avg_persite_depth"




args<-commandArgs(trailingOnly=TRUE)
dref<-read.csv(args[1],sep='\t',header = FALSE)
dcall1<-read.csv(args[2],sep='\t',header = FALSE)
dcall2<-read.csv(args[3],sep='\t',header = FALSE)
type1name<-args[4]
type2name<-args[5]
dcall2_rep<-args[6]
target_depth<-args[7]
gl<-args[8]
allow_singlechar_gt<-as.logical(args[9])
outfile1<-args[10]
outfile2<-args[11]

use_common_version<-args[12]

actual_depth_file<-args[13]
actual_depth<-read.csv(actual_depth_file,sep='\t',header = FALSE)$V1







###############################
# # dummy data:
# dref<-data.frame(Chrom=c(rep("chr1",5),rep("chr2",5)),Pos=c(1:5,2:6),Genotype=c("AA","GG","GC","TT","AA","GG","CC","TT","AA","TC"))
# dcall1<-data.frame(Chrom=c(rep("chr1",4),rep("chr2",3)),Pos=c(1:4,4:6),Genotype=c("AA","GG","CG","TC","CC","TT","A"))
# dcall2<-data.frame(Chrom=c(rep("chr1",4),rep("chr2",2)),Pos=c(1:4,5:6),Genotype=c("AA","GG","CG","CC","AT","AA"))
# type1name<-"type1"
# type2name<-"type2"
# allow_singlechar_gt<-TRUE
# target_depth<-0.1
# actual_depth<-0.998
# gl<-1
# dcall2_rep<-"3"
###############################


colnames(dref)<-c("Chrom","Pos","Genotype")
colnames(dcall1)<-c("Chrom","Pos","Genotype")
colnames(dcall2)<-c("Chrom","Pos","Genotype")


merge(dref,dcall1,by=c("Chrom","Pos"),suffixes=c("_ref","_call1"),all.x=TRUE,all.y=FALSE) -> merged1
merge(dref,dcall2,by=c("Chrom","Pos"),suffixes=c("_ref","_call2"),all.x=TRUE,all.y=FALSE) -> merged2
merge(merged1,merged2,by=c("Chrom","Pos","Genotype_ref"),suffixes=c("_1","_2"),all.x=TRUE,all.y=TRUE) -> merged


# nothing in genotype_ref should be NA
if(sum(is.na(merged$Genotype_ref))!=0){
  stop("Genotype_ref contains NA")
}


#### ----------------------------------

if(use_common_version=="v1"){
  # include all sites
  print("usecommon v1: include all sites")
}

if(use_common_version=="v2"){
# usecommon v2: only use the sites that are present in type1 calls
# type1 is the subsampled truth data, only use sites where type1 calls is nonNA
  merged<-merged %>%
    filter(!is.na(Genotype_call1))
}

if(use_common_version=="v3"){
# usecommon v3: only use sites common for both type1 and type2 calls
  merged<-merged %>%
  filter(!is.na(Genotype_call1)) %>%
  filter(!is.na(Genotype_call2))
}

#### ----------------------------------

if(sum(nchar(merged$Genotype_ref)==1)!=0){
  stop("Genotype_ref contains single character genotypes")
}

if(sum(nchar(merged$Genotype_call1)==1,na.rm=TRUE)!=0){
  if(allow_singlechar_gt){
    merged$Genotype_call1<-ifelse(nchar(merged$Genotype_call1)==1,paste0(merged$Genotype_call1,merged$Genotype_call1),merged$Genotype_call1)
  }else{
    stop("Genotype_call contains single character genotypes")
  }
}

if(sum(nchar(merged$Genotype_call2)==1,na.rm=TRUE)!=0){
  if(allow_singlechar_gt){
    merged$Genotype_call2<-ifelse(nchar(merged$Genotype_call2)==1,paste0(merged$Genotype_call2,merged$Genotype_call2),merged$Genotype_call2)
  }else{
    stop("Genotype_call contains single character genotypes")
  }
}

get_unordered_edit_distance<-function(ref,call){
  if(is.na(call)){
    return(NA)
  }
  adist(paste(sort(unlist(strsplit(ref, ""))), collapse = ""),paste(sort(unlist(strsplit(call, ""))), collapse = ""))
}


# adist(c("AA","GG"),c("GG","AA"))
# apply so that only the elements at the same index are compared
merged$Distance1<-apply(merged[,c("Genotype_ref","Genotype_call1")],1,function(x) get_unordered_edit_distance(x[1],x[2]))
merged$Distance2<-apply(merged[,c("Genotype_ref","Genotype_call2")],1,function(x) get_unordered_edit_distance(x[1],x[2]))

merged$isDiscordant1<-as.integer(merged$Distance1>0)
merged$isDiscordant2<-as.integer(merged$Distance2>0)

write.table(merged,outfile1,sep=',',row.names=FALSE,col.names=TRUE,quote=FALSE)


merged$Reftype<-ifelse(unlist(lapply(strsplit(merged$Genotype_ref,""),function(x) length(unique(x))))==1,"Hom","Het")
merged$Reftype<-ifelse(is.na(merged$Genotype_ref),NA,merged$Reftype)
merged$Calltype1<-ifelse(unlist(lapply(strsplit(merged$Genotype_call1,""),function(x) length(unique(x))))==1,"Hom","Het")
merged$Calltype1<-ifelse(is.na(merged$Genotype_call1),NA,merged$Calltype1)
merged$Calltype2<-ifelse(unlist(lapply(strsplit(merged$Genotype_call2,""),function(x) length(unique(x))))==1,"Hom","Het")
merged$Calltype2<-ifelse(is.na(merged$Genotype_call2),NA,merged$Calltype2)

ifelse(merged$isDiscordant1==1, paste0(merged$Reftype,"To",merged$Calltype1,"Discordant"),paste0(merged$Reftype,"To",merged$Calltype1,"Concordant")) -> merged$Type1
ifelse(merged$isDiscordant2==1, paste0(merged$Reftype,"To",merged$Calltype2,"Discordant"),paste0(merged$Reftype,"To",merged$Calltype2,"Concordant")) -> merged$Type2

# nCategories:
# Het Het Concordant
# Het Het Discordant
# Hom Hom Concordant
# Hom Hom Discordant
# Het Hom Discordant
# Hom Het Discordant

Categories=c(
  "nSitesHetInCall",
  "nSitesHomInCall",
  "DiscordanceRateOverCalls",
  "DiscordanceRateOverAll",
  "MissingnessRate",
  "nDiscordantSites",
  "nConcordantSites",
  "nSitesCompared",
  "nSitesInRefNotinCall",
  "nSitesInCallNotinRef",
  "nSitesInCall",
  "nHetToHetConcordant",
  "nHetToHetDiscordant",
  "nHomToHomConcordant",
  "nHomToHomDiscordant",
  "nHetToHomDiscordant",
  "nHomToHetDiscordant"
)


# shared for both types
summary0<-data.frame(
  Category=c("nSitesHetInRef","nSitesHomInRef","nSitesTotal","nSitesInRef"),
  Type=NA,
  Value=c(
    sum(merged$Reftype=="Het",na.rm=TRUE),
    sum(merged$Reftype=="Hom",na.rm=TRUE),
    nrow(merged),
    sum(!is.na(merged$Genotype_ref))
  ),
  Rep=NA,
  TargetDepth=target_depth,
  ActualDepth=actual_depth,
  GL=gl
)



summary0$Type<-type1name
summary0%>% rbind(
# type 1 
data.frame(
  Category=Categories,
  Type=type1name,
  Value=c(
    sum(merged$Calltype1=="Het",na.rm=TRUE),
    sum(merged$Calltype1=="Hom",na.rm=TRUE),
    sum(merged$isDiscordant1,na.rm=TRUE)/(sum(!merged$isDiscordant1,na.rm=TRUE)+sum(merged$isDiscordant1,na.rm=TRUE)),
    sum(merged$isDiscordant1,na.rm=TRUE)/nrow(merged),
    sum(is.na(merged$Genotype_call1))/nrow(merged),
    sum(merged$isDiscordant1,na.rm=TRUE),
    sum(!merged$isDiscordant1,na.rm=TRUE),
    sum(!is.na(merged$Genotype_ref) & !is.na(merged$Genotype_call1)),
    sum(!is.na(merged$Genotype_ref) & is.na(merged$Genotype_call1)),
    sum(is.na(merged$Genotype_ref) & !is.na(merged$Genotype_call1)),
    sum(!is.na(merged$Genotype_call1)),
    sum(merged$Type1=="HetToHetConcordant",na.rm=TRUE),
    sum(merged$Type1=="HetToHetDiscordant",na.rm=TRUE),
    sum(merged$Type1=="HomToHomConcordant",na.rm=TRUE),
    sum(merged$Type1=="HomToHomDiscordant",na.rm=TRUE),
    sum(merged$Type1=="HetToHomDiscordant",na.rm=TRUE),
    sum(merged$Type1=="HomToHetDiscordant",na.rm=TRUE)
  ),
  Rep=NA,
  TargetDepth=target_depth,
  ActualDepth=actual_depth,
  GL=gl
)
) -> summary1

summary0$Type<-type2name
summary0%>% rbind(
# type 2
data.frame(
  Category=Categories,
  Type=type2name,
  Value=c(
    sum(merged$Calltype2=="Het",na.rm=TRUE),
    sum(merged$Calltype2=="Hom",na.rm=TRUE),
    sum(merged$isDiscordant2,na.rm=TRUE)/(sum(!merged$isDiscordant2,na.rm=TRUE)+sum(merged$isDiscordant2,na.rm=TRUE)),
    sum(merged$isDiscordant2,na.rm=TRUE)/nrow(merged),
    sum(is.na(merged$Genotype_call2))/nrow(merged),
    sum(merged$isDiscordant2,na.rm=TRUE),
    sum(!merged$isDiscordant2,na.rm=TRUE),
    sum(!is.na(merged$Genotype_ref) & !is.na(merged$Genotype_call2)),
    sum(!is.na(merged$Genotype_ref) & is.na(merged$Genotype_call2)),
    sum(is.na(merged$Genotype_ref) & !is.na(merged$Genotype_call2)),
    sum(!is.na(merged$Genotype_call2)),
    sum(merged$Type2=="HetToHetConcordant",na.rm=TRUE),
    sum(merged$Type2=="HetToHetDiscordant",na.rm=TRUE),
    sum(merged$Type2=="HomToHomConcordant",na.rm=TRUE),
    sum(merged$Type2=="HomToHomDiscordant",na.rm=TRUE),
    sum(merged$Type2=="HetToHomDiscordant",na.rm=TRUE),
    sum(merged$Type2=="HomToHetDiscordant",na.rm=TRUE)
  ),
  Rep=dcall2_rep,
  TargetDepth=target_depth,
  ActualDepth=actual_depth,
  GL=gl
  )
) -> summary2

summary<-rbind(summary1,summary2)

write.table(summary,outfile2,sep=',',row.names=FALSE,col.names=FALSE,quote=FALSE)

#########################################################################