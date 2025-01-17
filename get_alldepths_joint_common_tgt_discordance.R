# USE COMMON SITES FROM MULTIPLE DEPTHS
rm(list=ls())
library(dplyr)
library(tidyr)

DEPTHS = c(1, 5, 10, 20, 100)

args<-commandArgs(trailingOnly=TRUE)

# get_unordered_edit_distance<-function(ref,call){
#   # if(is.na(call)){
#   #   return(NA)
#   # }
#   # print(ref)
#   # stop("A")
#   # adist(paste(sort(unlist(strsplit(ref, ""))), collapse = ""),paste(sort(unlist(strsplit(call, ""))), collapse = ""))
#   # diag(adist(sapply(strsplit(ref,""),function(x) paste(sort(x),collapse="")), sapply(strsplit(call,""),function(x) paste(sort(x),collapse=""))))
#   # ref<-sapply(strsplit(ref,""),function(x) paste(sort(x),collapse=""))
#   # ref[ref==""]<-NA
#   # call<-sapply(strsplit(call,""),function(x) paste(sort(x),collapse=""))
#   # call[call==""]<-NA
  
#   diag(adist(sapply(strsplit(ref,""),function(x){
#     # if(length(x)==0){
#     if(sum(is.na(x))>0){
#       return(NA)
#     }else{
#       return(paste(sort(x),collapse=""))
#     }
#   }), sapply(strsplit(call,""),function(x){
#     # print(sum(is.na(x)))
#     # print(x)
#     if(sum(is.na(x))>0){
#       return(NA)
#     }else{
#       return(paste(sort(x),collapse=""))
#     }
#   })))
# }


get_unordered_edit_distance<-function(ref,call){
  if(is.na(call)){
    return(NA)
  }
  adist(paste(sort(unlist(strsplit(ref, ""))), collapse = ""),paste(sort(unlist(strsplit(call, ""))), collapse = ""))
}


# x1=c("AA","AC","AA","AA","GG","TT","CT")
# x2=c("AA","AA", NA, "GT","GG", NA, "TC")
# get_unordered_edit_distance(x1,x2)

# ###############################
# # ~/.conda/envs/snakemake/bin/Rscript get_alldepths_joint_common_tgt_discordance.R sim_v1/results/HG00096_chr21_depth100_gl1_subsamplev1_gt_tgt.tsv sim_v1/results/HG00096_chr21_depth1_gl1_subsamplev1_gt_tgt.tsv sim_v1/results/HG00096_chr21_depth2_gl1_subsamplev1_gt_tgt.tsv sim_v1/results/HG00096_chr21_depth10_gl1_subsamplev1_gt_tgt.tsv sim_v1/results/HG00096_chr21_depth20_gl1_subsamplev1_gt_tgt.tsv sim_v1/results/HG00096_chr21_depth100_gl1_subsamplev1_gt_tgt.tsv sim_v1/results/HG00096_chr21_depth1_rep6_gl1_vcfglv3_gt_tgt.tsv sim_v1/results/HG00096_chr21_depth2_rep6_gl1_vcfglv3_gt_tgt.tsv sim_v1/results/HG00096_chr21_depth10_rep6_gl1_vcfglv3_gt_tgt.tsv sim_v1/results/HG00096_chr21_depth20_rep6_gl1_vcfglv3_gt_tgt.tsv sim_v1/results/HG00096_chr21_depth100_rep6_gl1_vcfglv3_gt_tgt.tsv 6 1 FALSE sim_v1/results/HG00096_chr21_alldepths_rep6_gl1_methods-subsamplev1-vcfglv3_usecommonv1_gt_tgt_discordance_persite.csv sim_v1/results/HG00096_chr21_alldepths_rep6_gl1_methods-subsamplev1-vcfglv3_usecommonv1_gt_tgt_discordance_summary.csv v1 sim_v1/depth/HG00096_chr21_subsample_depth1_gl1_bcftools.bcf_avg_persite_depth sim_v1/depth/HG00096_chr21_subsample_depth2_gl1_bcftools.bcf_avg_persite_depth sim_v1/depth/HG00096_chr21_subsample_depth10_gl1_bcftools.bcf_avg_persite_depth sim_v1/depth/HG00096_chr21_subsample_depth20_gl1_bcftools.bcf_avg_persite_depth sim_v1/depth/HG00096_chr21_subsample_depth100_gl1_bcftools.bcf_avg_persite_depth
# args[1]<-"sim_v1/results/HG00096_chr21_depth100_gl1_subsamplev1_gt_tgt.tsv"
# args[2]<-"sim_v1/results/HG00096_chr21_depth1_gl1_subsamplev1_gt_tgt.tsv"
# args[3]<-"sim_v1/results/HG00096_chr21_depth2_gl1_subsamplev1_gt_tgt.tsv"
# args[4]<-"sim_v1/results/HG00096_chr21_depth10_gl1_subsamplev1_gt_tgt.tsv"
# args[5]<-"sim_v1/results/HG00096_chr21_depth20_gl1_subsamplev1_gt_tgt.tsv"
# args[6]<-"sim_v1/results/HG00096_chr21_depth100_gl1_subsamplev1_gt_tgt.tsv"
# args[7]<-"sim_v1/results/HG00096_chr21_depth1_rep6_gl1_vcfglv3_gt_tgt.tsv"
# args[8]<-"sim_v1/results/HG00096_chr21_depth2_rep6_gl1_vcfglv3_gt_tgt.tsv"
# args[9]<-"sim_v1/results/HG00096_chr21_depth10_rep6_gl1_vcfglv3_gt_tgt.tsv"
# args[10]<-"sim_v1/results/HG00096_chr21_depth20_rep6_gl1_vcfglv3_gt_tgt.tsv"
# args[11]<-"sim_v1/results/HG00096_chr21_depth100_rep6_gl1_vcfglv3_gt_tgt.tsv"
# args[12]<-6
# args[13]<-1
# args[14]<-FALSE
# args[15]<-"sim_v1/results/HG00096_chr21_alldepths_rep6_gl1_methods-subsamplev1-vcfglv3_usecommonv1_gt_tgt_discordance_persite.csv"
# args[16]<-"sim_v1/results/HG00096_chr21_alldepths_rep6_gl1_methods-subsamplev1-vcfglv3_usecommonv1_gt_tgt_discordance_summary.csv"
# args[17]<-"v1"
# args[18]<-"sim_v1/depth/HG00096_chr21_subsample_depth1_gl1_bcftools.bcf_avg_persite_depth"
# args[19]<-"sim_v1/depth/HG00096_chr21_subsample_depth2_gl1_bcftools.bcf_avg_persite_depth"
# args[20]<-"sim_v1/depth/HG00096_chr21_subsample_depth10_gl1_bcftools.bcf_avg_persite_depth"
# args[21]<-"sim_v1/depth/HG00096_chr21_subsample_depth20_gl1_bcftools.bcf_avg_persite_depth"
# args[22]<-"sim_v1/depth/HG00096_chr21_subsample_depth100_gl1_bcftools.bcf_avg_persite_depth"





# ###############################
# # # dummy data:
# # dref<-data.frame(Chrom=c(rep("chr1",5),rep("chr2",5)),Pos=c(1:5,2:6),Genotype=c("AA","GG","GC","TT","AA","GG","CC","TT","AA","TC"))
# # dcall1<-data.frame(Chrom=c(rep("chr1",4),rep("chr2",3)),Pos=c(1:4,4:6),Genotype=c("AA","GG","CG","TC","CC","TT","A"))
# # dcall2<-data.frame(Chrom=c(rep("chr1",4),rep("chr2",2)),Pos=c(1:4,5:6),Genotype=c("AA","GG","CG","CC","AT","AA"))
# merged<-data.frame(Chrom=c(rep("chr1",5),rep("chr2",5)),Pos=c(1:5,2:6),
# Genotype_ref=c("AA","GG","GC","TT","AA","GG","CC","TT","AA","TC"),
# Genotype_subsamplev1_1=c("AA","GG","GC","TT","AA","GG","CC","TT","AA","TC"),
# Genotype_subsamplev1_2=c("AA","GG","GC","TT","AA","GG","CC","TT","AA","CT"),
# Genotype_subsamplev1_10=c("AA","GG","GC","TT","AA","GG","CC","TT","AA","TC"),
# Genotype_subsamplev1_20=c("AA","GG","GC","TT","AA","GG","CC","TT","AA","TC"),
# Genotype_subsamplev1_100=c("AA","GG","GC","TT","AA","GG","CC","TT","AA","TC"),
# Genotype_vcfglv3_1=c("AA","GG","GC","TT","AA","GG","CC","TT","AA","TC"),
# Genotype_vcfglv3_2=c("AA","GG","GC","TT","AA","GG","CC","TT","AA","TC"),
# Genotype_vcfglv3_10=c("AA","GG","GC","TT","AA","GG","CC","TT","AA","TC"),
# Genotype_vcfglv3_20=c("AA","GG","GC","TT","AA","GG","CC","TC","CC","TC"),
# Genotype_vcfglv3_100=c("AA","GG","GC","TT","AA","GG","CC","TT","AA","CC"))
# use_common_version="v1"
# type1name<-"type1"
# type2name<-"type2"
# allow_singlechar_gt<-TRUE
# # target_depth<-0.1
# # actual_depth<-0.998
# gl<-1
# dcall2_rep<-"3"



merged<-read.csv(args[1],sep='\t',header = FALSE)

colnames(merged)<-c("Chrom","Pos","Genotype_ref")


for (i in 2:6){

  fname<-args[i]
  typename<-strsplit(basename(fname),"_")[[1]][length(strsplit(basename(fname),"_")[[1]])-2]
  dcall1i<-read.csv(fname,sep='\t',header = FALSE)
# get the depth from filename
#"sim_v1/results/HG00096_chr21_depth{depth}_gl{{gl}}_subsamplev1_gt_tgt.tsv"
  depth<-as.numeric(sub('.*depth([0-9]+).*', '\\1', fname))
  if(depth!=DEPTHS[i-1]){
    stop("depth mismatch")
  }
  colnames(dcall1i)<-c("Chrom","Pos",paste0("Genotype_",typename,"_",depth))
  merged<-merge(merged,dcall1i,by=c("Chrom","Pos"),all.x=TRUE,all.y=TRUE)
  apply(merged[,c("Genotype_ref",paste0("Genotype_",typename,"_",depth))],1,function(x) get_unordered_edit_distance(x[1],x[2])) -> merged[[paste0("Distance_",typename,"_",depth)]]
}

for(i in 7:11){
  fname<-args[i]
  typename<-strsplit(basename(fname),"_")[[1]][length(strsplit(basename(fname),"_")[[1]])-2]
  dcall2i<-read.csv(fname,sep='\t',header = FALSE)
  # "sim_v1/results/HG00096_chr21_depth{depth}_rep{{rep}}_gl{{gl}}_vcfgl{{vcfglv}}_gt_tgt.tsv"
  depth<-as.numeric(sub('.*depth([0-9]+).*', '\\1', fname))
  if(depth!=DEPTHS[i-6]){
    stop("depth mismatch")
  }
  colnames(dcall2i)<-c("Chrom","Pos",paste0("Genotype_",typename,"_",depth))
  merged<-merge(merged,dcall2i,by=c("Chrom","Pos"),all.x=TRUE,all.y=TRUE)
  apply(merged[,c("Genotype_ref",paste0("Genotype_",typename,"_",depth))],1,function(x) get_unordered_edit_distance(x[1],x[2])) -> merged[[paste0("Distance_",typename,"_",depth)]]
}

dcall2_rep<-args[12]
gl<-args[13]
allow_singlechar_gt<-as.logical(args[14])
outfile1<-args[15]
outfile2<-args[16]

use_common_version<-args[17]

actual_depths<-data.frame(Depth=NULL,ActualDepth=NULL)
for(i in 18:length(args)){
  fname<-args[i]
  actual_depthi<-read.csv(fname,sep='\t',header = FALSE)
  depth<-as.numeric(sub('.*depth([0-9]+).*', '\\1', fname))
  if(depth!=DEPTHS[i-17]){
    stop("depth mismatch")
  }
  actual_depths<-rbind(actual_depths,data.frame(Depth=depth,ActualDepth=actual_depthi$V1))
}




# actual_depths<-data.frame(Depth=c(1,2,10,20,100),ActualDepth=c(1.1,2.2,10.1,20.01,100.2))
# ###############################










# nothing in genotype_ref should be NA
if(sum(is.na(merged$Genotype_ref))!=0){
  stop("Genotype_ref contains NA")
}



#### ----------------------------------


if(use_common_version=="v1"){
  print("usecommon v1: use sites with nonNA in all Genotypes")
  # no NA in any column with name starting with Genotype_*
  merged<-merged[rowSums(is.na(merged[,grep("^Genotype_",colnames(merged))]))==0,]
}


#### ----------------------------------

if(sum(nchar(merged$Genotype_ref)==1)!=0){
  stop("Genotype_ref contains single character genotypes")
}

if(!allow_singlechar_gt){
  for(i in 1:ncol(merged)){
    if(grepl("^Genotype_",colnames(merged)[i])){
      if(sum(nchar(merged[,i])==1,na.rm=TRUE)!=0){
        stop(paste0("Genotype_",i," contains single character genotypes"))
      }
    }
  }

}




for(i in 1:ncol(merged)){
  if(grepl("^Distance_",colnames(merged)[i])){
    #get type and depth from column name
    typei<-strsplit(colnames(merged)[i],"_")[[1]][2]
    if(typei=="ref"){
      next
    }
    depthi<-strsplit(colnames(merged)[i],"_")[[1]][3]
    merged[[paste0("isDiscordant_",typei,"_",depthi)]]<-as.integer(merged[[colnames(merged)[i]]]>0)
  }
}

write.table(merged,outfile1,sep=',',row.names=FALSE,col.names=TRUE,quote=FALSE)


merged$Reftype<-ifelse(unlist(lapply(strsplit(merged$Genotype_ref,""),function(x) length(unique(x))))==1,"Hom","Het")
merged$Reftype<-ifelse(is.na(merged$Genotype_ref),NA,merged$Reftype)


for(i in 1:ncol(merged)){
  if(grepl("^Genotype_",colnames(merged)[i])){
    typei<-strsplit(colnames(merged)[i],"_")[[1]][2]
    if(typei=="ref"){
      next
    }
    depthi<-strsplit(colnames(merged)[i],"_")[[1]][3]
    merged[[paste0("Calltype_",typei,"_",depthi)]]<-ifelse(unlist(lapply(strsplit(merged[,i],""),function(x) length(unique(x))))==1,"Hom","Het")
    merged[[paste0("Calltype_",typei,"_",depthi)]]<-ifelse(is.na(merged[,i]),NA,merged[[paste0("Calltype_",typei,"_",depthi)]])
    merged[[paste0("Type_",typei,"_",depthi)]]<-paste0(merged$Reftype,"To",merged[[paste0("Calltype_",typei,"_",depthi)]],ifelse(merged[[paste0("isDiscordant_",typei,"_",depthi)]]==1,"Discordant","Concordant"))
  }
}


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
  TargetDepth=NA,
  ActualDepth=NA,
  GL=gl
)


#TODO HERE
summary<-NULL
for(i in 1:ncol(merged)){

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
    TargetDepth=NA,
    ActualDepth=NA,
    GL=gl
  )

  if(grepl("^Genotype_",colnames(merged)[i])){
    typei<-strsplit(colnames(merged)[i],"_")[[1]][2]
    if(typei=="ref"){
      next
    }
    if(grepl("^vcfgl",typei)){
      repi<-dcall2_rep
    }else{
      repi<-NA
    }
    depthi<-strsplit(colnames(merged)[i],"_")[[1]][3]
    summary0$Type<-typei
    summary0$TargetDepth<-depthi
    summary0$ActualDepth<-actual_depths[actual_depths$Depth==as.numeric(depthi),"ActualDepth"]

    summary0%>% rbind(
      data.frame(
        Category=Categories,
        Type=typei,
        Value=c(
          sum(merged[[paste0("Calltype_",typei,"_",depthi)]]=="Het",na.rm=TRUE),
          sum(merged[[paste0("Calltype_",typei,"_",depthi)]]=="Hom",na.rm=TRUE),
          sum(merged[[paste0("isDiscordant_",typei,"_",depthi)]],na.rm=TRUE)/(sum(!merged[[paste0("isDiscordant_",typei,"_",depthi)]],na.rm=TRUE)+sum(merged[[paste0("isDiscordant_",typei,"_",depthi)]],na.rm=TRUE)),
          sum(merged[[paste0("isDiscordant_",typei,"_",depthi)]],na.rm=TRUE)/nrow(merged),
          sum(is.na(merged[[i]])/nrow(merged)),
          sum(merged[[paste0("isDiscordant_",typei,"_",depthi)]],na.rm=TRUE),
          sum(!merged[[paste0("isDiscordant_",typei,"_",depthi)]],na.rm=TRUE),
          sum(!is.na(merged$Genotype_ref) & !is.na(merged[[i]])),
          sum(!is.na(merged$Genotype_ref) & is.na(merged[[i]])),
          sum(is.na(merged$Genotype_ref) & !is.na(merged[[i]])),
          sum(!is.na(merged[[i]])),
          sum(merged[[paste0("Type_",typei,"_",depthi)]]=="HetToHetConcordant",na.rm=TRUE),
          sum(merged[[paste0("Type_",typei,"_",depthi)]]=="HetToHetDiscordant",na.rm=TRUE),
          sum(merged[[paste0("Type_",typei,"_",depthi)]]=="HomToHomConcordant",na.rm=TRUE),
          sum(merged[[paste0("Type_",typei,"_",depthi)]]=="HomToHomDiscordant",na.rm=TRUE),
          sum(merged[[paste0("Type_",typei,"_",depthi)]]=="HetToHomDiscordant",na.rm=TRUE),
          sum(merged[[paste0("Type_",typei,"_",depthi)]]=="HomToHetDiscordant",na.rm=TRUE)
        ),
        Rep=repi,
        TargetDepth=depthi,
        ActualDepth=actual_depths[actual_depths$Depth==as.numeric(depthi),"ActualDepth"],
        GL=gl
      )
    ) -> summaryi
    summary <- rbind(summary,summaryi)
  }
}
write.table(summary,outfile2,sep=',',row.names=FALSE,col.names=FALSE,quote=FALSE)

#########################################################################