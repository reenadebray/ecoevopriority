## Code for analyzing evolved bacterial populations
# Reena Debray
# 2022


# Identify sites that differ from reference in evolved populations but not in ancestor:
## remove sites present in ancestor AND all evolved pops
t<-table(EEP_breseq$Position)
sites_all<-names(t[t==7])
EEP_breseq<-EEP_breseq[!EEP_breseq$Position%in%sites_all,]


## variants present in ancestor
RIFANC_sites<-unlist(EEP_breseq[EEP_breseq$Sample=="RIFANC","Position"])
## remove variants present in ancestor
EEP_breseq_excludeANC<-EEP_breseq[!EEP_breseq$Position%in%RIFANC_sites,]

## to avoid false positives from highly polymorphic regions, filter variant calls with 5 or more sites within a 50-bp sliding window
EEP_breseq_excludeANC_filtered<-data.frame()
N=4
for (sample in unique(EEP_breseq_excludeANC$Sample)){
  pre_filter<-EEP_breseq_excludeANC[EEP_breseq_excludeANC$Sample==sample,]
  pre_filter_POS<-sort(pre_filter$Position) 
  sites_to_remove<-c()
  for (pos in pre_filter_POS){
    ## sliding window
    min=pos-50
    while (min<pos){
      max<-min+50
      if (length(pre_filter_POS[pre_filter_POS>=min & pre_filter_POS<max])>N) {sites_to_remove<-unique(c(sites_to_remove,pre_filter_POS[pre_filter_POS>=min & pre_filter_POS<max]))} # if a window contains more than 3 mutations, remove all of them
      min<-min+1
    }
  }
  sites_to_remove<-unique(sites_to_remove)
  ## filter sites
  post_filter<-EEP_breseq_excludeANC[EEP_breseq_excludeANC$Sample==sample & !EEP_breseq_excludeANC$Position%in%sites_to_remove,]
  EEP_breseq_excludeANC_filtered<-rbind(EEP_breseq_excludeANC_filtered,post_filter)
}

## filter for coding mutations
EEP_breseq_excludeANC_coding<-EEP_breseq_excludeANC_filtered[!grepl("intergenic",EEP_breseq_excludeANC_filtered$Annotation,fixed=TRUE) & !is.na(EEP_breseq_excludeANC_filtered$Annotation),]
for (i in seq(1,nrow(EEP_breseq_excludeANC_coding))){
  amino<-unlist(strsplit(unlist(EEP_breseq_excludeANC_coding[i,"Annotation"]),split="[(]"))[1]
  ref<-substr(amino,1,1)
  alt<-substr(amino,nchar(amino)-1,nchar(amino)-1)
  if (ref==alt & !is.na(ref) & !is.na(alt)){EEP_breseq_excludeANC_coding<-EEP_breseq_excludeANC_coding[-c(i),]}
}

## subset to evolved populations and calculate metapopulation statistics (mean allele frequency, occurrence across populations)
EEP_breseq_evo<-EEP_breseq_excludeANC_coding
for (i in seq(1,nrow(EEP_breseq_evo))){
  pos<-unlist(EEP_breseq_evo[i,"Position"])
  EEP_breseq_evo[i,"pop_freq"]<-nrow(EEP_breseq_evo[EEP_breseq_evo$Position==pos,])
  EEP_breseq_evo[i,"mean_AF"]<-mean(unlist(EEP_breseq_evo[EEP_breseq_evo$Position==pos,"Frequency"]))
}


# Identify sites that differ from reference in ancestor but not in evolved populations:
## to avoid false positives from highly polymorphic regions, filter variant calls with 5 or more sites within a 50-bp sliding window
RIFANC<-EEP_breseq[EEP_breseq$Sample=="RIFANC",]
N=4
pre_filter<-RIFANC
pre_filter_POS<-sort(pre_filter$Position) 
sites_to_remove<-c()
for (pos in pre_filter_POS){
  ## sliding window
  min=pos-50
  while (min<pos){
    max<-min+50
    if (length(pre_filter_POS[pre_filter_POS>=min & pre_filter_POS<max])>N) {sites_to_remove<-unique(c(sites_to_remove,pre_filter_POS[pre_filter_POS>=min & pre_filter_POS<max]))} # if a window contains more than 3 mutations, remove all of them
    min<-min+1
  }
}
sites_to_remove<-unique(sites_to_remove)
## filter sites
ANC_filtered<-RIFANC[!RIFANC$Position%in%sites_to_remove,]

## filter for fixed sites and coding mutations
ANC_fixed<-ANC_filtered[ANC_filtered$Frequency==1,]
ANC_coding<-ANC_fixed[!grepl("intergenic",ANC_fixed$Annotation,fixed=TRUE) & !is.na(ANC_fixed$Annotation),]
for (i in seq(1,nrow(ANC_coding))){
  amino<-unlist(strsplit(unlist(ANC_coding[i,"Annotation"]),split="[(]"))[1]
  ref<-substr(amino,1,1)
  alt<-substr(amino,nchar(amino)-1,nchar(amino)-1)
  if (ref==alt & !is.na(ref) & !is.na(alt)){ANC_coding<-ANC_coding[-c(i),]}
}

## identify mutations not in evolved populations and add to data frame
ANC_coding<-ANC_coding[,1:7]
RIFEVO1_rev_sites<-setdiff(ANC_coding$Position,unlist(EEP_breseq[EEP_breseq$Sample=="RIFEVO1","Position"]))
for (pos in RIFEVO1_rev_sites){
  info<-ANC_coding[ANC_coding$Position==pos,]
  info[1]<-"rev_ref"
  EEP_breseq_evo[nrow(EEP_breseq_evo)+1,1:7]<-info
  EEP_breseq_evo[nrow(EEP_breseq_evo),8]<-"RIFEVO1"
}
RIFEVO4_rev_sites<-setdiff(ANC_coding$Position,unlist(EEP_breseq[EEP_breseq$Sample=="RIFEVO4","Position"]))
for (pos in RIFEVO1_rev_sites){
  info<-ANC_coding[ANC_coding$Position==pos,]
  info[1]<-"rev_ref"
  EEP_breseq_evo[nrow(EEP_breseq_evo)+1,1:7]<-info
  EEP_breseq_evo[nrow(EEP_breseq_evo),8]<-"RIFEVO4"
}
RIFEVO10_rev_sites<-setdiff(ANC_coding$Position,unlist(EEP_breseq[EEP_breseq$Sample=="RIFEVO10","Position"]))
for (pos in RIFEVO1_rev_sites){
  info<-ANC_coding[ANC_coding$Position==pos,]
  info[1]<-"rev_ref"
  EEP_breseq_evo[nrow(EEP_breseq_evo)+1,1:7]<-info
  EEP_breseq_evo[nrow(EEP_breseq_evo),8]<-"RIFEVO10"
}
RIFEVO11_rev_sites<-setdiff(ANC_coding$Position,unlist(EEP_breseq[EEP_breseq$Sample=="RIFEVO11","Position"]))
for (pos in RIFEVO1_rev_sites){
  info<-ANC_coding[ANC_coding$Position==pos,]
  info[1]<-"rev_ref"
  EEP_breseq_evo[nrow(EEP_breseq_evo)+1,1:7]<-info
  EEP_breseq_evo[nrow(EEP_breseq_evo),8]<-"RIFEVO11"
}
RIFEVO15_rev_sites<-setdiff(ANC_coding$Position,unlist(EEP_breseq[EEP_breseq$Sample=="RIFEVO15","Position"]))
for (pos in RIFEVO1_rev_sites){
  info<-ANC_coding[ANC_coding$Position==pos,]
  info[1]<-"rev_ref"
  EEP_breseq_evo[nrow(EEP_breseq_evo)+1,1:7]<-info
  EEP_breseq_evo[nrow(EEP_breseq_evo),8]<-"RIFEVO15"
}
RIFEVO23_rev_sites<-setdiff(ANC_coding$Position,unlist(EEP_breseq[EEP_breseq$Sample=="RIFEVO23","Position"]))
for (pos in RIFEVO1_rev_sites){
  info<-ANC_coding[ANC_coding$Position==pos,]
  info[1]<-"rev_ref"
  EEP_breseq_evo[nrow(EEP_breseq_evo)+1,1:7]<-info
  EEP_breseq_evo[nrow(EEP_breseq_evo),8]<-"RIFEVO23"
}
