library(lavaan)
library(semTools)

library(future)
library(future.apply)
plan(multisession, workers=parallel::detectCores()-100)

heartdata = read.csv('/home1/suwenjing/20250320zhouzairen/mediation/data/Heart_fill_re.csv')

braindata_or = read.csv('/home1/suwenjing/20250320zhouzairen/mediation/data/agegap.csv')
braindata_or <- braindata_or[,c('eid','delta_corrected')]

covadata = read.csv('/home1/suwenjing/20250320zhouzairen/mediation/data/physical_merge.csv')
covadata <- covadata[,c('eid','age2','Sex','interval')]

usedata <- merge(braindata_or, covadata, by=c('eid'))
usedata$age2 <- scale(usedata$age2)
usedata$delta_corrected <- scale(usedata$delta_corrected)


#ill_heart = read.csv('G:/SU_otherwork/20241103周在人/1122_brainage_heartIDPs/4_SEM/ill_relevant.csv')
#brain_heart = read.csv('G:/SU_otherwork/20241103周在人/1122_brainage_heartIDPs/4_SEM/brain_relevant.csv')
common_df = read.csv('/home1/suwenjing/20250320zhouzairen/mediation/data/processV2.csv')


targetfiles <- list.files('/home1/suwenjing/20250320zhouzairen/mediation/data/Targets_RAW/Targets_RAW/',pattern = '.csv')
for(j in seq(1,length(targetfiles))){
  targetfile <- targetfiles[j]
  if(targetfile=="F0.csv"|targetfile=="F1.csv"|targetfile=="F2.csv"|targetfile=="F4.csv"|targetfile=="F5.csv"|
     targetfile=="G0.csv"|targetfile=="G1.csv"|targetfile=="G2.csv"|targetfile=="G3.csv"|targetfile=="G4.csv"){
    targetdata <- read.csv(paste0('/home1/suwenjing/20250320zhouzairen/mediation/data/Targets_RAW/Targets_RAW/',targetfile),sep=',',header=T)
    targetdata <- targetdata[,c('eid','target_y','BL2Target_yrs')]
    targetname <- gsub('.csv','',targetfile)
    colnames(targetdata) <- c('eid',targetname,'BL2Target_yrs')
    usedata <- merge(usedata, targetdata, by=c('eid'))
    
    usedata$time <- usedata$BL2Target_yrs-usedata$interval 
    usedata[[targetname]][usedata$time > 0] <- 0
    
    usedata <- usedata[, !(names(usedata) %in% c("BL2Target_yrs", "time"))]
  }
  if(j%%10==0){
    print(j)
  }
}


CFI <- c()
RMSEA <- c()
SRMR <- c()

list_disease <- c()
list_cate <- c()
list_heart <- c()
list_direction <- c()
num_list <- c()

list_a <- c()
list_b <- c()
list_ab <- c()
list_cp <- c()
list_d <- c()
list_ep <- c()
list_abd <- c()
list_total <- c()
list_dtotal <- c()
list_all <- c()

list_ap <- c()
list_bp <- c()
list_abp <- c()
list_cpp <- c()
list_dp <- c()
list_epp <- c()
list_abdp <- c()
list_totalp <- c()
list_dtotalp <- c()
list_allp <- c()

FDR_a <- c()
FDR_b <- c()
FDR_ab <- c()
FDR_cp <- c()
FDR_d <- c()
FDR_ep <- c()
FDR_abd <- c()
FDR_total <- c()
FDR_dtotal <- c()
FDR_all <- c()

fID <- c()

lianxu_factor <- c('X22033.2.0', 'X22189.0.0', 'X4079', 'X4080', 'X48.2.0','X1269.2.0','X1279.2.0')
ill_list <- c('F0','F2','F4','F5','G0')
heart_cate <- c('left atrium', 'left ventricle','right atrium','right ventricle')
heart_direction <- c('positive','negative')

for(facid in seq(1,length(lianxu_factor))){ 
  use_fac <- lianxu_factor[facid]
  print(use_fac)
  m_factor = read.csv('/home1/suwenjing/20250320zhouzairen/mediation/data/modifiable_factors_meanbloodpressure.csv')
  m_factor_use <- m_factor[,c('eid',use_fac)]
  colnames(m_factor_use) <- c('eid','lifestyle')
  
  if (use_fac == 'X22189.0.0') {
    m_factor_use <- m_factor_use[!is.na(m_factor_use$lifestyle),]
  } else {
    m_factor_use <- m_factor_use[!is.na(m_factor_use$lifestyle) & m_factor_use$lifestyle!=-3 & m_factor_use$lifestyle!=-1,]
  }

  usedata_singlefac <- merge(usedata, m_factor_use, by=c('eid'))
  usedata_singlefac$lifestyle <- scale(usedata_singlefac$lifestyle)
  
  list_ap1 <- c()
  list_bp1 <- c()
  list_abp1 <- c()
  list_cpp1 <- c()
  list_dp1 <- c()
  list_epp1 <- c()
  list_abdp1 <- c()
  list_totalp1 <- c()
  list_dtotalp1 <- c()
  list_allp1 <- c()
  
  for(i in seq(1,length(ill_list))){
    print(i)
    ill_use = ill_list[i]
    use_ill_heart <- common_df[common_df$target == ill_use, ]
    for(t in seq(1,length(heart_cate))){
      cate_use = heart_cate[t]
      print(cate_use)
      use_ill_heart_cate <- use_ill_heart[use_ill_heart$Category == cate_use, ]
      for(k in seq(1,length(heart_direction))){
        direction_use = heart_direction[k]
        print(direction_use)
        if(k==1){
          use_ill_heart_cate_dire  <- use_ill_heart_cate[use_ill_heart_cate$diseases>0, ]
        }
        if(k==2){
          use_ill_heart_cate_dire  <- use_ill_heart_cate[use_ill_heart_cate$diseases<0, ]
        }
        if(nrow(use_ill_heart_cate_dire) > 0){
          list_disease <- c(list_disease, ill_use)
          list_cate <- c(list_cate, cate_use)
          fID <- c(fID, use_fac)
          list_direction <- c(list_direction, direction_use)
          
          cols_to_extract <- use_ill_heart_cate_dire $heart_index
          cols_to_extract <- c("eid", use_ill_heart_cate_dire $heart_index)
          heart_cate_choose <- heartdata[, cols_to_extract]
          
          usedata_merge <- merge(heart_cate_choose, usedata_singlefac, by=c('eid'))
          num_list <- c(num_list, nrow(usedata_merge))
          
          
          if(nrow(use_ill_heart_cate_dire ) == 1){
            heart_name = use_ill_heart_cate_dire $heart_index 
            usedata_merge$mean <- usedata_merge[[heart_name]]
          }
          if(nrow(use_ill_heart_cate_dire ) > 1){
            usedata_merge$mean <- rowMeans(usedata_merge[, use_ill_heart_cate_dire $heart_index], na.rm = TRUE)
          }
          list_heart <- c(list_heart, paste(use_ill_heart_cate_dire$heart_index, collapse = ","))
          usedata_merge$mean <- scale(usedata_merge$mean)
          
          model <- paste(
            "mean ~ d*lifestyle+age2+Sex",
            "delta_corrected ~ a*mean+age2+Sex",
            paste0(ill_use, " ~ b*delta_corrected+age2+Sex"),
            paste0(ill_use, " ~ cp*mean+age2+Sex"),
            paste0(ill_use, " ~ ep*lifestyle+age2+Sex"),
            "ab := a * b",
            "abd := a*b*d",
            "total := cp + ab",
            "dtotal := d*total",
            "all:=ep+dtotal",sep='\n'
          )

          sem_result <- sem(
            model,
            data = usedata_merge,
            ordered = c(ill_use),    # 确保变量名用引号包裹，且存在于数据中
            estimator = "DWLS",       # 或 estimator = "DWLS"
            se = "bootstrap",
            bootstrap = 1000
          )
          result <- summary(sem_result, fit.measures = TRUE, standardized=TRUE)
          
          CFI <- c(CFI, result[["fit"]][["cfi"]])
          RMSEA <- c(RMSEA, result[["fit"]][["rmsea"]])
          SRMR <- c(SRMR, result[["fit"]][["srmr"]])
          
          list_a <- c(list_a, result[["pe"]][["est"]][4])
          list_b <- c(list_b, result[["pe"]][["est"]][7])
          list_ab <- c(list_ab, result[["pe"]][["est"]][27])
          list_cp <- c(list_cp, result[["pe"]][["est"]][10])
          list_d <- c(list_d, result[["pe"]][["est"]][1])
          list_ep <- c(list_ep, result[["pe"]][["est"]][11])
          list_abd <- c(list_abd, result[["pe"]][["est"]][28])
          list_total <- c(list_total, result[["pe"]][["est"]][29])
          list_dtotal <- c(list_dtotal, result[["pe"]][["est"]][30])
          list_all <- c(list_all, result[["pe"]][["est"]][31])
          
          list_ap <- c(list_ap, result[["pe"]][["pvalue"]][4])
          list_bp <- c(list_bp, result[["pe"]][["pvalue"]][7])
          list_abp <- c(list_abp, result[["pe"]][["pvalue"]][27])
          list_cpp <- c(list_cpp, result[["pe"]][["pvalue"]][10])
          list_dp <- c(list_dp, result[["pe"]][["pvalue"]][1])
          list_epp <- c(list_epp, result[["pe"]][["pvalue"]][11])
          list_abdp <- c(list_abdp, result[["pe"]][["pvalue"]][28])
          list_totalp <- c(list_totalp, result[["pe"]][["pvalue"]][29])
          list_dtotalp <- c(list_dtotalp, result[["pe"]][["pvalue"]][30])
          list_allp <- c(list_allp, result[["pe"]][["pvalue"]][31])
          
          list_ap1 <- c(list_ap1, result[["pe"]][["pvalue"]][4])
          list_bp1 <- c(list_bp1, result[["pe"]][["pvalue"]][7])
          list_abp1 <- c(list_abp1, result[["pe"]][["pvalue"]][27])
          list_cpp1 <- c(list_cpp1, result[["pe"]][["pvalue"]][10])
          list_dp1 <- c(list_dp1, result[["pe"]][["pvalue"]][1])
          list_epp1 <- c(list_epp1, result[["pe"]][["pvalue"]][11])
          list_abdp1 <- c(list_abdp1, result[["pe"]][["pvalue"]][28])
          list_totalp1 <- c(list_totalp1, result[["pe"]][["pvalue"]][29])
          list_dtotalp1 <- c(list_dtotalp1, result[["pe"]][["pvalue"]][30])
          list_allp1 <- c(list_allp1, result[["pe"]][["pvalue"]][31])
        }
      }
    }
  }
  FDR_a <- c(FDR_a, p.adjust(list_ap1,method='BH'))
  FDR_b <- c(FDR_b, p.adjust(list_bp1,method='BH'))
  FDR_ab <- c(FDR_ab, p.adjust(list_abp1,method='BH'))
  FDR_cp <- c(FDR_cp, p.adjust(list_cpp1,method='BH'))
  FDR_d <- c(FDR_d, p.adjust(list_dp1,method='BH'))
  FDR_ep <- c(FDR_ep, p.adjust(list_epp1,method='BH'))
  FDR_abd <- c(FDR_abd, p.adjust(list_abdp1,method='BH'))
  FDR_total <- c(FDR_total, p.adjust(list_totalp1,method='BH'))
  FDR_dtotal <- c(FDR_dtotal, p.adjust(list_dtotalp1,method='BH'))
  FDR_all <- c(FDR_all, p.adjust(list_allp1,method='BH'))
}
    
    

resultframe <- data.frame(diseases=list_disease, heart_features=list_heart, heart_cate=list_cate, direction=list_direction,
                          num=num_list, CFI=CFI, RMSEA=RMSEA, SRMR=SRMR, a=list_a, a_pvalues=list_ap, FDR_a=FDR_a, 
                          b=list_b, b_pvalues=list_bp, FDR_b=FDR_b, ab=list_ab, ab_pvalues=list_abp, FDR_ab=FDR_ab,
                          cp=list_cp, cp_pvalues=list_cpp, FDR_cp=FDR_cp, d=list_d, d_pvalues=list_dp, FDR_d=FDR_d,
                          ep=list_ep, ep_pvalues=list_epp, FDR_ep=FDR_ep, abd=list_abd, abd_pvalues=list_abdp, FDR_abd=FDR_abd,
                          total=list_total, total_pvalues=list_totalp, FDR_total=FDR_total,
                          dtotal=list_dtotal, dtotal_pvalues=list_dtotalp, FDR_dtotal=FDR_dtotal,
                          all=list_all, all_pvalues=list_allp,FDR_all=FDR_all, fID=fID)

fID_note <- read.table('/home1/suwenjing/20250320zhouzairen/mediation/data/lianxu.csv',sep=',',header=T)
resultframe <- merge(resultframe, fID_note, by='fID', all.x = T)

diseaseCode <- read.table('/home1/suwenjing/20250320zhouzairen/mediation/data/Target_code.csv',sep=',',header=T)
diseaseCode <- diseaseCode[,c('Disease_code','Disease')]
colnames(diseaseCode) <- c('diseases','Disease_name')
resultframe <- merge(resultframe, diseaseCode, by='diseases', all.x = T)

# heartnote <- read.table('G:/SU_otherwork/20241103周在人/1112heart分类/1_heart_name_cate.csv',sep=',',header=T)
# heartnote <- heartnote[,c('fid','Full_name')]
# colnames(heartnote) <- c('heart_features','Heart_full_name')
# resultframe <- merge(resultframe, heartnote, by='heart_features', all.x = T)

write.table(resultframe, '/home1/suwenjing/20250320zhouzairen/factors/result/1_lianxuV1_V2.csv',sep=',',row.names = F)

#######################################################################
use_fac <- "Education"
print(use_fac)
m_factor = read.csv('/home1/suwenjing/UKB/cov/cova_data_more_withPC_all2visits.csv')
m_factor_use <- m_factor[,c('eid', 'Education')]
colnames(m_factor_use) <- c('eid','lifestyle')
m_factor_use <- m_factor_use[!is.na(m_factor_use$lifestyle),]
usedata_singlefac <- merge(usedata, m_factor_use, by=c('eid'))
usedata_singlefac$lifestyle <- scale(usedata_singlefac$lifestyle)
list_ap1 <- c()
list_bp1 <- c()
list_abp1 <- c()
list_cpp1 <- c()
list_dp1 <- c()
list_epp1 <- c()
list_abdp1 <- c()
list_totalp1 <- c()
list_dtotalp1 <- c()
list_allp1 <- c()

for(i in seq(1,length(ill_list))){
  print(i)
  ill_use = ill_list[i]
  use_ill_heart <- common_df[common_df$target == ill_use, ]
  for(t in seq(1,length(heart_cate))){
    cate_use = heart_cate[t]
    print(cate_use)
    use_ill_heart_cate <- use_ill_heart[use_ill_heart$Category == cate_use, ]
    for(k in seq(1,length(heart_direction))){
      direction_use = heart_direction[k]
      print(direction_use)
      if(k==1){
        use_ill_heart_cate_dire  <- use_ill_heart_cate[use_ill_heart_cate$diseases>0, ]
      }
      if(k==2){
        use_ill_heart_cate_dire  <- use_ill_heart_cate[use_ill_heart_cate$diseases<0, ]
      }
      if(nrow(use_ill_heart_cate_dire) > 0){
        list_disease <- c(list_disease, ill_use)
        list_cate <- c(list_cate, cate_use)
        fID <- c(fID, use_fac)
        list_direction <- c(list_direction, direction_use)
        
        cols_to_extract <- use_ill_heart_cate_dire $heart_index
        cols_to_extract <- c("eid", use_ill_heart_cate_dire $heart_index)
        heart_cate_choose <- heartdata[, cols_to_extract]
        
        usedata_merge <- merge(heart_cate_choose, usedata_singlefac, by=c('eid'))
        num_list <- c(num_list, nrow(usedata_merge))
        
        
        if(nrow(use_ill_heart_cate_dire ) == 1){
          heart_name = use_ill_heart_cate_dire $heart_index 
          usedata_merge$mean <- usedata_merge[[heart_name]]
        }
        if(nrow(use_ill_heart_cate_dire ) > 1){
          usedata_merge$mean <- rowMeans(usedata_merge[, use_ill_heart_cate_dire $heart_index], na.rm = TRUE)
        }
        list_heart <- c(list_heart, paste(use_ill_heart_cate_dire$heart_index, collapse = ","))
        usedata_merge$mean <- scale(usedata_merge$mean)
        
        model <- paste(
          "mean ~ d*lifestyle+age2+Sex",
          "delta_corrected ~ a*mean+age2+Sex",
          paste0(ill_use, " ~ b*delta_corrected+age2+Sex"),
          paste0(ill_use, " ~ cp*mean+age2+Sex"),
          paste0(ill_use, " ~ ep*lifestyle+age2+Sex"),
          "ab := a * b",
          "abd := a*b*d",
          "total := cp + ab",
          "dtotal := d*total",
          "all:=ep+dtotal",sep='\n'
        )
        
        sem_result <- sem(
          model,
          data = usedata_merge,
          ordered = c(ill_use),    # 确保变量名用引号包裹，且存在于数据中
          estimator = "DWLS",       # 或 estimator = "DWLS"
          se = "bootstrap",
          bootstrap = 1000
        )
        result <- summary(sem_result, fit.measures = TRUE, standardized=TRUE)
        
        CFI <- c(CFI, result[["fit"]][["cfi"]])
        RMSEA <- c(RMSEA, result[["fit"]][["rmsea"]])
        SRMR <- c(SRMR, result[["fit"]][["srmr"]])
        
        list_a <- c(list_a, result[["pe"]][["est"]][4])
        list_b <- c(list_b, result[["pe"]][["est"]][7])
        list_ab <- c(list_ab, result[["pe"]][["est"]][27])
        list_cp <- c(list_cp, result[["pe"]][["est"]][10])
        list_d <- c(list_d, result[["pe"]][["est"]][1])
        list_ep <- c(list_ep, result[["pe"]][["est"]][11])
        list_abd <- c(list_abd, result[["pe"]][["est"]][28])
        list_total <- c(list_total, result[["pe"]][["est"]][29])
        list_dtotal <- c(list_dtotal, result[["pe"]][["est"]][30])
        list_all <- c(list_all, result[["pe"]][["est"]][31])
        
        list_ap <- c(list_ap, result[["pe"]][["pvalue"]][4])
        list_bp <- c(list_bp, result[["pe"]][["pvalue"]][7])
        list_abp <- c(list_abp, result[["pe"]][["pvalue"]][27])
        list_cpp <- c(list_cpp, result[["pe"]][["pvalue"]][10])
        list_dp <- c(list_dp, result[["pe"]][["pvalue"]][1])
        list_epp <- c(list_epp, result[["pe"]][["pvalue"]][11])
        list_abdp <- c(list_abdp, result[["pe"]][["pvalue"]][28])
        list_totalp <- c(list_totalp, result[["pe"]][["pvalue"]][29])
        list_dtotalp <- c(list_dtotalp, result[["pe"]][["pvalue"]][30])
        list_allp <- c(list_allp, result[["pe"]][["pvalue"]][31])
        
        list_ap1 <- c(list_ap1, result[["pe"]][["pvalue"]][4])
        list_bp1 <- c(list_bp1, result[["pe"]][["pvalue"]][7])
        list_abp1 <- c(list_abp1, result[["pe"]][["pvalue"]][27])
        list_cpp1 <- c(list_cpp1, result[["pe"]][["pvalue"]][10])
        list_dp1 <- c(list_dp1, result[["pe"]][["pvalue"]][1])
        list_epp1 <- c(list_epp1, result[["pe"]][["pvalue"]][11])
        list_abdp1 <- c(list_abdp1, result[["pe"]][["pvalue"]][28])
        list_totalp1 <- c(list_totalp1, result[["pe"]][["pvalue"]][29])
        list_dtotalp1 <- c(list_dtotalp1, result[["pe"]][["pvalue"]][30])
        list_allp1 <- c(list_allp1, result[["pe"]][["pvalue"]][31])
      }
    }
  }
}
FDR_a <- c(FDR_a, p.adjust(list_ap1,method='BH'))
FDR_b <- c(FDR_b, p.adjust(list_bp1,method='BH'))
FDR_ab <- c(FDR_ab, p.adjust(list_abp1,method='BH'))
FDR_cp <- c(FDR_cp, p.adjust(list_cpp1,method='BH'))
FDR_d <- c(FDR_d, p.adjust(list_dp1,method='BH'))
FDR_ep <- c(FDR_ep, p.adjust(list_epp1,method='BH'))
FDR_abd <- c(FDR_abd, p.adjust(list_abdp1,method='BH'))
FDR_total <- c(FDR_total, p.adjust(list_totalp1,method='BH'))
FDR_dtotal <- c(FDR_dtotal, p.adjust(list_dtotalp1,method='BH'))
FDR_all <- c(FDR_all, p.adjust(list_allp1,method='BH'))

resultframe <- data.frame(diseases=list_disease, heart_features=list_heart, heart_cate=list_cate, direction=list_direction,
                          num=num_list, CFI=CFI, RMSEA=RMSEA, SRMR=SRMR, a=list_a, a_pvalues=list_ap, FDR_a=FDR_a, 
                          b=list_b, b_pvalues=list_bp, FDR_b=FDR_b, ab=list_ab, ab_pvalues=list_abp, FDR_ab=FDR_ab,
                          cp=list_cp, cp_pvalues=list_cpp, FDR_cp=FDR_cp, d=list_d, d_pvalues=list_dp, FDR_d=FDR_d,
                          ep=list_ep, ep_pvalues=list_epp, FDR_ep=FDR_ep, abd=list_abd, abd_pvalues=list_abdp, FDR_abd=FDR_abd,
                          total=list_total, total_pvalues=list_totalp, FDR_total=FDR_total,
                          dtotal=list_dtotal, dtotal_pvalues=list_dtotalp, FDR_dtotal=FDR_dtotal,
                          all=list_all, all_pvalues=list_allp,FDR_all=FDR_all, fID=fID)

fID_note <- read.table('/home1/suwenjing/20250320zhouzairen/mediation/data/lianxu.csv',sep=',',header=T)
resultframe <- merge(resultframe, fID_note, by='fID', all.x = T)

diseaseCode <- read.table('/home1/suwenjing/20250320zhouzairen/mediation/data/Target_code.csv',sep=',',header=T)
diseaseCode <- diseaseCode[,c('Disease_code','Disease')]
colnames(diseaseCode) <- c('diseases','Disease_name')
resultframe <- merge(resultframe, diseaseCode, by='diseases', all.x = T)

# heartnote <- read.table('G:/SU_otherwork/20241103周在人/1112heart分类/1_heart_name_cate.csv',sep=',',header=T)
# heartnote <- heartnote[,c('fid','Full_name')]
# colnames(heartnote) <- c('heart_features','Heart_full_name')
# resultframe <- merge(resultframe, heartnote, by='heart_features', all.x = T)

write.table(resultframe, '/home1/suwenjing/20250320zhouzairen/factors/result/1_lianxuV2_V2.csv',sep=',',row.names = F)

#######################################################################

