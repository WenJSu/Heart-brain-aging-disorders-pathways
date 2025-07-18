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

diseaseCode <- read.table('/home1/suwenjing/20250320zhouzairen/mediation/data/Target_code.csv',sep=',',header=T)
diseaseCode <- diseaseCode[,c('Disease_code','Disease')]
colnames(diseaseCode) <- c('target','Disease')

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

#write.table(usedata, 'G:/SU_otherwork/20241103周在人/1122_brainage_heartIDPs/整理/5_mediate_analysis/use_data_illV2.csv',sep=',',row.names = F)
num_list <- c()
CFI <- c()
RMSEA <- c()
SRMR <- c()

list_disease <- c()
list_cate <- c()
list_heart <- c()
list_direction <- c()

list_a <- c()
list_b <- c()
list_c <- c()
list_ab <- c()
list_total <- c()

list_apvalues <- c()
list_bpvalues <- c()
list_cpvalues <- c()
list_abpvalues <- c()
list_totalpvalues <- c()


ill_list <- c('F0','F2','F4','F5','G0','G1','G2','G4')
heart_cate <- c('left atrium', 'left ventricle','right atrium','right ventricle')
heart_direction <- c('positive','negative')

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
        use_ill_heart_cate_dire <- use_ill_heart_cate[use_ill_heart_cate$diseases>0, ]
      }
      if(k==2){
        use_ill_heart_cate_dire <- use_ill_heart_cate[use_ill_heart_cate$diseases<0, ]
      }
      if(nrow(use_ill_heart_cate_dire) > 0){
        list_disease <- c(list_disease, ill_use)
        list_cate <- c(list_cate, cate_use)
        list_direction <- c(list_direction, direction_use)
        
        cols_to_extract <- use_ill_heart_cate_dire$heart_index
        cols_to_extract <- c("eid", use_ill_heart_cate_dire$heart_index)
        heart_cate_choose <- heartdata[, cols_to_extract]
        
        usedata_merge <- merge(heart_cate_choose, usedata, by=c('eid'))
        num_list <- c(num_list, nrow(usedata_merge))
        
        
        if(nrow(use_ill_heart_cate_dire) == 1){
          heart_name = use_ill_heart_cate_dire$heart_index 
          usedata_merge$mean <- usedata_merge[[heart_name]]
        }
        if(nrow(use_ill_heart_cate_dire) > 1){
          usedata_merge$mean <- rowMeans(usedata_merge[, use_ill_heart_cate_dire$heart_index], na.rm = TRUE)
        }
        usedata_merge$mean <- scale(usedata_merge$mean)
        list_heart <- c(list_heart, paste(use_ill_heart_cate_dire$heart_index, collapse = ","))
        
          model <- paste(
          "delta_corrected ~ a*mean+age2+Sex",
          paste0(ill_use, " ~ b*delta_corrected+age2+Sex"),
          paste0(ill_use, " ~ cp*mean+age2+Sex"),
          "ab := a * b",
          "total := cp + ab",sep='\n')
        #sem_result = sem(model, usedata_merge, se = "bootstrap", bootstrap = 1000)
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
        
        list_a <- c(list_a, result[["pe"]][["est"]][1])
        list_b <- c(list_b, result[["pe"]][["est"]][4])
        list_c <- c(list_c, result[["pe"]][["est"]][7])
        list_ab <- c(list_ab, result[["pe"]][["est"]][21])
        list_total <- c(list_total, result[["pe"]][["est"]][22])
        
        list_apvalues <- c(list_apvalues, result[["pe"]][["pvalue"]][1])
        list_bpvalues <- c(list_bpvalues, result[["pe"]][["pvalue"]][4])
        list_cpvalues <- c(list_cpvalues, result[["pe"]][["pvalue"]][7])
        list_abpvalues <- c(list_abpvalues, result[["pe"]][["pvalue"]][21])
        list_totalpvalues <- c(list_totalpvalues, result[["pe"]][["pvalue"]][22])
      }
    }
  }
}



Padj_FDR_a <- p.adjust(list_apvalues,method='BH')
Padj_FDR_b <- p.adjust(list_bpvalues,method='BH')
Padj_FDR_c <- p.adjust(list_cpvalues,method='BH')
Padj_FDR_ab <- p.adjust(list_abpvalues,method='BH')
Padj_FDR_total <- p.adjust(list_totalpvalues,method='BH')



resultframe <- data.frame(target=list_disease, heart_features=list_heart, heart_cate=list_cate, direction=list_direction,
                          num=num_list, CFI=CFI, RMSEA=RMSEA, SRMR=SRMR, a_estimate=list_a,
                          a_pvalues=list_apvalues, a_Padj_FDR=Padj_FDR_a, b_estimate=list_b,
                          b_pvalues=list_bpvalues, b_Padj_FDR=Padj_FDR_b, c_estimate=list_c,
                          c_pvalues=list_cpvalues, c_Padj_FDR=Padj_FDR_c,
                          ab_estimate=list_ab,
                          ab_pvalues=list_abpvalues, ab_Padj_FDR=Padj_FDR_ab, total_estimate=list_total,
                          total_pvalues=list_totalpvalues, total_Padj_FDR=Padj_FDR_total)

resultframe <- merge(resultframe, diseaseCode, by='target', all.x = T)
# heartnote <- read.table('G:/SU_otherwork/20241103周在人/1112heart分类/1_heart_name_cate.csv',sep=',',header=T)
# heartnote <- heartnote[,c('fid','Full_name')]
# colnames(heartnote) <- c('heart_features','Full_name')
# resultframe <- merge(resultframe, heartnote, by='heart_features', all.x = T)

write.table(resultframe, '/home1/suwenjing/20250320zhouzairen/mediation/result/heart_brain_diseases_cate_mediate_0320.csv',sep=',',row.names = F)
print("Finish!!!")
