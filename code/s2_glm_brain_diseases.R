library("fastDummies")

braindata = read.csv('G:/SU_otherwork/20241103周在人/1122_brainage_heartIDPs/preidict_brain_age/result/agegap.csv')
braindata <- braindata[,c('eid','delta_corrected')]

covadata = read.csv('G:/SU_otherwork/20241103周在人/1122_brainage_heartIDPs/0_datapreprocess/physical_merge.csv')


targetfiles <- list.files('E:/UKB/diseases/Targets_RAW/Targets_RAW/',pattern = '.csv')
diseaseCode <- read.table('E:/UKB/diseases/Target_code.csv',sep=',',header=T)
diseaseCode <- diseaseCode[,c('Disease_code','Disease')]
colnames(diseaseCode) <- c('target','Disease')

list_agingfile <- c()
list_targetfile <- c()
list_numberAll <- c()
list_numberTarget <- c()
list_coef <- c()
list_std <- c()
list_tvalues <- c()
list_pvalues <- c()
list_FDR_separate <- c()


pvalues <- c()
for(j in seq(1,length(targetfiles))){
  targetfile <- targetfiles[j]
  
  list_agingfile <- c(list_agingfile, "Brain_GM_age_gap")
  list_targetfile <- c(list_targetfile, targetfile)
  # convert variables
  covadata_used <- covadata[,c('eid','age2','Sex','interval','TDI','smoking2','drinking2','BMI2')]
  covadata_used <- covadata_used[!is.na(covadata_used$smoking2) & !is.na(covadata_used$drinking2) & covadata_used$smoking2!=-3 & covadata_used$drinking2!=-3,]
  covadata_used <- fastDummies::dummy_cols(covadata_used, select_columns = "drinking2", remove_first_dummy = TRUE)
  covadata_used <- fastDummies::dummy_cols(covadata_used, select_columns = "smoking2", remove_first_dummy = TRUE)
  tempdata <- merge(braindata, covadata_used, by=c('eid'))
  
  targetdata <- read.csv(paste0('E:/UKB/diseases/Targets_RAW/Targets_RAW/',targetfile),sep=',',header=T)
  useddata <- merge(tempdata, targetdata, by=c('eid'))
  
  useddata$time <- useddata$BL2Target_yrs-useddata$interval 
  useddata$target_y[useddata$time > 0] <- 0
  
  list_numberAll <- c(list_numberAll, nrow(useddata))
  list_numberTarget <- c(list_numberTarget, nrow(useddata[useddata$target_y==1,]))
  
  # cross-section analysis
  cleaned_data <- useddata[,c('age2','delta_corrected','Sex','target_y','smoking2_1','smoking2_2','drinking2_1','drinking2_2','BMI2','TDI')]
  colnames(cleaned_data) <- c('age2','brain','Sex','status', 'smoking2_1','smoking2_2','drinking2_1','drinking2_2','BMI2','TDI')
  
  model <- glm(status ~ brain + age2 + Sex+smoking2_1+smoking2_2+drinking2_1+drinking2_2+BMI2+TDI, data = cleaned_data)
  summary_model <- summary(model)
  coefficients <- summary_model$coefficients
  
  list_coef <- c(list_coef, coefficients[2,1])
  list_std <- c(list_std, coefficients[2,2])
  list_tvalues <- c(list_tvalues, coefficients[2,3])
  list_pvalues <- c(list_pvalues, coefficients[2,4])
  
  pvalues <- c(pvalues, coefficients[2,4])
  
  if(j%%10==0){
    print(j)
  }
  
}
list_FDR_separate <- c(list_FDR_separate,p.adjust(pvalues, method = 'BH'))


resultframe <- data.frame(braint_index=list_agingfile, target=list_targetfile, numberAll=list_numberAll, 
                          numberTarget=list_numberTarget, coef=list_coef, std=list_std,
                          tvalue=list_tvalues, 
                          pvalue=list_pvalues, separateFDR=list_FDR_separate)
resultframe$target <- sapply(resultframe$target, function(x) gsub('.csv','',x))
resultframe <- merge(resultframe, diseaseCode, by='target', all.x = T)
resultframe$Padj_FDR_overall <- p.adjust(resultframe$pvalue,method='BH')
resultframe$Padj_Bon_overall <- p.adjust(resultframe$pvalue,method='bonferroni')
write.table(resultframe, 'G:/SU_otherwork/20241103周在人/1122_brainage_heartIDPs/preidict_brain_age/2_glm_brain_ill/glm_brain_ill.csv',sep=',',row.names = F)
print("Finish!")