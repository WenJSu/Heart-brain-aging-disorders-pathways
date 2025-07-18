library("fastDummies")
setwd('G:/SU_otherwork/20241103周在人/data/')

heartdata = read.csv('Heart_fill_re.csv')
heart_index <- colnames(heartdata)
covadata = read.csv('physical_merge.csv')


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

for(i in seq(2,length(heart_index))){

heart_single <- heartdata[c(1,i)]
pvalues <- c()
for(j in seq(1,length(targetfiles))){
  targetfile <- targetfiles[j]
  if(targetfile=="F0.csv"|targetfile=="F1.csv"|targetfile=="F2.csv"|targetfile=="F3.csv"|targetfile=="F4.csv"|targetfile=="F5.csv"|
     targetfile=="G0.csv"|targetfile=="G1.csv"|targetfile=="G2.csv"|targetfile=="G3.csv"|targetfile=="G4.csv"){
  list_agingfile <- c(list_agingfile, heart_index[i])
  list_targetfile <- c(list_targetfile, targetfile)
  # convert variables
  covadata_used <- covadata[,c('eid','age2','Sex','interval','TDI','smoking2','drinking2','BMI2')]
  covadata_used <- covadata_used[!is.na(covadata_used$smoking2) & !is.na(covadata_used$drinking2) & covadata_used$smoking2!=-3 & covadata_used$drinking2!=-3,]
  covadata_used <- fastDummies::dummy_cols(covadata_used, select_columns = "drinking2", remove_first_dummy = TRUE)
  covadata_used <- fastDummies::dummy_cols(covadata_used, select_columns = "smoking2", remove_first_dummy = TRUE)
  tempdata <- merge(heart_single, covadata_used, by=c('eid'))
  
  targetdata <- read.csv(paste0('E:/UKB/diseases/Targets_RAW/Targets_RAW/',targetfile),sep=',',header=T)
  useddata <- merge(tempdata, targetdata, by=c('eid'))
  
  useddata$time <- useddata$BL2Target_yrs-useddata$interval 
  useddata$target_y[useddata$time > 0] <- 0
  
  list_numberAll <- c(list_numberAll, nrow(useddata))
  list_numberTarget <- c(list_numberTarget, nrow(useddata[useddata$target_y==1,]))
  
  # cross-section analysis
  cleaned_data <- useddata[,c('age2',heart_index[i],'Sex','target_y','smoking2_1','smoking2_2','drinking2_1','drinking2_2','BMI2','TDI')]
  colnames(cleaned_data) <- c('age2','heart','Sex','status', 'smoking2_1','smoking2_2','drinking2_1','drinking2_2','BMI2','TDI')
  
  model <- glm(status ~ heart + age2 + Sex+smoking2_1+smoking2_2+drinking2_1+drinking2_2+BMI2+TDI, data = cleaned_data, family=binomial)
  summary_model <- summary(model)
  coefficients <- summary_model$coefficients
  
  list_coef <- c(list_coef, coefficients[2,1])
  list_std <- c(list_std, coefficients[2,2])
  list_tvalues <- c(list_tvalues, coefficients[2,3])
  list_pvalues <- c(list_pvalues, coefficients[2,4])
  
  pvalues <- c(pvalues, coefficients[2,4])
  }
}
print(i)
list_FDR_separate <- c(list_FDR_separate,p.adjust(pvalues, method = 'BH'))
}

resultframe <- data.frame(heart_index=list_agingfile, target=list_targetfile, numberAll=list_numberAll, 
                          numberTarget=list_numberTarget, coef=list_coef, std=list_std,
                          tvalue=list_tvalues, 
                          pvalue=list_pvalues, separateFDR=list_FDR_separate)
resultframe$target <- sapply(resultframe$target, function(x) gsub('.csv','',x))
resultframe <- merge(resultframe, diseaseCode, by='target', all.x = T)
resultframe$Padj_FDR_overall <- p.adjust(resultframe$pvalue,method='BH')
write.table(resultframe, 'G:/SU_otherwork/20241103周在人/1122_brainage_heartIDPs/整理/4_glm_heart_ill/glm_heat_ill_prevalent_0317_V2.csv',sep=',',row.names = F)
print("Finish!")