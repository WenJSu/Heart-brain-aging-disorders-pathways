# survial analysis
################################## with more covariates
library(survival)
library("fastDummies")
library(data.table)

heartdata = read.csv('C:/Users/admin/Desktop/zhouzairen/V2_height/Heart_fill_re.csv')
heart_index <- colnames(heartdata)

targetfiles <- list.files('C:/Users/admin/Desktop/BA/data/20240312diseases/Targets_RAW/Targets_RAW/',pattern = '.csv')
covafile <- read.table('C:/Users/admin/Desktop/zhouzairen/V2_height/physical_merge.csv',sep=',',header=T)
covadata <- covafile[,c('eid','age2','Sex','interval','TDI','smoking2','drinking2','BMI2','Education','eTIV','site')]

diseaseCode <- read.table('C:/Users/admin/Desktop/image_BA/0908renew/code/生存曲线/result/Target_code.csv',sep=',',header=T)
diseaseCode <- diseaseCode[,c('Disease_code','Disease')]
colnames(diseaseCode) <- c('target','Disease')

# cox harzard model, with delta as continuous variable, in all subjs, with more covariates
list_agingfile <- c()
list_targetfile <- c()
list_numberAll <- c()
list_numberTarget <- c()
list_pvalues <- c()
list_zvalues <- c()
list_coef <- c()
list_expcoef <- c()
list_secoef <- c()
list_FDR_separate <- c()
list_up <- c()
list_down <- c()
for(i in seq(2,length(heart_index))){

  heart_single <- heartdata[c(1,i)]
  pvalues <- c()
  for(j in seq(1,length(targetfiles))){
    targetfile <- targetfiles[j]
    list_agingfile <- c(list_agingfile, heart_index[i])
    list_targetfile <- c(list_targetfile, targetfile)
    # convert variables
    covadata_used <- covadata[,c('eid','age2','Sex','interval','TDI','smoking2','drinking2','BMI2')]
    covadata_used <- covadata_used[!is.na(covadata_used$smoking2) & !is.na(covadata_used$drinking2) & covadata_used$smoking2!=-3 & covadata_used$drinking2!=-3,]
    covadata_used <- fastDummies::dummy_cols(covadata_used, select_columns = "drinking2", remove_first_dummy = TRUE)
    covadata_used <- fastDummies::dummy_cols(covadata_used, select_columns = "smoking2", remove_first_dummy = TRUE)
    tempdata <- merge(heart_single, covadata_used, by=c('eid'))

    targetdata <- read.csv(paste0('C:/Users/admin/Desktop/BA/data/20240312diseases/Targets_RAW/Targets_RAW/',targetfile),sep=',',header=T)
    useddata <- merge(tempdata, targetdata, by=c('eid'))

    useddata$time <- useddata$BL2Target_yrs-useddata$interval 

    useddata <- useddata[useddata$time>0,]
    list_numberAll <- c(list_numberAll, nrow(useddata))
    list_numberTarget <- c(list_numberTarget, nrow(useddata[useddata$target_y==1,]))
    # survival analysis

    cleaned_data <- useddata[,c('age2',heart_index[i],'Sex','target_y','time','smoking2_1','smoking2_2','drinking2_1','drinking2_2','BMI2','TDI')]
    colnames(cleaned_data) <- c('age2','heart','Sex','status','years','smoking2_1','smoking2_2','drinking2_1','drinking2_2','BMI2','TDI')
    
    cox_fit <- coxph(Surv(years, status) ~ heart+Sex+age2+smoking2_1+smoking2_2+drinking2_1+drinking2_2+BMI2+TDI, data=cleaned_data)

    result <- summary(cox_fit)$coefficients
    CI_range <- summary(cox_fit)$conf.int[1,c(3,4)]
    # Check for the possibility of infinite coefficients
    if (any(is.infinite(coefficients(cox_fit)))) {
      print(j)
      print("Warning: Model coefficients may be infinite.")
    }
    list_pvalues <- c(list_pvalues, result[1,5])
    list_zvalues <- c(list_zvalues, result[1,4])
    list_coef <- c(list_coef, result[1,1])
    list_expcoef <- c(list_expcoef, result[1,2])
    list_secoef <- c(list_secoef, result[1,3])
    list_up <- c(list_up, CI_range[2])
    list_down <- c(list_down, CI_range[1])
    pvalues <- c(pvalues, result[1,5])
  }
  print(i)
  list_FDR_separate <- c(list_FDR_separate,p.adjust(pvalues, method = 'BH'))
}

resultframe <- data.frame(aging=list_agingfile, target=list_targetfile, numberAll=list_numberAll, 
                          numberTarget=list_numberTarget, expcoef=list_expcoef, coef=list_coef, 
                          down=unname(list_down), up=unname(list_up), secoef=list_secoef, zvalue=list_zvalues, 
                          pvalue=list_pvalues, separateFDR=list_FDR_separate)
resultframe$target <- sapply(resultframe$target, function(x) gsub('.csv','',x))
resultframe <- merge(resultframe, diseaseCode, by='target', all.x = T)
resultframe$Padj_FDR_overall <- p.adjust(resultframe$pvalue,method='BH')
write.table(resultframe, 'C:/Users/admin/Desktop/zhouzairen/V2_height/cox_heat_ill_V2.csv',sep=',',row.names = F)
print("Finish!")