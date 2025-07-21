# survial analysis
################################## with more covariates
library(survival)
library("fastDummies")
library(data.table)

heartdata = read.csv('Heart_fill_re.csv')
heart_index <- colnames(heartdata)
targetfiles <- list.files('Targets_RAW/',pattern = '.csv')
covafile <- read.table('physical_merge.csv',sep=',',header=T)
covadata <- covafile[,c('eid','age2','Sex','interval','TDI','smoking2','drinking2','BMI2','Education','eTIV','site')]
diseaseCode <- read.table('Target_code.csv',sep=',',header=T)
diseaseCode <- diseaseCode[,c('Disease_code','Disease')]
colnames(diseaseCode) <- c('target','Disease')

# cox harzard model, with delta as continuous variable, in all subjs, with more covariates
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
    targetdata <- read.csv(paste0('Targets_RAW/',targetfile),sep=',',header=T)
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
  }
}