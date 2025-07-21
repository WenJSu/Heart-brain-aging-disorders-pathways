library("fastDummies")
heartdata = read.csv('Heart_fill_re.csv')
heart_index <- colnames(heartdata)
covadata = read.csv('physical_merge.csv')
targetfiles <- list.files('Targets_RAW/',pattern = '.csv')
diseaseCode <- read.table('Target_code.csv',sep=',',header=T)
diseaseCode <- diseaseCode[,c('Disease_code','Disease')]
colnames(diseaseCode) <- c('target','Disease')


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
    
    targetdata <- read.csv(paste0(targetfile),sep=',',header=T)
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
    }
  }
}