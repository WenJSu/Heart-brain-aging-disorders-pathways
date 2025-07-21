library(lavaan)
library(semTools)
library(future)
library(future.apply)

heartdata = read.csv('Heart_fill_re.csv')
braindata_or = read.csv('agegap.csv')
braindata_or <- braindata_or[,c('eid','delta_corrected')]
covadata = read.csv('physical_merge.csv')
covadata <- covadata[,c('eid','age2','Sex','interval')]
usedata <- merge(braindata_or, covadata, by=c('eid'))
usedata$age2 <- scale(usedata$age2)
usedata$delta_corrected <- scale(usedata$delta_corrected)
diseaseCode <- read.table('Target_code.csv',sep=',',header=T)
diseaseCode <- diseaseCode[,c('Disease_code','Disease')]
colnames(diseaseCode) <- c('target','Disease')
common_df = read.csv('processV2.csv')

targetfiles <- list.files('./',pattern = '.csv')
for(j in seq(1,length(targetfiles))){
  targetfile <- targetfiles[j]
  if(targetfile=="F0.csv"|targetfile=="F1.csv"|targetfile=="F2.csv"|targetfile=="F4.csv"|targetfile=="F5.csv"|
     targetfile=="G0.csv"|targetfile=="G1.csv"|targetfile=="G2.csv"|targetfile=="G3.csv"|targetfile=="G4.csv"){
    targetdata <- read.csv(paste0('./',targetfile),sep=',',header=T)
    targetdata <- targetdata[,c('eid','target_y','BL2Target_yrs')]
    targetname <- gsub('.csv','',targetfile)
    colnames(targetdata) <- c('eid',targetname,'BL2Target_yrs')
    usedata <- merge(usedata, targetdata, by=c('eid'))
    
    usedata$time <- usedata$BL2Target_yrs-usedata$interval 
    usedata[[targetname]][usedata$time > 0] <- 0
    
    usedata <- usedata[, !(names(usedata) %in% c("BL2Target_yrs", "time"))]
  }
}

num_list <- c()
CFI <- c()
RMSEA <- c()
SRMR <- c()

list_disease <- c()
list_cate <- c()
list_heart <- c()
list_direction <- c()

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
          ordered = c(ill_use),
          estimator = "DWLS",
          se = "bootstrap",
          bootstrap = 1000
        )
        result <- summary(sem_result, fit.measures = TRUE, standardized=TRUE)
      }
    }
  }
}