library(lavaan)
library(semTools)

library(future)
library(future.apply)
plan(multisession, workers=parallel::detectCores()-100)

heartdata = read.csv('Heart_fill_re.csv')

braindata_or = read.csv('agegap.csv')
braindata_or <- braindata_or[,c('eid','delta_corrected')]

covadata = read.csv('physical_merge.csv')
covadata <- covadata[,c('eid','age2','Sex','interval')]

usedata <- merge(braindata_or, covadata, by=c('eid'))
usedata$age2 <- scale(usedata$age2)
usedata$delta_corrected <- scale(usedata$delta_corrected)
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

fID <- c()

lianxu_factor <- c('X22033.2.0', 'X22189.0.0', 'X4079', 'X4080', 'X48.2.0','X1269.2.0','X1279.2.0')
ill_list <- c('F0','F2','F4','F5','G0')
heart_cate <- c('left atrium', 'left ventricle','right atrium','right ventricle')
heart_direction <- c('positive','negative')

for(facid in seq(1,length(lianxu_factor))){ 
  use_fac <- lianxu_factor[facid]
  print(use_fac)
  m_factor = read.csv('modifiable_factors_meanbloodpressure.csv')
  m_factor_use <- m_factor[,c('eid',use_fac)]
  colnames(m_factor_use) <- c('eid','lifestyle')
  
  if (use_fac == 'X22189.0.0') {
    m_factor_use <- m_factor_use[!is.na(m_factor_use$lifestyle),]
  } else {
    m_factor_use <- m_factor_use[!is.na(m_factor_use$lifestyle) & m_factor_use$lifestyle!=-3 & m_factor_use$lifestyle!=-1,]
  }

  usedata_singlefac <- merge(usedata, m_factor_use, by=c('eid'))
  usedata_singlefac$lifestyle <- scale(usedata_singlefac$lifestyle)
  
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
            ordered = c(ill_use),
            estimator = "DWLS",
            se = "bootstrap",
            bootstrap = 1000
          )
          result <- summary(sem_result, fit.measures = TRUE, standardized=TRUE)
}

#######################################################################
use_fac <- "Education"
print(use_fac)
m_factor = read.csv('cova.csv')
m_factor_use <- m_factor[,c('eid', 'Education')]
colnames(m_factor_use) <- c('eid','lifestyle')
m_factor_use <- m_factor_use[!is.na(m_factor_use$lifestyle),]
usedata_singlefac <- merge(usedata, m_factor_use, by=c('eid'))
usedata_singlefac$lifestyle <- scale(usedata_singlefac$lifestyle)

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
          ordered = c(ill_use), 
          estimator = "DWLS", 
          se = "bootstrap",
          bootstrap = 1000
        )
        result <- summary(sem_result, fit.measures = TRUE, standardized=TRUE)
        
        CFI <- c(CFI, result[["fit"]][["cfi"]])
        RMSEA <- c(RMSEA, result[["fit"]][["rmsea"]])
        SRMR <- c(SRMR, result[["fit"]][["srmr"]])
      }
    }
  }
}
