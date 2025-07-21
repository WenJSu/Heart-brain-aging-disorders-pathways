library(dplyr)
library(readr)
library(caret)
library(glmnet)
library(Metrics)

# Define paths
cova_data_path <- "cova.csv"
cova_data <- read.csv(cova_data_path)
eid_data_path <- "BA_health_eid.csv"
allSiteInfo <- read.table('AllSiteInfo.csv', sep=',' ,header=T, check.names=F)
HeightInfo <- read.table("AllHeightInfo.csv", sep=',', header=T, check.names=F)
HeightInfo <- HeightInfo[,c('eid','50-0.0')]
colnames(HeightInfo) <- c('eid', 'Height')
# Define organs and their respective data files
organs <- c('Brain_WM')
########################## model traning and testing on normal controls
for (i in 1:length(organs)) {
  organ <- organs[i]
  
  # Load data
  df_GM <- read.csv("Brain_GM_clean.csv")
  df_WM <- read.csv("Brain_WM_clean.csv")
  df <- merge(df_GM, df_WM, by=c('eid'))
  
  # using only normal controls
  eid <- read.csv(eid_data_path)
  df <- inner_join(df, eid['eid'], by = 'eid')
  age <- read.csv(cova_data_path)
  age <- age[,c('eid','age2')]
  df <- inner_join(df, age, by = 'eid')
  df <- rename(df, Age = age2)

  # get confounding info
  df <- merge(df, allSiteInfo[,c('eid','54-2.0')])
  colnames(df)[which(colnames(df) == "54-2.0")] <- "Site"
  df$Site <- as.character(df$Site)
  df <- fastDummies::dummy_cols(df, select_columns = "Site", remove_first_dummy = TRUE)
  SiteCols <- unlist(sapply(colnames(df),function(x) if(grepl('Site_',x)) x),use.names = F)
  TIV <- cova_data
  TIV <- TIV[,c('eid','eTIV','Sex')]
  df <- merge(df, TIV, by = 'eid')
  df <- na.omit(df)
  table(df$Site)

  # Split data into training and testing sets
  set.seed(12345)
  folds <- createFolds(df$Age, k = 10)
  list_mae <- list()
  list_r <- list()
  list_models <- list()
  list_correctionModel <- list()
  list_lambda <- list()
  list_scaler <- list()
  list_plots <- list()
  list_coef <- list()
  list_covModles <- list()
  
  alldata <- df
  for (k in 1:10) {
    print(k)
    test_indices <- folds[[k]]
    train_indices <- setdiff(1:nrow(alldata), test_indices)
    
    x_train_data <- alldata[train_indices,]
    x_test_data <- alldata[test_indices,]
    data_matrix <- x_train_data
    data_matrix_test <- x_test_data
    
    residuals_matrix <- matrix(NA, nrow = nrow(data_matrix), ncol = ncol(data_matrix))
    residuals_matrix_test <- matrix(NA, nrow = nrow(data_matrix_test), ncol = ncol(data_matrix_test))
    colnames(residuals_matrix) <- colnames(data_matrix)
    colnames(residuals_matrix_test) <- colnames(data_matrix_test)
    covModels <- list()
    for (cc in 1:ncol(data_matrix)) {
      temp_data <- cbind(dependent_variable = data_matrix[, cc], x_train_data[,c('eTIV',SiteCols)])
      testing_data <- cbind(dependent_variable = data_matrix_test[, cc], x_test_data[,c('eTIV',SiteCols)])
      model <- lm(as.formula(paste0('dependent_variable ~ eTIV +',paste(SiteCols,collapse ='+'))), data = temp_data)
      residuals_matrix[, cc] <- residuals(model)
      # Predict the effect of covariates on the testing data using the training model
      predicted_values_test <- predict(model, newdata = testing_data)
      # Calculate residuals for the testing data
      residuals_matrix_test[,cc] <- testing_data$dependent_variable - predicted_values_test
      covModels[[cc]] <- model
    }
    list_covModles[[k]] <- covModels
    X_train_fold <- residuals_matrix
    y_train_fold <- alldata[train_indices, 'Age']
    X_test_fold <- residuals_matrix_test
    y_test_fold <- alldata[test_indices, 'Age']
    
    # normalizaion and standardziation
    library(caret)
    normParam <- preProcess(X_train_fold)
    X_train_fold <- predict(normParam, X_train_fold)
    X_test_fold <- predict(normParam, X_test_fold)
    list_scaler[[k]] <- normParam
    
    set.seed(1010)
    lasso_reg = glmnet(as.matrix(X_train_fold), y_train_fold, alpha = 1)
    cv_lasso <- cv.glmnet(as.matrix(X_train_fold), y_train_fold, lambda = 10^seq(5, -8, by = -.1), alpha = 1, type.measure = 'mae',nfolds = 10,parallel = T)
    plot(cv_lasso)
    list_plots[[k]] <- recordPlot()
    optimal_lambda <- cv_lasso$lambda.min
    lasso_model = predict(lasso_reg, type = "coefficients", s = optimal_lambda)
    lasso_model_nonzero <- lasso_model[(lasso_model[,1]!= 0),]
    list_lambda[[k]] <- optimal_lambda
    lasso_predict <- predict(lasso_reg, newx = as.matrix(X_test_fold), s = optimal_lambda)
    
    df[test_indices, 'pre_age'] <- as.numeric(lasso_predict)
    
    # Age bias correction
    train_predicted_y <- predict(lasso_reg, newx = as.matrix(X_train_fold), s = optimal_lambda)
    model <- lm(train_predicted_y ~ y_train_fold)
    test <- data.frame(TrueAge = as.numeric(y_test_fold), age_LASSO = as.numeric(lasso_predict))
    test$age_correct_LASSO <- (test$age_LASSO - coef(model)[1]) / coef(model)[2]
    
    df[test_indices, 'delta_corrected'] <- unlist(test$age_correct_LASSO - test$TrueAge)
    df[test_indices, 'pre_age_corrected'] <- test$age_correct_LASSO
    
    list_mae[[k]] <- mae(test$TrueAge, test$age_LASSO)
    list_r[[k]] <- cor(test$TrueAge, test$age_LASSO)
  }
  
  best_model_index <- which.min(list_mae)
  
  ##################### Model evaluation
  r <- cor(df$Age, df$pre_age)
  r_corrected <- cor(df$Age, df$pre_age_corrected)
  mae <- mae(df$Age, df$pre_age)
  mae_corrected <- mae(df$Age, df$pre_age_corrected)
  corrdata <- data.frame(df$Age, df$pre_age)
  write.table(corrdata, 'brain_age_prediction_evaluation.csv', sep=',', row.names=F)

  # Get weights of the predictors
  brain_traits <- data.frame(rownames(best_coef),as.matrix(best_coef))
  colnames(brain_traits) <- c('field','coef')
  write_csv(brain_traits, 'Age_LASSO_weights.csv')
  
  ##################### Apply the trained model to other subjects (excluding normal controls)
  df <- merge(df_GM, df_WM, by=c('eid'))
  
  # exclude only normal controls
  eid <- read.csv(eid_data_path)
  df <- df[!(df$eid %in% eid$eid),]
  cnt <- rowSums(!is.na(df))
  df <- df[cnt == ncol(df), ]

  # apply the covariate regression model to other subjs
  x_test_data <- df
  data_matrix_test <- x_test_data[,unlist(sapply(colnames(x_test_data),function(x) grepl('^X',x)), use.names = F)]
  # Initialize a matrix to store residuals
  residuals_matrix_test <- matrix(NA, nrow = nrow(data_matrix_test), ncol = ncol(data_matrix_test))
  colnames(residuals_matrix_test) <- colnames(data_matrix_test)
  for (cc in 1:ncol(data_matrix_test)) {
    testing_data <- cbind(dependent_variable = data_matrix_test[, cc], x_test_data[,c('eTIV',SiteCols)])
    # Predict the effect of covariates on the testing data using the training model
    predicted_values_test <- predict(best_covModel[[cc]], newdata = testing_data)
    # Calculate residuals for the testing data
    residuals_matrix_test[,cc] <- testing_data$dependent_variable - predicted_values_test
  }
  df[,unlist(sapply(colnames(df),function(x) grepl('^X',x)), use.names = F)] <- residuals_matrix_test
  df <- df[, !(names(df) %in% c(unlist(sapply(colnames(df),function(x) {if(grepl('Site',x)) x}), use.names = F),'eTIV','Height','Sex'))]
  
  
  # normalization and standardization
  df[,2:ncol(df)-1] <- predict(best_scaler, df[,2:ncol(df)-1])
  X_all <- df[, -c(1, ncol(df))]
  y_all_pred <- as.numeric(predict(best_model, newx = as.matrix(X_all), s = best_lambda))
  y_all_pred_corrected <- (y_all_pred - coef(best_correction_model)[1]) / coef(best_correction_model)[2]
  all_delta_corrected <- y_all_pred_corrected - df$Age
  df$pre_age <- y_all_pred
  df$pre_age_corrected <- y_all_pred_corrected
  df$delta_corrected <- all_delta_corrected
  df_all <- rbind(df_HC, df)
  write_csv(data.frame(df_all), 'brain_age_prediction_all.csv')
}
selected_df_agegap <-  df_all[, c("eid", "Age", "pre_age", "pre_age_corrected", "delta_corrected")]
write_csv(data.frame(selected_df_agegap), 'agegap.csv')


