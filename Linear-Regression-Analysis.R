multRegData = read.table("MultRegData.txt",header = TRUE,sep = " ")


fmax <- function(data, na_rm = TRUE) {
  
  if (!is.numeric(data)) {
    stop("Input data must be numeric.")
  }
  
  if (length(data) == 0) {
    return(NA)
  }
  
  data_vec <- unlist(data)
  
  if (na_rm) {
    data_vec <- data_vec[!is.na(data_vec)]
  }
  
  if (length(data_vec) == 0) {
    return(NA)
  }
  
  max_val <- data_vec[1]
  
  for (i in 2:length(data_vec)) {
    if (data_vec[i] > max_val) {
      max_val <- data_vec[i]
    }
  }
  
  return(max_val)
}

fmin <- function(data, na_rm = TRUE) {
  
  if (!is.numeric(data)) {
    stop("Input data must be numeric.")
  }
  
  if (length(data) == 0) {
    return(NA)
  }
  
  data_vec <- unlist(data)
  
  if (na_rm) {
    data_vec <- data_vec[!is.na(data_vec)]
  }
  
  if (length(data_vec) == 0) {
    return(NA)
  }
  
  min_val <- data_vec[1]
  
  for (i in 2:length(data_vec)) {
    if (data_vec[i] < min_val) {
      min_val <- data_vec[i]
    }
  }
  
  return(min_val)
}

fmedian <- function(x) {
  n <- length(x)
  x <- x[!is.na(x)]  # Remove NA values
  n <- length(x)
  sorted_x <- sort(x)
  if (n %% 2 == 0) {
    mid <- n / 2
    median_value <- (sorted_x[mid] + sorted_x[mid + 1]) / 2
  } else {
    median_value <- sorted_x[(n + 1) / 2]
  }
  return(median_value)
}

fsum <- function(x) {
  total <- 0
  count <- 0
  for (i in 1:length(x)) {
    if (!is.na(x[i])) {
      total <- total + x[i]
      count <- count + 1
    }
  }
  return(total)
}

fmean <- function(data,na.rm=TRUE){
  if (!is.numeric(data)) {
    stop("Input must be a numeric vector")
  }
  if(na.rm){
    data <- data[!is.na(data)]
  }
  return(fsum(data)/length(data))
}

linear_regression <- function(y,...){
  
  x_list <- list(...)
  x <- do.call(cbind, x_list)
  
  x <- as.matrix(x)
  
  x <- cbind(1,x)
  
  beta <- solve(t(x) %*% x) %*% t(x) %*% y
  
  y_pred <- x %*% beta
  
  residuals <- y - y_pred
  
  mean_y <- fmean(y)
  
  tss <- fsum((y-mean_y)^2)
  rmss <- fsum((y_pred-mean_y)^2)
  rss <- fsum((y-y_pred)^2)
  r_squared <- rmss/tss
  
  return(list(
    beta=beta,
    residuals=residuals,
    ypred=y_pred,
    totalSumOfSquare=tss,
    regressionModelOfSquare=rmss,
    residualSumOfSquare=rss,
    rSquared=r_squared
  ))
}


model_selection <- function(y, ...) {
  x_list <- list(...)
  x_matrix <- do.call(cbind, x_list)
  n <- ncol(x_matrix)
  
  combinations <- unlist(lapply(1:n, function(m) {
    combn(1:n, m, simplify=FALSE)
  }), recursive=FALSE)
  
  models <- list()
  
  for (comb in combinations) {
    x_subset <- x_matrix[, comb, drop=FALSE]
    result <- linear_regression(y, x_subset)
    models[[length(models) + 1]] <- list(
      variable_name=paste("X",comb,collapse=",",sep=""),
      TSS=result$totalSumOfSquare,
      RMSS=result$regressionModelOfSquare,
      RSS=result$residualSumOfSquare,
      R_Square=result$rSquare
    )
  }
  
  models_sorted <- models[order(sapply(models, function(m) m$R_Square), decreasing=TRUE)]
  
  number_of_variables <- c()
  variable_name <- c()
  TSS <- c()
  RMSS <- c()
  RSS <- c()
  R_Square <- c()
  
  for (model in models_sorted) {
    number_of_variables <- c(number_of_variables, sapply(strsplit(model$variable_name, ","), length))
    variable_name <- c(variable_name, model$variable_name)
    TSS <- c(TSS, model$TSS)
    RMSS <- c(RMSS, model$RMSS)
    RSS <- c(RSS, model$RSS)
    R_Square <- c(R_Square, model$R_Square)
  }
  
  
  models_ <- cbind(number_of_variables,variable_name,TSS,RMSS,RSS,R_Square)
  
  return(models_)
}
