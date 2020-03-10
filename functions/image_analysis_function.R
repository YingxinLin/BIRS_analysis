# image function
library(raster)
library(ggplot2)
library(ggthemes)

plotExpression <- function(raster_obj, sce_cellLabel, exprsMat,
                           gene, scale = TRUE) {
  
  
  ddf <- rasterToPoints(raster_obj)
  ddf <- data.frame(ddf)
  colnames(ddf) <- c("X", "Y")
  cell_label <- raster::values(r)
  
  notInLabel <- unique(cell_label[!cell_label %in% sce_cellLabel])
  
  # exprsMat <- assay(sce, exprs_value)
  
  new_values <- mapvalues((cell_label), 
                          from = notInLabel, 
                          to = rep(min(exprsMat[gene, ]),
                                   length(notInLabel)))
  
  
  new_values <- mapvalues((new_values), 
                          from = sce_cellLabel, 
                          to = exprsMat[gene, ])
  
  if (scale) {
    new_values <- zero_one_scale(new_values)
  }
  
  
  
  g <- ggplot(NULL) + 
    geom_raster(data = ddf, aes(X, Y, fill = new_values)) +
    theme_minimal() +
    scale_fill_viridis_c() +
    coord_quickmap() +
    theme(aspect.ratio = 1) +
    labs(fill = gene)
  g
  #return(g)
}


zero_one_scale <- function(mat) {
  (mat - min(mat))/(max(mat) - min(mat))
}

sig_pcf <- function(x, idx, E_pcf) {
  if (x > E_pcf$hi[idx]) {
    return(1)
  } else if (x < E_pcf$lo[idx]) {
    return(-1)
  } else {
    return(0)
  }
}

extract_pcf_results <- function(E_pcf, results_name = NULL) {
  
  if (all(is.na(E_pcf$obs))) {
    res <- rep(NA, 11)
    if (!is.null(results_name)) {
      names(res) <-  paste(results_name, c("max_pcf", "max_r", 
                                           "max_sig", "r50_pcf", "r50_pcf_sig","r100_pcf", "r100_pcf_sig",
                                           "r150_pcf", "r150_pcf_sig","r200_pcf", "r200_pcf_sig"
      ), sep = "_")
    } else {
      names(res) <-  c("max_pcf", "max_r", 
                                           "max_sig", "r50_pcf", "r50_pcf_sig","r100_pcf", "r100_pcf_sig",
                                           "r150_pcf", "r150_pcf_sig","r200_pcf", "r200_pcf_sig"
      )
    }

    return(res)
  }
  
  
  max_idx <- which.max(E_pcf$obs[!is.infinite(E_pcf$obs)])
  max_pcf <- E_pcf$obs[!is.infinite(E_pcf$obs)][max_idx]
  max_r <- E_pcf$r[!is.infinite(E_pcf$obs)][max_idx]
  max_sig <- sig_pcf(max_pcf, max_idx, E_pcf)
  
  r50_idx <- which.max(E_pcf$r == 50)
  r50_pcf <- E_pcf$obs[r50_idx]
  r50_pcf_sig <- sig_pcf(r50_pcf, r50_idx, E_pcf)
  
  r100_idx <- which.max(E_pcf$r == 100)
  r100_pcf <- E_pcf$obs[r100_idx]
  r100_pcf_sig <- sig_pcf(r100_pcf, r100_idx, E_pcf)
  
  r150_idx <- which.max(E_pcf$r == 150)
  r150_pcf <- E_pcf$obs[r150_idx]
  r150_pcf_sig <- sig_pcf(r150_pcf, r150_idx, E_pcf)
  
  r200_idx <- which.max(E_pcf$r == 200)
  r200_pcf <- E_pcf$obs[r200_idx]
  r200_pcf_sig <- sig_pcf(r200_pcf, r200_idx, E_pcf)
  
  res <- c(max_pcf = max_pcf, max_r = max_r, 
           max_sig = max_sig,
           r50_pcf = r50_pcf, r50_pcf_sig = r50_pcf_sig,
           r100_pcf = r100_pcf, r100_pcf_sig = r100_pcf_sig,
           r150_pcf = r150_pcf, r150_pcf_sig = r150_pcf_sig,
           r200_pcf = r200_pcf, r200_pcf_sig = r200_pcf_sig
  )
  
  if (!is.null(results_name)) {
    names(res) <- paste(results_name, names(res), sep = "_")
  }
  
  return(res)
  
  
}


runSVM <- function(X, y, K = 5, num_repeated = 100, ncores = 2) {
  
  n <- length(y)
  y <- droplevels(y)
  
  
  cv_pred_res <- pbmcapply::pbmclapply(1:num_repeated, function(i) {
    # if (i %% 10 == 0) { print(i) }
    cvSets = cvTools::cvFolds(n, K)  # permute all the data, into 5 folds
    
    cv_acc = NA  # initialise results vector
    predict_res <- c()
    for (j in 1:K) {
      test_id <- cvSets$subsets[cvSets$which == j]
      X_test <- X[test_id, ]
      X_train <- X[-test_id, ]
      y_test <- y[test_id]
      y_train <- y[-test_id]
      svm_res <- e1071::svm(x = X_train, 
                            y = as.factor(droplevels(y_train)),
                            kernel = "linear")
      fit <- predict(svm_res, X_test)
      pred <- as.character(fit)
      cv_acc[j] = sum(pred == as.character(y_test))/length(y_test)
      names(pred) <- names(fit)
      predict_res <- append(predict_res, pred)
    }
    
    predict_res <- predict_res[rownames(X)]
    
    return(list(predict_res = predict_res, accuracy = sum(predict_res == as.character(y))/length(y)))
  }, mc.cores = ncores)
  
  predict_res_mat <- do.call(rbind, lapply(cv_pred_res, function(x) x$predict_res))
  cv_acc <- do.call(c, lapply(cv_pred_res, function(x) x$accuracy))
  
  return(list(predict_res_mat = predict_res_mat, cv_acc = cv_acc))
}


runRF <- function(X, y, K = 5, num_repeated = 100, ncores = 2) {
  
  n <- length(y)
  y <- droplevels(y)
  
  
  cv_pred_res <- pbmcapply::pbmclapply(1:num_repeated, function(i) {
    # if (i %% 10 == 0) { print(i) }
    cvSets = cvTools::cvFolds(n, K)  # permute all the data, into 5 folds
    
    cv_acc = NA  # initialise results vector
    predict_res <- c()
    for (j in 1:K) {
      test_id <- cvSets$subsets[cvSets$which == j]
      X_test <- X[test_id, ]
      X_train <- X[-test_id, ]
      y_test <- y[test_id]
      y_train <- y[-test_id]
      rf_res <- randomForest::randomForest(x = X_train, y = as.factor(droplevels(y_train)))
      fit <- predict(rf_res, X_test)
      pred <- as.character(fit)
      cv_acc[j] = sum(pred == as.character(y_test))/length(y_test)
      names(pred) <- names(fit)
      predict_res <- append(predict_res, pred)
    }
    
    predict_res <- predict_res[rownames(X)]
    
    return(list(predict_res = predict_res, accuracy = sum(predict_res == as.character(y))/length(y)))
  }, mc.cores = ncores)
  
  predict_res_mat <- do.call(rbind, lapply(cv_pred_res, function(x) x$predict_res))
  cv_acc <- do.call(c, lapply(cv_pred_res, function(x) x$accuracy))
  
  return(list(predict_res_mat = predict_res_mat, cv_acc = cv_acc))
}


