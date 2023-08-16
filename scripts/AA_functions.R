### Some of these functions are adapted from Caetano et al. 2022 PLoS Biol. 20, e3001544. Please refer to that publication for the original, unaltered versions

library(tidyverse)

# function for even sampling across folds
even_folds <- function (pred, cat, CV_folds1, CV_folds2, seed) {
  set.seed(seed)  # Set the random seed for reproducibility
  nfolds <- CV_folds2 + 1  # Total number of folds (including the validation fold)
  class <- unique(pred$class)  # Unique classes in the 'pred' data frame
  biome <- unique(pred$biome)  # Unique biomes in the 'pred' data frame
  nclass <- length(class)  # Number of unique classes
  nbiome <- length(biome)  # Number of unique biomes
  exclude <- list()  # Initialize an empty list to store excluded samples
  folds <- vector("list", nfolds)  # Create a list to store the folds
  
  # Loop through each fold and initialize it as an empty list
  folds <- lapply(1:CV_folds1, function(z) folds)
  
  # Loop through each class and biome combination
  for (i in class) {
    for (j in biome) {
      # Subset the 'pred' data frame based on class, biome, and 'cat' value
      index <- which(((pred$class == i) * (pred$biome == j) * (pred[, cat] == "no")) == 1)
      nindex <- length(index)  # Number of samples in the subset
      
      if (nindex < CV_folds1 && nindex > 0) {
        # If the number of samples is less than the desired number of folds but greater than zero,
        # create a temporary list and add it to 'exclude' list with a specific name
        tmp <- list(a = index)
        names(tmp) <- paste0(i, ".", j, ".no")
        exclude <- c(exclude, tmp)
      } else {
        # If there are enough samples, randomly sample from them without replacement
        index <- sample(index, size = nindex, replace = FALSE)
        n <- floor(nindex / CV_folds1)  # Number of samples per fold
        
        if (n > 0) {
          # Distribute the samples evenly across the folds
          for (z in 1:CV_folds1) {
            folds[[z]][[1]] <- c(folds[[z]][[1]], index[(z - 1) * n + c(1:n)])
          }
        }
        
        n2 <- nindex - n * CV_folds1  # Remaining samples
        if (n2 > 0) {
          # Randomly assign the remaining samples to a fold
          tmp <- sample(1:CV_folds1, size = n2, replace = FALSE)
          for (z in 1:n2) {
            folds[[tmp[z]]][[1]] <- c(folds[[tmp[z]]][[1]], index[CV_folds1 * n + z])
          }
        }
        
        for (z in 1:CV_folds1) {
          remain <- index[!index %in% folds[[z]][[1]]]  # Remaining samples after the initial fold assignment
          nremain <- length(remain)  # Number of remaining samples
          remain <- sample(remain, size = nremain, replace = FALSE)
          n <- floor(nremain / CV_folds2)  # Number of samples per fold in the second level of folding
          
          if (n > 0) {
            # Distribute the remaining samples evenly across the folds (second level of folding)
            for (a in 2:(CV_folds2 + 1)) {
              folds[[z]][[a]] <- c(folds[[z]][[a]], remain[(a - 2) * n + c(1:n)])
            }
          }
          
          n2 <- nremain - n * CV_folds2  # Remaining samples in the second level of folding
          if (n2 > 0) {
            # Randomly assign the remaining samples to a fold (second level of folding)
            tmp <- sample(2:(CV_folds2 + 1), size = n2, replace = FALSE)
            for (a in 1:n2) {
              folds[[z]][[tmp[a]]] <- c(folds[[z]][[tmp[a]]], remain[CV_folds2 * n + a])
            }
          }
        }
      }
      
      # Repeat the same process for samples with 'cat' value equal to "yes"
      index <- which(((pred$class == i) * (pred$biome == j) * (pred[, cat] == "yes")) == 1)
      nindex <- length(index)
      
      if (nindex < CV_folds1 && nindex > 0) {
        tmp <- list(a = index)
        names(tmp) <- paste0(i, ".", j, ".yes")
        exclude <- c(exclude, tmp)
      } else {
        index <- sample(index, size = nindex, replace = FALSE)
        n <- floor(nindex / CV_folds1)
        
        if (n > 0) {
          for (z in 1:CV_folds1) {
            folds[[z]][[1]] <- c(folds[[z]][[1]], index[(z - 1) * n + c(1:n)])
          }
        }
        
        n2 <- nindex - n * CV_folds1
        if (n2 > 0) {
          tmp <- sample(1:CV_folds1, size = n2, replace = FALSE)
          for (z in 1:n2) {
            folds[[tmp[z]]][[1]] <- c(folds[[tmp[z]]][[1]], index[CV_folds1 * n + z])
          }
        }
        
        for (z in 1:CV_folds1) {
          remain <- index[!index %in% folds[[z]][[1]]]
          nremain <- length(remain)
          remain <- sample(remain, size = nremain, replace = FALSE)
          n <- floor(nremain / CV_folds2)
          
          if (n > 0) {
            for (a in 2:(CV_folds2 + 1)) {
              folds[[z]][[a]] <- c(folds[[z]][[a]], remain[(a - 2) * n + c(1:n)])
            }
          }
          
          n2 <- nremain - n * CV_folds2
          if (n2 > 0) {
            tmp <- sample(2:(CV_folds2 + 1), size = n2, replace = FALSE)
            for (a in 1:n2) {
              folds[[z]][[tmp[a]]] <- c(folds[[z]][[tmp[a]]], remain[CV_folds2 * n + a])
            }
          }
        }
      }
    }
  }
  
  # Return the 'exclude' list and the 'folds' list as the output of the function
  list(exclude = exclude, folds = folds)
}


# function for hyperparameter tuning (using 5-fold cross-validation)
hp_tuning <- function(X2, y2, folds2, interaction_constraints=NULL, seed) {
  # set seed for reproducibility
  set.seed(seed)
  
  # Create 10,000 random samples of hyperparameter values
  par_list <- list()
  
  for (iter in 1:10000) {
    # Generate random values for each hyperparameter
    par_list[[iter]] <- list(
      booster = "gbtree",
      objective = "binary:logistic",
      eta = runif(1, .01, .3),                      # Learning rate
      max_depth = sample(3:10, 1),                   # Maximum tree depth
      subsample = runif(1, .7, 1),                   # Subsample ratio of the training instances
      colsample_bytree = runif(1, .6, 1),            # Subsample ratio of columns when constructing each tree
      min_child_weight = runif(1, .5, 1.5),          # Minimum sum of instance weight needed in a child
      scale_pos_weight = runif(1, .5, 1.5),          # Control the balance of positive and negative weights
      gamma = runif(1, 0, 1),                        # Minimum loss reduction required to make a further partition on a leaf node of the tree
      alpha = runif(1, 0, 1),                        # L1 regularization term on weights
      lambda = runif(1, .5, 1.5)                     # L2 regularization term on weights
    )
  }
  
  # Run 10,000 parameter sets in parallel
  cl <- parallel::makeCluster(7)
  parallel::clusterExport(cl, c("X", "y", "par_list", "folds2", "interaction_constraints"), envir = environment())
  
  res_list <- pbapply::pblapply(cl = cl, 1:length(par_list), function(j) {
    set.seed(1)
    
    # Create empty variable to store error rate
    err <- NULL
    
    for (i in 1:length(folds2)) {
      # Partition the data into training and validation sets
      X_train <- X[unlist(folds2[-i]), ]
      y_train <- y[unlist(folds2[-i])]
      X_val <- X[unlist(folds2[i]), ]
      y_val <- y[unlist(folds2[i])]
      
      # Convert the data to xgb.DMatrix format
      if (is.data.frame(X_train)) {
        dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X_train), label = y_train, missing = NA)
        dval <- xgboost::xgb.DMatrix(data = as.matrix(X_val), label = y_val, missing = NA)
      } else if (is.matrix(X_train) || (inherits(X_train, 'dgCMatrix'))) {
        dtrain <- xgboost::xgb.DMatrix(data = X_train, label = y_train, missing = NA)
        dval <- xgboost::xgb.DMatrix(data = X_val, label = y_val, missing = NA)
      } else {
        stop(simpleError("X must be either a data.frame or a (sparse-) matrix"))
      }
      
      # Configure the XGBoost parameters
      params_xgboost <- par_list[[j]]
      params_xgboost[['watchlist']] <- list(train = dtrain, validation = dval)
      params_xgboost[['data']] <- dtrain
      params_xgboost[['interaction_constraints']] <- interaction_constraints
      params_xgboost[['print_every_n']] <- 10
      params_xgboost[['early_stopping_rounds']] <- 10
      params_xgboost[['maximize']] <- FALSE
      params_xgboost[['eval_metric']] <- 'error'
      params_xgboost[['nrounds']] <- 200
      params_xgboost[['prediction']] <- TRUE
      
      # Train the XGBoost model with the current parameter set
      bst <- suppressWarnings(do.call(xgboost::xgb.train, params_xgboost))
      
      # Calculate the classification error rate
      err <- 1 - bst$evaluation_log[bst$best_iteration, 3]
      
      gc()  # Perform garbage collection to free up memory
      
    }
    
    list(err = err)
  })
  
  parallel::stopCluster(cl)
  
  # Return the hyperparameter set with the highest predictive capability based on binary classification error rate
  acc <- sapply(res_list, function(i) i$err)
  idx <- which.max(acc)[1]
  param <- par_list[[idx]]
  
  list(param = param)
}


# function for automatic assessment:
auto_ass <- function(X, y, path, folds, seed, interaction_constraints=NULL) {
  # set seed for reproducibility
  set.seed(seed)
  
  # Initialize variables and containers
  xgb.fit <- list()
  param.list <- list()
  CV_folds1 <- length(folds[[2]])
  err <- numeric(CV_folds1)
  nexclude <- length(folds[[1]])
  err.exclude <- matrix(NA,CV_folds1,nexclude)
  colnames(err.exclude) <- names(folds[[1]])
  
  # Perform cross-validation
  for (i in 1:CV_folds1) {
    # Split data into fitting and testing sets
    X_fit <- X[unlist(folds[[2]][[i]][-1]),]
    y_fit <- y[unlist(folds[[2]][[i]][-1])]
    X_test <- X[c(unlist(folds[[1]]),unlist(folds[[2]][[i]][1])),]
    y_test <- y[c(unlist(folds[[1]]),unlist(folds[[2]][[i]][1]))]
    folds2 <- folds[[2]][[i]][-1]
    
    # Perform hyperparameter tuning using the fitting set
    out <- hp_tuning(X=X, y=y, folds2=folds2, interaction_constraints=interaction_constraints, seed=seed)
    param.list[[i]] <- out$param
    
    # Convert data to xgb.DMatrix format
    dfit <- xgboost::xgb.DMatrix(data=X_fit, label=y_fit, missing=NA)
    dtest <- xgboost::xgb.DMatrix(data=X_test, label=y_test, missing=NA)
    
    # Prepare excluded data for evaluation (if applicable)
    if (nexclude > 0) {
      dexclude <- list()
      for (j in 1:nexclude) {
        X_exclude <- X[folds[[1]][[j]],]
        if (!is.matrix(X_exclude)) {
          X_exclude <- matrix(X_exclude, nrow=1)
        }
        y_exclude <- y[folds[[1]][[j]]]
        tmp <- list(a = xgboost::xgb.DMatrix(data=X_exclude, label=y_exclude, missing=NA))
        names(tmp) <- names(folds[[1]])[j]
        dexclude <- c(dexclude, tmp)
      }
    }
    
    # Set the parameters for xgboost training
    params_xgboost <- param.list[[i]]
    if (length(folds[[1]]) > 0) {
      params_xgboost[['watchlist']] <- c(dexclude, list(train=dfit, validation=dtest))
    } else {
      params_xgboost[['watchlist']] <- list(train=dfit, validation=dtest)
    }
    params_xgboost[['data']] <- dfit
    params_xgboost[['interaction_constraints']] <- interaction_constraints
    params_xgboost[['print_every_n']] <- 10
    params_xgboost[['early_stopping_rounds']] <- 10
    params_xgboost[['maximize']] <- FALSE
    params_xgboost[['eval_metric']] <- 'error'
    params_xgboost[['nrounds']] <- 5000
    
    # Train the xgboost model with the chosen parameters
    xgb.fit[[i]] <- suppressWarnings(do.call(xgboost::xgb.train, params_xgboost))
    
    # Calculate the error rate on the validation dataset
    tmp <- as.matrix(xgb.fit[[i]]$evaluation_log)
    err[i] <- tmp[xgb.fit[[i]]$best_iteration, nexclude+3]
    
    if (nexclude > 0) {
      err.exclude[i,] <- tmp[xgb.fit[[i]]$best_iteration, 2:(nexclude+1)]
    }
  }
  
  # Compute the average hyperparameter values chosen in folds2 for the final model
  bst.param <- list(
    booster = "gbtree",
    objective = "binary:logistic",
    eta = mean(sapply(param.list, function(i) i$eta)),
    max_depth = round(mean(sapply(param.list, function(i) i$max_depth))),
    subsample = mean(sapply(param.list, function(i) i$subsample)),
    colsample_bytree = mean(sapply(param.list, function(i) i$colsample_bytree)),
    min_child_weight = mean(sapply(param.list, function(i) i$min_child_weight)),
    scale_pos_weight = mean(sapply(param.list, function(i) i$scale_pos_weight)),
    gamma = mean(sapply(param.list, function(i) i$gamma)),
    alpha = mean(sapply(param.list, function(i) i$alpha)),
    lambda = mean(sapply(param.list, function(i) i$lambda))
  )
  
  # Convert the entire dataset to xgb.DMatrix format
  dall <- xgboost::xgb.DMatrix(data=X, label=y, missing=NA)
  
  # Set the parameters for the final xgboost model
  params_xgboost <- bst.param
  params_xgboost[['watchlist']] <- list(train=dall)
  params_xgboost[['data']] <- dall
  params_xgboost[['interaction_constraints']] <- interaction_constraints
  params_xgboost[['print_every_n']] <- 10
  params_xgboost[['early_stopping_rounds']] <- 10
  params_xgboost[['maximize']] <- FALSE
  params_xgboost[['eval_metric']] <- 'error'
  params_xgboost[['nrounds']] <- 5000
  params_xgboost[['prediction']] <- TRUE
  
  # Train the final xgboost model
  xgb.final <- suppressWarnings(do.call(xgboost::xgb.train, params_xgboost))
  
  # Write feature importance to a CSV file
  importance <- data.frame(xgboost::xgb.importance(colnames(X), model=xgb.final))
  write.csv(importance, file=paste0(path, seed, cat, "_importance.csv"))
  
  # Generate and export partial dependence plots
  pdp <- lapply(xgb.final$feature_names, function(x) pdp::partial(
    xgb.final,
    pred.var = x,
    ice = T,
    center = T,
    plot = T,
    alpha = .1,
    plot.engine = "ggplot2",
    train = X[, xgb.final$feature_names]
  ))
  
  ggsave(
    paste0(path, seed, cat, ".pdf"),
    cowplot::plot_grid(plotlist = pdp),
    scale = 1,
    dpi = 300
  )
  
  # Plot individual trees and ensemble tree
  trees_plot <- xgb.plot.tree(model = xgb.final, render = FALSE)
  DiagrammeR::export_graph(trees_plot, paste0(path, seed, cat, "_trees_plot.pdf"))
  
  ensemble_plot <- xgb.plot.multi.trees(model = xgb.final, render = FALSE)
  DiagrammeR::export_graph(ensemble_plot, paste0(path, seed, cat, "_ensemble_plot.pdf"))
  
  # Return the results as a list
  return(list(
    xgb.fit = xgb.fit,
    param.list = param.list,
    err = err,
    err.exclude = err.exclude,
    bst.param = bst.param,
    xgb.final = xgb.final
  ))
  
  # xgb.fit: Models fit in folds1 with hyperparameters chosen in folds2
  # param.list: Chosen hyperparameters in folds2
  # err: Error in the validation dataset in folds1. Used to assess model fit
  # err.exclude: Error of fitted models in folds1 to the biome_class_table with too few species
  # bst.param: Hyperparameters used for the final model
  # xgb.final: Final model fit to all species. Used for prediction
}


# function for prediction
pred_ass <- function(X, model.threatened = NULL, model.nt = NULL, model.en_cr = NULL, model.cr = NULL) {
  # Predict threatened vs. non-threatened species
  
  if (is.null(model.threatened)) {
    error("Need to provide a model to predict threatened species")
  } else {
    # Get the number of samples
    n <- dim(X)[1]
    
    # Create DMatrix for prediction
    dpred <- xgboost::xgb.DMatrix(data = X, missing = NA)
    
    # Predict threatened species
    tmp <- predict(model.threatened, dpred)
    ind.unthreatened <- which(tmp <= 0.5)
    ind.threatened <- which(tmp > 0.5)
    
    # Initialize category as "LC" (Least Concern)
    cat <- rep("LC", n)
    
    # Predict non-threatened species (NT)
    if (!is.null(model.nt) && length(ind.unthreatened) > 0) {
      dpred <- xgboost::xgb.DMatrix(data = X[ind.unthreatened, ], missing = NA)
      tmp <- predict(model.nt, dpred)
      cat[ind.unthreatened[tmp > 0.5]] <- "NT"
    }
    
    # Predict endangered (EN), vulnerable (VU), and critically endangered (CR) species
    if (!is.null(model.en_cr) && length(ind.threatened) > 0) {
      dpred <- xgboost::xgb.DMatrix(data = X[ind.threatened, ], missing = NA)
      tmp <- predict(model.en_cr, dpred)
      
      # Set category as EN for predicted endangered species
      cat[ind.threatened[tmp > 0.5]] <- "EN"
      
      # Set category as VU for predicted vulnerable species
      cat[ind.threatened[tmp <= 0.5]] <- "VU"
      
      ind.en_cr <- which(tmp > 0.5)
      
      # Predict critically endangered species (CR)
      if (!is.null(model.cr) && length(ind.en_cr) > 0) {
        dpred <- xgboost::xgb.DMatrix(data = X[ind.threatened[ind.en_cr], ], missing = NA)
        tmp <- predict(model.cr, dpred)
        cat[ind.threatened[ind.en_cr[tmp > 0.5]]] <- "CR"
      }
    }
  }
  
  # Return the predicted categories
  return(cat)
}