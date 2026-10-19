# density of survival times with censoring
f_td <-  function(shape, scale, t, d) {
  ((dweibull(t,shape=shape, scale = scale)^d)) * (pweibull(t, shape = shape, scale =scale,lower.tail = F)^((1-d))) 
}


survform=~CEP+male
# gammas_params <- theta$gamma_matrix
density_weibull <- function(shape_params, scale_params, data_wide) {
  
  N <- nrow(data_wide)
  k <- nrow(as.matrix(shape_params))
  
  
  f_tds <- matrix(0,nrow = N, ncol = k)
  # w <- model.matrix(~ -1+X3+X4+X5, data=data_wide)
  delta <- data_wide$delta
  # initial density function for Weibull
  
  for (j in 1:k) {
    # j <- 1
    shape <- shape_params[j]
    scale <- scale_params[j]
    # linear_pred <- exp()
    delta <- data_wide$delta
    # fixed effects in survival component
    # wxi <- as.matrix(w)
    ti <- data_wide$event_time
    for(i in 1:N){
      f_tds[i, j] <- f_td(shape, scale, ti[i], delta[i])
    }
  }
  return(f_tds)
}

# ftds <- density_weibull(shape_params = theta$shape, scale_params = theta$scale, gammas_params = gammas_params, data_wide, ~CEP+male)

# data_long <- datalong
marginal_density_yis <- function(data_long,XB_in, theta, long_form_g,ranform) {
  # parms prep
  betas_matrix <- theta$betas
  XB <- theta$XB
  sigmas <- theta$sigma_e
  # D_matrices <- theta$D
  
  # data prep
  y<- model.frame(long_form_g,data_long)[,1]
  Xg <- model.matrix(long_form_g, data_long)
  Zg <- model.matrix(ranform, data_long)
  data_long$y <- y
  
  N <- length(unique(data_long$ID)) #;print(N) # number of unique subjects under study.
  m <- nrow(betas_matrix) #;print(m);print(betas_matrix) # number of unique groups for mixture model
  f_yi <- matrix(0,nrow = N, ncol = m)
  nis<-table(data_long$ID)
  IDs<-unique(data_long$ID)
  # length(nis)
  
  D_matrices <- tcrossprod(theta$cholesky);# print(D_matrices)
  for (k in 1:m) {
    sigma_g <- sigmas # error variances 
    D_g <- D_matrices # random effect variance-covariance matrices
    beta_g <- betas_matrix[k,] # fixed effects
    
    for(i in seq_len(length(IDs))){
      
      # i<-2
      # n_i <- nis[i]
      # k <- 1
      ids <- IDs[i] 
      # data_long <- data;i <- 1
      subset <- data_long[data_long$ID == ids, ] #! ensure that IDs are consistently named across datasets
      # subset<-sim[sim$ID==i,]
      Z_ig <- Zg[data_long$ID == ids, , drop = FALSE]    # Random Effects Design Matrix
      X_ig <- Xg[data_long$ID==ids,,drop=FALSE ]# Cluster-specific Fixed Effects Design Matrix
      
      p <- ncol(X_ig)
      Xi_beta_k <- matrix(XB_in[data_long$ID==ids,k,drop=FALSE]) # Fixed Effects Linear Predictor
      # Fixed Effects Linear Predictor
      
      q <- ncol(Z_ig)
      Y_i <- subset$y # Longitudinal recordings
      
      #? calculating the distributional parameters for the group specific random effects
      # D_g_inv <- solve(D_g)
      # 3-by-2 dimensions
      V_ig <- as.matrix(Z_ig %*% (D_g) %*% t(Z_ig) + sigma_g^2 * diag(nrow(X_ig)))
      # simgas[i] <- list(Sigma_i)
      f_yi[i, k] <- mvnfast::dmvn(t(Y_i), Xi_beta_k, V_ig)
    }
  }
  return(f_yi)
}



plot_long <- function(fitted_nnjlcm, type=c("marg", "cond", "both"), ylabs="Longitudinal Biomarker", plotSurvScores=T){
  
  
  # extract data from fitted object 
  datalong <- fitted_nnjlcm$data_long
  data_wide <- fitted_nnjlcm$data_wide
  var.time <- "measure_time"
  X <- fitted_nnjlcm$fixef_design
  P <- fitted_nnjlcm$postprob
  theta <- fitted_nnjlcm$theta
  ng <- k_tofit
  tt <- datalong[[var.time]]
  y  <- fitted_nnjlcm$data_long$y
  IDs <- fitted_nnjlcm$IDs
  Z <- X
  ranef_est <- fitted_nnjlcm$ranef
  break.times <- quantile(tt, prob = seq(0, 1, length.out = 10))
  ntps <- length(break.times) - 1
  # 
  idx <- match(datalong$ID, data_wide$ID)
  if (anyNA(idx)) stop("some datalong IDs are not in data_wide")
  W  <- P[idx, , drop = FALSE]                        # n_obs x ng
  XB <- theta$XB                         # n_obs x ng
  
  Cl   <- max.col(P, ties.method = "first")[idx]
  maxT <- vapply(1:ng, function(g) if (any(Cl == g)) max(tt[Cl == g]) else -Inf, 0)
  bin <- findInterval(tt, break.times, rightmost.closed = TRUE)
  
  # keep[, g]: observation counts for class g (same truncation rule at maxT[g])
  keep <- sapply(1:ng, function(g)
    tt < maxT[g] | (bin == ntps & tt == maxT[g] & maxT[g] > break.times[ntps]))
  Wk <- W * keep
  
  # sums per bin, filled into all ntps rows
  binsum <- function(M) {
    M <- as.matrix(M) + 0
    out <- matrix(0, ntps, ncol(M))
    r <- rowsum(M, bin)
    out[as.integer(rownames(r)), ] <- r
    out
  }
  # # print(dim(Wk)); print(dim(y))
  y <-  as.numeric(y)
  sw  <- binsum(Wk);      swp <- binsum(Wk * XB);
  swy <- binsum(Wk *y);    sw2 <- binsum(Wk^2)
  n   <- binsum(keep);      sy  <- binsum(keep * y);  syy <- binsum(keep * y^2)
  
  sdy <- sqrt((syy - sy^2 / n) / (n - 1)); sdy[n < 2] <- NA   # unweighted sd(y), as before
  pred <- swp / sw
  obs  <- swy / sw
  se   <- sdy * sqrt(sw2) / sw
  
  pred_times0 <- unname(cbind((break.times[-length(break.times)] + break.times[-1]) / 2,
                              pred, obs, obs - 1.96 * se, obs + 1.96 * se))
  pred_times0[rowSums(n) == 0, -1] <- NA;
  colnames(pred_times0) <- c(var.time, paste0("pred.class", 1:ng), paste0("obs.class", 1:ng),
                             paste0("lower.obs.class", 1:ng), paste0("upper.obs.class", 1:ng))
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  # Adjust the time column name if it's not literally `meas`
  pred_long <- as.data.frame(pred_times0) |>
    pivot_longer(
      cols = -measure_time,
      names_to = c(".value", "class"),
      names_pattern = "(.+)\\.(class\\d+)"
    ) |>
    filter(!is.na(pred))   # drop the rows where pred/obs are NA so lines don't break
  
  p_lme <- ggplot(pred_long, aes(x = measure_time, colour = class, fill = class)) +
    geom_ribbon(aes(ymin = lower.obs, ymax = upper.obs),
                alpha = 0.2, colour = NA) +
    geom_point(aes(y = obs)) +
    geom_line(aes(y = pred)) +
    facet_grid(~class)+
    labs(x = "Measurement time", y = ylabs,
         colour = "Class", fill = "Class") +
    theme_minimal()
  if(type=="marginal"){
    print(p_lme)
  }
  yhat <- theta$yhat
  rownames(P) <-data_wide$ID
  preds <- rowSums(yhat * P[(datalong$ID),])
  local_MSE_y_in <- Metrics::mse(actual = y, predicted = preds)
  #   
  predinfo <- data.frame(
    MSE = local_MSE_y_in
    # IBS = ibs_score
  )
  #   
  print(predinfo, row.names=FALSE)
  var.time <- "measure_time"
  ng <- k_tofit
  tt <- datalong[[var.time]]
  y  <- datalong$y
  
  break.times <- quantile(tt, prob = seq(0, 1, length.out = 10))
  ntps <- length(break.times) - 1
  
  # # posterior weights per observation, without merge()
  idx <- match(datalong$ID, data_wide$ID)
  if (anyNA(idx)) stop("some datalong IDs are not in data_wide")
  W  <- P[idx, , drop = FALSE]                        # n_obs x ng
  XB <- yhat                         # n_obs x ng
  
  # # max follow-up time per modal class (same as your maxT)
  Cl   <- max.col(P, ties.method = "first")[idx]
  maxT <- vapply(1:ng, function(g) if (any(Cl == g)) max(tt[Cl == g]) else -Inf, 0)
  
  # # time bin of each observation; last bin closed on the right
  bin <- findInterval(tt, break.times, rightmost.closed = TRUE)
  
  # # keep[, g]: observation counts for class g (same truncation rule at maxT[g])
  keep <- sapply(1:ng, function(g)
    tt < maxT[g] | (bin == ntps & tt == maxT[g] & maxT[g] > break.times[ntps]))
  Wk <- W * keep
  
  # # sums per bin, filled into all ntps rows
  binsum <- function(M) {
    M <- as.matrix(M) + 0
    out <- matrix(0, ntps, ncol(M))
    r <- rowsum(M, bin)
    out[as.integer(rownames(r)), ] <- r
    out
  }
  y <- as.numeric(y)
  sw  <- binsum(Wk);        swp <- binsum(Wk * XB);
  swy <- binsum(Wk * y);    sw2 <- binsum(Wk^2)
  n   <- binsum(keep);      sy  <- binsum(keep * y);  syy <- binsum(keep * y^2)
  
  sdy <- sqrt((syy - sy^2 / n) / (n - 1)); sdy[n < 2] <- NA   # unweighted sd(y), as before
  pred <- swp / sw
  obs  <- swy / sw
  se   <- sdy * sqrt(sw2) / sw
  
  pred_times0 <- unname(cbind((break.times[-length(break.times)] + break.times[-1]) / 2,
                              pred, obs, obs - 1.96 * se, obs + 1.96 * se))
  pred_times0[rowSums(n) == 0, -1] <- NA
  colnames(pred_times0) <- c(var.time, paste0("pred.class", 1:ng), paste0("obs.class", 1:ng),
                             paste0("lower.obs.class", 1:ng), paste0("upper.obs.class", 1:ng))
  
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  # Adjust the time column name if it's not literally `meas`
  pred_long <- as.data.frame(pred_times0) |>
    pivot_longer(
      cols = -measure_time,
      names_to = c(".value", "class"),
      names_pattern = "(.+)\\.(class\\d+)"
    ) |>
    filter(!is.na(pred))   # drop the rows where pred/obs are NA so lines don't break
  
  p_lme2 <- ggplot(pred_long, aes(x = measure_time, colour = class, fill = class)) +
    geom_ribbon(aes(ymin = lower.obs, ymax = upper.obs),
                alpha = 0.2, colour = NA) +
    geom_point(aes(y = obs)) +
    geom_line(aes(y = pred)) +
    facet_grid(~class)+
    labs(x = "Measurement time", y = ylabs,
         colour = "Class", fill = "Class") +
    theme_minimal()
  if(type=="conditional"){
    print(p_lme2)
  }
  
  if(type=="both"){
    library(patchwork)
    item <- p_lme + p_lme2 + plot_layout(ncol = 1, heights = c(1, 1)) + plot_annotation(title = paste0("Train Set: Weighted Marginal Predictions vs Weighted Subject-Specific Predictions for k=",k_tofit), subtitle = "Weighted by posterior probabilities of class membership", caption = "Solid = predicted, dashed = observed, ribbon = obs 95% CI")
    print(item)
  }
  
}

# # example usage
# ysu <- marginal_density_yis(
#   data_long = datalong, 
# theta=theta,
# long_form_g = long_form_g,
# long_form_o = long_form_o,
# ranform = ranform)

#' IBS with Parallel Processing
#'
#' IBS is  an integrated version of the Brier which is used to calculate the integration of the Brier Score.
#' The Brier Score is the mean square difference between the true classes and the predicted probabilities.
#' Basically, the IBS is an integrated
#' weighted squared distance between the estimated survival function and the empirical
#' survival function. The inverse probability censoring weighting(IPCW) is used to adjust for censoring.
#'
#' The percentage of censored observations increases in time, and this will surely affect the dispersion of the empirical Brier Score.
#' The question of how censoring in finite samples acts on the distribution of our measures of inaccuracy is an interesting subject.
#' Our recommendation is to choose t* in a way that censoring is not too heavy (for example, the median follow-up time).
#' We also prefer measures with integrated loss functions since they will reflect inaccuracy over an interval rather than just at one point in time.
#' In addition, the corresponding empirical measures are likely to have lower dispersion, because censored observations contribute their estimated event-free probabilities to the integrand until the censoring occurs.
#'
#' @param object object of class \code{Surv} in the testing set created by Surv function.
#' @param sp_matrix a matrix or data.frame of predicted values of survival probabilities for the testing set.
#' @param IBSrange a vector contains all discrete time points corresponding to the predicted probability in sp_matrix.
#' Or the scale you want to get the IBS; and if it is a single point the return value will be the Brier Score at the timepoint.
#' @param parallel logical, whether to use parallel processing (default: TRUE)
#' @param n_cores integer, number of cores to use for parallel processing. If NULL, uses all available cores minus 1.
#' @param plot logical, whether to create a plot of Brier scores over time (default: FALSE)
#' @param return_brier logical, whether to return individual Brier scores along with IBS (default: FALSE)
#'
#' @return If return_brier = FALSE, returns the integration of brierscore (IBS). 
#' If return_brier = TRUE, returns a list containing:
#' \itemize{
#'   \item IBS: the integrated Brier score
#'   \item brier_scores: vector of individual Brier scores at each time point
#'   \item time_points: corresponding time points
#' }
#'
#' @author Hanpu Zhou \email{zhouhanpu@csu.edu.cn}
#' @references
#' HooraMoradian, DenisLarocque, & FranoisBellavance. (2017). \\(l_1\\) splitting rules in survival forests. Lifetime Data Analysis, 23(4), 671–691.
#'
#' Graf, Erika, Schmoor, Claudia, Sauerbrei, & Willi, et al. (1999). Assessment and comparison of prognostic classification schemes for survival data. Statist. Med., 18(1718), 2529-2545.
#'
#' Brier, G. W. . (1950). Verification of forecasts expressed in terms of probability. Monthly Weather Review, 78.
#'
#' Gneiting, T. , &  Raftery, A. E. . (2007). Strictly Proper Scoring Rules, Prediction, and Estimation.

#' @examples
#' library(survival)
#' library(SurvMetrics)
#' library(parallel)
#' set.seed(123)
#' N <- 100
#' mydata <- SDGM4(N, p = 20, c_step = -0.5)
#' index.train <- sample(1:N, 2 / 3 * N)
#' data.train <- mydata[index.train, ]
#' data.test <- mydata[-index.train, ]
#'
#' time_interest <- sort(data.train$time[data.train$status == 1])
#' sp_matrix <- matrix(sort(runif(nrow(data.test) * length(time_interest)),
#'   decreasing = TRUE
#' ), nrow = nrow(data.test))
#' object <- Surv(data.test$time, data.test$status)
#'
#' # the default time points with parallel processing and plotting
#' result <- IBS(object, sp_matrix, time_interest, parallel = TRUE, plot = TRUE, return_brier = TRUE)
#' # a time range with parallel processing and custom cores
#' IBS(object, sp_matrix, c(18:100), parallel = TRUE, n_cores = 4, plot = TRUE)
#'
#' @importFrom survival Surv
#' @importFrom parallel makeCluster stopCluster parLapply detectCores clusterEvalQ clusterExport
#' @importFrom graphics plot lines points legend title
#' @importFrom grDevices dev.new
#'
#'
#' @export
IBS <- function(object, sp_matrix, IBSrange = c(-2, -1), parallel = TRUE, n_cores = NULL, plot = FALSE, return_brier = FALSE) {
  # case1、coxph AND testing set
  if (inherits(object, "coxph")) {
    obj <- object
    test_data <- sp_matrix
    
    # the interesting times of training set
    distime <- sort(unique(as.vector(obj$y[obj$y[, 2] == 1])))
    
    mat_coxph <-
      pec::predictSurvProb(obj, test_data, distime) # get the survival probability matrix
    object_coxph <- Surv(test_data$time, test_data$status)
    
    object <- object_coxph
    sp_matrix <- mat_coxph
    if (max(IBSrange) <= 0) {
      IBSrange <- range(object[, 1])
    } # the fixed time point
  }
  
  
  # case2、RSF AND testing set
  if (inherits(object, c("rfsrc"))) {
    obj <- object
    test_data <- sp_matrix
    
    mat_rsf <-
      predict(obj, test_data)$survival # get the survival probability matrix
    object_rsf <- Surv(test_data$time, test_data$status)
    
    object <- object_rsf
    sp_matrix <- mat_rsf
    if (max(IBSrange) <= 0) {
      IBSrange <- range(object[, 1])
    } # the fixed time point
  }
  
  
  # case3 survreg AND testing set
  if (inherits(object, c("survreg"))) {
    obj <- object
    test_data <- sp_matrix
    
    # the interesting times of training set
    distime <- sort(unique(as.vector(obj$y[obj$y[, 2] == 1])))
    
    object <- Surv(test_data$time, test_data$status)
    sp_matrix <- pec::predictSurvProb.survreg(obj, test_data, distime)
    if (max(IBSrange) <= 0) {
      IBSrange <- range(object[, 1])
    } # the fixed time point
  }
  
  if (!is.numeric(IBSrange)) {
    stop("The class of the IBSrange must be numeric! or use the default setting")
  }
  
  if (any(is.na(IBSrange))) {
    stop("Cannot calculate IBS in the interval containing NA")
  }
  
  #default time range for IBS()
  if(max(IBSrange) <= 0){
    IBSrange = range(object[, 1])
  }
  
  if (!inherits(object, "Surv")) {
    stop("object is not of class Surv")
  }
  
  if (missing(object)) {
    stop("The survival object of the testing set is missing")
  }
  
  if (missing(sp_matrix)) {
    stop("The prediction of the survival probability matrix is missing")
  }
  
  if (any(is.na(object))) {
    stop("The input vector cannot have NA")
  }
  
  if (any(is.na(sp_matrix))) {
    stop("The input probability matrix cannot have NA")
  }
  
  if (length(object) != nrow(sp_matrix)) {
    stop(
      "number of rows of the sp_matrix and the survival object have different lengths"
    )
  }
  
  if (any(IBSrange <= 0)) {
    stop("The integration interval must be positive")
  }
  
  if (any(diff(IBSrange) <= 0)) {
    stop("The integral interval value must increase")
  }
  
  # Helper function for plotting Brier scores
  plot_brier_scores <- function(time_points, brier_scores, ibs_value) {
    # Create the plot
    plot(time_points, brier_scores, 
         type = "l", 
         col = "blue", 
         lwd = 2,
         xlab = "Time", 
         ylab = "Brier Score",
         main = paste("Brier Scores Over Time\n(IBS =", round(ibs_value, 4), ")"),
         cex.main = 1.2,
         cex.lab = 1.1)
    
    # Add points
    points(time_points, brier_scores, col = "blue", pch = 19, cex = 0.8)
    
    # Add horizontal line for mean Brier score (IBS approximation)
    mean_brier <- mean(brier_scores)
    abline(h = mean_brier, col = "red", lty = 2, lwd = 2)
    
    # Add grid for better readability
    grid(col = "lightgray", lty = 3)
    
    # Add legend
    legend("topright", 
           legend = c("Brier Score", "Mean Brier Score"),
           col = c("blue", "red"),
           lty = c(1, 2),
           lwd = c(2, 2),
           pch = c(19, NA),
           cex = 0.9,
           bg = "white")
  }
  # Helper function to return results
  format_results <- function(ibs_value, time_points = NULL, brier_scores = NULL, return_brier = FALSE, plot = FALSE) {
    if (plot) {
      plot_brier_scores(time_points, brier_scores, ibs_value)
    }
    
    if (return_brier) {
      return(list(
        IBS = ibs_value,
        brier_scores = brier_scores,
        time_points = time_points
      ))
    } else {
      return(ibs_value)
    }
  }
  compute_brier_parallel <- function(indices, object, sp_matrix, IBSrange) {
    sapply(indices, function(i) {
      pre_sp <- sp_matrix[, i]
      t_star <- IBSrange[i]
      Brier(object, pre_sp, t_star)
    })
  }
  
  # Determine the input dimension
  if (ncol(data.frame("sp" = sp_matrix)) == 1 &
      length(IBSrange) == 1) {
    bs <- Brier(object, sp_matrix, IBSrange)
    names(bs) <- "Brier Score"
    return(round(bs, 6))
  } else if ((ncol(data.frame("sp" = sp_matrix)) - 1) * (length(IBSrange) -
                                                         1) == 0) {
    stop("The length is illegal")
  } else if (ncol(sp_matrix) == length(IBSrange)) {
    IBSrange <- sort(IBSrange)
    
    # Parallel computation of Brier scores
    if (parallel && length(IBSrange) > 1) {
      # Set up parallel processing
      if (is.null(n_cores)) {
        n_cores <- max(1, parallel::detectCores() - 1)
      }
      
      # Only use parallel processing if we have multiple cores and multiple time points
      if (n_cores > 1 && length(IBSrange) > 2) {
        cl <- parallel::makeCluster(n_cores)
        
        # Export necessary objects and functions to cluster
        parallel::clusterExport(cl, c("object", "sp_matrix", "IBSrange", "Brier"), 
                                envir = environment())
        
        # Load required packages on each worker
        parallel::clusterEvalQ(cl, {
          library(survival)
          library(SurvMetrics)
        })
        
        # Split indices for parallel processing
        indices_list <- split(1:length(IBSrange), 
                              cut(1:length(IBSrange), n_cores, labels = FALSE))
        
        # Compute Brier scores in parallel
        t_brier_list <- parallel::parLapply(cl, indices_list, 
                                            compute_brier_parallel,
                                            object = object, 
                                            sp_matrix = sp_matrix,
                                            IBSrange = IBSrange)
        
        # Stop cluster
        parallel::stopCluster(cl)
        
        # Combine results
        t_brier <- unlist(t_brier_list)
      } else {
        # Fall back to sequential processing
        t_brier <- sapply(1:length(IBSrange), function(i) {
          pre_sp <- sp_matrix[, i]
          t_star <- IBSrange[i]
          Brier(object, pre_sp, t_star)
        })
      }
    } else {
      # Sequential computation
      t_brier <- sapply(1:length(IBSrange), function(i) {
        pre_sp <- sp_matrix[, i]
        t_star <- IBSrange[i]
        Brier(object, pre_sp, t_star)
      })
    }
    
    t_IBS <- 0
    for (i in 1:(length(IBSrange) - 1)) {
      t_IBS <- t_IBS + (IBSrange[i + 1] - IBSrange[i]) * t_brier[i]
    }
    t_IBS <- t_IBS / (range(IBSrange)[2] - range(IBSrange)[1])
    names(t_IBS) <- "IBS"
    ret <- format_results(t_IBS,time_points = IBSrange, brier_scores=t_brier, return_brier = return_brier, plot=plot)
    return(ret)
  } else {
    t_IBSrange <- range(IBSrange)
    p <- ncol(sp_matrix)
    # Simulation to generate discrete time series
    IBSrange <- seq(t_IBSrange[1], t_IBSrange[2], length = p)
    
    # Parallel computation of Brier scores
    if (parallel && p > 1) {
      # Set up parallel processing
      if (is.null(n_cores)) {
        n_cores <- max(1, parallel::detectCores() - 1)
      }
      
      # Only use parallel processing if we have multiple cores and multiple time points
      if (n_cores > 1 && p > 2) {
        cl <- parallel::makeCluster(n_cores)
        
        # Export necessary objects and functions to cluster
        parallel::clusterExport(cl, c("object", "sp_matrix", "IBSrange", "Brier"), 
                                envir = environment())
        
        # Load required packages on each worker
        parallel::clusterEvalQ(cl, {
          library(survival)
          library(SurvMetrics)
        })
        
        # Split indices for parallel processing
        indices_list <- split(1:p, cut(1:p, n_cores, labels = FALSE))
        
        # Compute Brier scores in parallel
        t_brier_list <- parallel::parLapply(cl, indices_list, 
                                            compute_brier_parallel,
                                            object = object, 
                                            sp_matrix = sp_matrix,
                                            IBSrange = IBSrange)
        
        # Stop cluster
        parallel::stopCluster(cl)
        
        # Combine results
        t_brier <- unlist(t_brier_list)
      } else {
        # Fall back to sequential processing
        t_brier <- sapply(1:p, function(i) {
          pre_sp <- sp_matrix[, i]
          t_star <- IBSrange[i]
          Brier(object, pre_sp, t_star)
        })
      }
    } else {
      # Sequential computation
      t_brier <- sapply(1:p, function(i) {
        pre_sp <- sp_matrix[, i]
        t_star <- IBSrange[i]
        Brier(object, pre_sp, t_star)
      })
    }
    
    t_IBS <- 0
    for (i in 1:(length(IBSrange) - 1)) {
      t_IBS <- t_IBS + (IBSrange[i + 1] - IBSrange[i]) * t_brier[i]
    }
    t_IBS <- t_IBS / (range(IBSrange)[2] - range(IBSrange)[1])
    names(t_IBS) <- "IBS"
    
    return(list(t_IBS, t_brier, IBSrange))
  }
}