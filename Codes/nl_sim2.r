nl_sim2 <- function(seed, n, censor_rate, long_form_g=y~1+ measure_time + I(measure_time^2), ran_form=~1+ measure_time, k_tofit=3, estimate=TRUE) {
  #%%
  library("lcmm")
  library("NormPsy")
  # library(tictoc)
  paquidS <- paquid[which(paquid$agedem> paquid$age_init),]
  paquidS$normMMSE <-normMMSE(paquidS$MMSE)
  paquidS$age65 <- (paquidS$age - 65)/10
  # 
  # mj1 <- jlcmm(fixed = normMMSE ~ age65 +I(age65^2),
  #              random =~ age65+I(age65^2),
  #              survival = Surv(age_init, agedem,dem)~1, hazard = "Weibull",
  #              subject = 'ID', data = paquidS, verbose = F, logscale = T, maxiter = 1000)
  # summary(mj1)
  # 
  # 
  # tic()
  #  mj3b <- gridsearch(rep = 75, maxiter = 30, minit = mj1, cl = parallel::detectCores() - 1,
  #                    jlcmm(fixed = normMMSE ~ age65 + I(age65^2),
  #                          random = ~ age65 + I(age65^2),
  #                          mixture = ~ age65 + I(age65^2),
  #                          ng = 3,
  #                          survival = Surv(age_init, agedem, dem) ~ 1,
  #                          hazard = "Weibull",
  #                          subject = 'ID', data = paquidS,
  #                          B = mj1, logscale = T)
  #                    )
  # 
  # toc()
  # plot(mj3b, which = "fit", var.time = "age65")
  # plot(mj3b, which="survival", shades=T)
  # save(mj3b,file= "init3model.RData")
  load("/home/teddy/Documents/Codes/init4model.RData")
  summary(mj4b)
  vcov(mj4b)
  init <- estimates(mj4b, cholesky=F)
  D <- VarCovRE(mj4b)
  
  D <- init[which(grepl('varcov', names(init)))]
  
  cholesky <- metafor::vec2mat(D, diag=T)
  
  cholesky[3,1] <-cholesky[2,2] 
  cholesky[2,2] <- cholesky[1,3]
  cholesky[1,3] <- cholesky[3,1]
  
  corrs <- cov2cor(cholesky)
  
  
  idWeib_shape <- which(grepl('Weibull2', names(init)))
  idWeib_scale <- which(grepl('Weibull1', names(init)))
  
  
  
  corrs[3,2] <- 0.36
  corrs[2,3] <- 0.36
  
  corrs[1,2] <- -0.55
  corrs[2,1] <- -0.55
  corrs[1,3] <- 0.23
  corrs[3,1] <- 0.23
  covs <- round(MBESS::cor2cov(corrs, sd=c(4.862681/3, 6.644651/5, 2.145413/4)), 3)
  covs <- Matrix::nearPD(covs)
  covs <- as.matrix(covs$mat); print(covs)
  is.positive.definite(as.matrix(covs))
  cholesky <- t(chol(covs))
  # cov2cor(covs)
  (scales <- (exp(init[idWeib_scale]))^(-1/exp(init[idWeib_shape])))
  (shapes <- exp(init[idWeib_shape]))
  
  scales <- c(44, 50, 55, 60)
  shapes <- c(3, 6, 15, 25)
  
  # scales <- scales[-4]
  # shapes <- shapes[-4]
  #  matrix(D[,1], nrow=3, ncol=3, byrow = F)
  # cholesky <- metafor::vec2mat(D[,1], diag=T)
  # cholesky[3,1] <-cholesky[2,2] 
  # cholesky[2,2] <- cholesky[1,3]
  # cholesky[1,3] <- cholesky[3,1]
  # cholesky <- round(cholesky)
  # init[4:11] <- exp(init[4:11])
  # )
  
  idxs_slope <- which(grepl(c('intercept', "age65 ", "age65)")[2], names(init)))
  
  b1 <- init[(idxs_slope-3)]
  b2 <- init[idxs_slope]
  b3 <- init[idxs_slope+3]
  
  df <- df2 <- (data.frame(int=b1,mt=b2, mt2=b3))
  betas_matrix <- as.matrix(df2)
  # betas_matrix <- betas_matrix[-4,]
  residual_err <- 2
  betas_matrix
  
  # Load required libraries
  library(dplyr)
  library(ggplot2)
  #%%
  # Set parameters for simulation
  # set.seed(123)  # for reproducibility
  # seed <- 1;n <- 500; 
  censor_rate <- censor_rate
  # seed <- 2
  set.seed(seed)
  # source("NNEM/helper_functions/util.R")
  
  # example 3: Xue & Yao 2022
  n_samples <- n 
  x <- runif(n_samples, -5, 12)
  # X <- cbind(1, x) 
  x <- sort(x)
  true_mixing_prop <- function(x) {
    pi1 <- (2 * exp(-0.1 * x^4)) / (1 + exp(-0.1 * x^4)) * as.numeric(x < 3)
    pi2 <- (1 - (2 * exp(-0.1 * x^4)) / (1 + exp(-0.1 * x^4))) * as.numeric(x < 0)
    pi3 <- (2 * exp(-0.1 * (x - 6)^4)) / (1 + exp(-0.1 * (x - 6)^4)) * as.numeric(x >= 3)
    pi4 <- 1 - (pi1 + pi2+ pi3)
    return(cbind(pi1, pi2, pi3, pi4))
  }
  library(dplyr)
  true_pi <- pi_init <- true_mixing_prop(x)
  
  df <- data.frame(x,pi_init)
  df<-df|>arrange(x)
  plot(df$x,df$pi1, type="l", col="blue")
  lines(df$x, df$pi2)
  lines(df$x,df$pi3, col="green")
  lines(df$x,df$pi4, col="red")
  component_assignments <- apply(pi_init, 1, function(p) sample(1:k_tofit, 1, prob = p))
  true_C <- component_assignments
  # component_assignments <- sample(1:3,size=n, prob=c(.3,.4, .3), replace = T)
  # probs <- matrix(rep(c(.3,.4, .3), n), ncol = 3, byrow = T)
  # pi_init <- probs
  shapes <- round(shapes)
  scales <- round(scales)
  #  scales <- scales^(-1/shapes)
  # parms <- get_parms("weibulli")
  # Generate data
  times <- numeric(n)
  # i <- 1
  # lam_D <- parms$lam_D
  time_L <- runif(n, min=0, max=1)
  Nsub <- n
  m <- length(shapes)
  g <- components <- component_assignments; 
  table(components)/n
  # slopes <- parms$slopes
  library(flexsurv)
  # print(-1*slopes)
  #   g <- components
  # ebx <- exp(rowSums(slopes[g,] * cbind(X3,X4, X5))) # proportional hazards term
  for (i in 1:m) {
    component_indices <- which(components == i)
    if (length(component_indices) > 0) {
      # Calculate linear predictor
      shape_g <- shapes[i]
      scale_g <- scales[i]
      # adjusted_scales <- scales[i] * exp(linear_pred[component_indices])
      #     times[component_indices] <- gen_model_survival(ebx[component_indices],
      # dist = "weibulli", 
      # scale=scale_g, 
      # shape=shape_g)
      times[component_indices] <- rweibull(length(component_indices), shape=shape_g, scale=scale_g)
    }
  }
  # hist(times)
  # scales <- 
  # censor_rate <- "0.2"; #!!NOTE censor_rate=0.2 is more like 90% and the if you don't multiply by 10 and it's about 50% if you multply lambda by 10; censor_rate=0.5 yields an actual censor rate of 70% if you multiply by 10 and an actual censor rate of 0.25% if you don't multiply 
  # r <- 0.05 # 5%
  # # r <- 0.2 # 25%
  # # r <- 0.5 # 50%
  # tmp_parms <- parms
  time_L <- numeric(Nsub)
  
  if (censor_rate==0){ 
    time_C <- Inf
  } else if(censor_rate==0.05){
    time_C <- rweibull(n, shape =5, scale =100)
  }else if(censor_rate==0.25){
    time_C <- rweibull(n, shape = 3, scale = 125)
  }else if(censor_rate==0.5){
    time_C <- rweibull(n, shape = 2, scale = 100)
  }
  delta <- as.numeric(time_C  >= times); print(table(delta)/n);
  time_Y <- pmin(times, time_C)
  # hist(time_Y)
  # hist(time_C)
  
  # Model parameters
  n_subjects <- n     # number of individuals
  # n_timepoints <-     # number of measure_time points per individual
  # 
  # hist(table(paquidS$ID))
  # median(table(paquidS$ID))
  # Generate the data structure (Inspired by Pinheiro & Bates, 2000, Ch. 6)
  data <- data.frame()
  tot <- c(0,6,12,18,24,30,36,42,48,54,60) # include a baseline recording at 0 
  for(i in seq_len(n_subjects)){
    # i <- 1
    
    measure_times <- tot
    
    data_i <- expand.grid(
      ID = i,
      # measure_time = seq(from=0,to=(n_timepoints), by=6)
      # measure_time =runif()
      
      measure_time = measure_times
      
    )
    data <- rbind(data,data_i)}
  
  
  # Generate individual random effects (random intercepts)
  # random_effects <- data.frame(
  #   ID = 1:n_subjects,
  #   u_i = rnorm(n_subjects, mean = 0, sd = sigma_u)
  # )
  
  datalong <- data.frame()
  ctr <- 0
  for(i in 1:n){
    # i<-1
    data_tmp1 <- data[data$ID==i,]
    data_tmp <- data_tmp1[(data_tmp1$measure_time)<time_Y[i],]; if(nrow(data_tmp1)!=nrow(data_tmp)){ctr <- ctr+1}
    data_tmp$event_time <- time_Y[i]
    data_tmp$delta <- delta[i]
    datalong <- rbind(datalong, data_tmp)
  }
  data <- datalong
  sd_e <- residual_err
  datalong$measure_time <- datalong$measure_time/10
  library(mvtnorm)
  mu_g <- matrix(0, 2, 2)
  # rho <- 0.5
  # cov <- rho*sqrt(9)*sqrt(4)
  cholesky
  # D_g_og <- diag(2)
  # prop_0 <- 1
  # prop_1 <- 0.5
  # prop_2 <- 1.5
  # prop_3 <- 2.5
  #? for the proportions, one must be very careful with the ordering particularly for initialization otherwise you will run into issues with estimation
  # props <- c(prop_3, prop_2, prop_0)
  # D_matrix2 <- (prop_2^2)*D_g_og
  # D_matrix3 <- (prop_3^2)*D_g_og
  #  i <- 1
  
  
  
  
  datalong$measure_time2 <- datalong$measure_time^2
  X <- model.matrix(~1+ measure_time +(measure_time2), data=datalong)
  Z <- X
  data_new <- time_Y
  
  # for(i in 1:m){
  # component_indices <- which(components == i)
  # if (length(component_indices) > 0) {
  #     D <- D_g_og
  #     # ranef <- pi_init[,i]*
  #     ranef <- pi_init[,i]*rmvnorm(n, sigma = D)
  #     bi <- bi+ranef
  # # init_ranef <- ranef + init_ranef
  # }
  # }
  
  
  # 
  
  # # fixef <- c(0.5,2,1)
  # betas_matrix <- matrix(c(0.5,2, 
  #                         3, 2,
  #                         2, 1), ncol = 2, byrow=T)
  # # betas_matrix <- matrix(c(fixef,fixef), nrow = 3, ncol=2)
  pseudo_g <- g[datalong$ID]
  # # ranefsx <- bi[[pseudo_g]]
  # ranefs <- ranef[datalong$ID,]
  
  # Zb <- ranefs[,1] + ranefs[,2]*datalong$measure_time
  # y <- betas_matrix[,1][pseudo_g]+betas_matrix[,2][pseudo_g]*datalong$measure_time+ Zb + rnorm(nrow(datalong), sd=sd_e) 
  library(MASS)
  library(matrixcalc)
  n_timepoints <- table(data$ID)
  # Create design matrix
  # Z <- model.matrix(~measure_time, data = data)
  # X <- Z
  sigma_e <- residual_err
  D <-covs
  V_blocks <- list()
  
  # nonlinear_mean_function <- function(X_i, g){
  #   if(g==1){
  #   return())
  # } else if(g==2){
  #   return(-3 + 3 *cos(X_i[,2]-2))
  # } else {
  #   return(10 + 2*cos(0.5*X_i[,2]-4))
  # }
  # }
  
  #%%
  nonlinear_mean_function <- function(X_j, g){
    if(g==1){
      return(30*((exp(-0.3*X_j[,2])) -exp(-0.8502*X_j[,2])))
    } else if(g==2){
      return(5*((exp(-0.1*X_j[,2])) -exp(-1*X_j[,2])))
    } else if(g==3){
      return(50*((exp(-0.2*X_j[,2])) -exp(-0.5*X_j[,2])))
    } else {
      return(15*((exp(-0.5*X_j[,2])) -exp(-2*X_j[,2])))
    }}
  
  xb <- c()
  for(i in seq_len(k_tofit)){
    XB <- nonlinear_mean_function(X, i)
    xb <- cbind(xb, XB)
  }
  
  trueXB <- XB_true <- xb
  
  # omega <- props
  sigma_e <- 2
  Z <- Z
  
  # D <- matrix(c(27,4,4,1.2), nrow = 2, byrow = T)
  # chol(D)
  # cov2cor(D)
  # save(D, file="covmat-sim1.RData")
  # load("covmat-sim1.RData")
  D <- covs
  mu <- c()
  for (i in 1:n) {
    # i <- 1
    idx <- which(data$ID == i)
    g <- component_assignments[i]
    Z_i <- Z[idx, , drop = FALSE]
    V_i <- Z_i %*% (D) %*% t(Z_i) + (sigma_e^2) * diag(length(idx))
    X_i <- X[idx, ,drop=FALSE]
    XBg <- nonlinear_mean_function(X_i, g)
    V_blocks[[i]] <- V_i
    mu <- c(mu, XBg)
  }
  V <- Matrix::bdiag(V_blocks)
  library(mvnfast)
  y <- t(rmvn(1, mu, V))
  datalong$y <- y
  head(datalong,10)
  
  head(mu,10)
  # Clean up intermediate columns (optional): don't need to do this if we generate from a poison
  # data <- data %>%
  #   select(ID, measure_time, y) %>%
  #   arrange(ID, measure_time)
  
  # Display first few rows
  head(datalong, 15)
  datalong$g <-pseudo_g 
  datalong$mu <- mu
  # Summary statistics
  
  cat("Summary of generated data:\n")
  cat("Number of subjects:", n_subjects, "\n")
  cat("Number of observations per subject:", mean(as.numeric(n_timepoints)), "\n")
  cat("Total observations:", nrow(data), "\n")
  cat("Mean response:", round(mean(datalong$y), 2), "\n")
  cat("SD of response:", round(sd(datalong$y), 2), "\n")
  
  # Plot a sample of trajectories
  sample_subjects <- sample(1:n_subjects, min(50, n_subjects))
  data_sample <- datalong 
  
  # p <- ggplot(data_sample, aes(x = measure_time, y =  y, group = ID, color = factor(g))) +
  #   geom_line() +
  #   labs(title = "Sample Individual Trajectories",
  #        x = "Time",
  #        y = "Response (Y)",
  #        color = "Subject ID") +
  #   theme_minimal()
  # print(p)
  
  
  
  pal4 <- c("1" = "#2a78d6",   
            "2" = "#eb6834",  
            "3" = "#1baf7a",
            "4" = "#6250d6")
  
  p <- ggplot(data_sample, aes(measure_time, y, group = ID)) +
    geom_line(colour = "grey75", alpha = 0.5, linewidth = 0.35) +
    geom_line(data = data_sample, aes(measure_time, mu, colour = factor(g)),
              inherit.aes = FALSE, linewidth = 1.1) +
    facet_grid(~ g, labeller = label_both) +
    scale_colour_manual(values = pal4, guide = "none") +
    labs(x = "Time", y = "Response (Y)") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  
  print(p)
  
  #%%
  # Verify the model by fitting it back to the generated data
  # library(lme4)
  # model_fit <- lmer(y ~ measure_time + (measure_time+measure_time2| ID), data = datalong)
  # print(summary(model_fit))
  
  
  data <- datalong
  
  m <- nrow(betas_matrix)
  D_g <- D
  # alphas <- props
  # alphas
  # gammas_matrix <- slopes[1:m,]
  # mu_g <- mu_g
  nis<-table(data$ID)
  
  
  library('data.table')
  sim <- data.table(data)
  temp2 <- subset(sim, select = c(ID, event_time, delta)) # baseline
  data_wide <- unique(temp2)
  data_wide$x <- x
  X <- model.matrix(~1+ measure_time +(measure_time2), data=datalong)
  n <- nrow(data_wide)
  Z <- X
  
  
  
  library(nnet)
  cholesky <- t(chol(D))
  tcrossprod(cholesky)
  # nis <- table(num_measure)
  #? compute class specific random effects 
  # k_tofit <- 4
  # compute initial estimates for a G=1 model
  # D <- VarCovRE(mj4b)
  
  theta_orig <- list(pj=pi_init, scale=scales, shape=shapes, 
                     betas=betas_matrix,  
                     cholesky=cholesky, D=D, sigma_e=sd_e)
  save(theta_orig,file = "theta_orig_nnjlcm_genk3_p2.RData")  # expects theta_orig
  
  
  if(estimate){
    
    mj1 <- hlme(fixed = y~1+ measure_time + I(measure_time^2),
                #  mixture=~1,
                # classmb = ~ X1 + X2,
                random = ~1+ measure_time + I(measure_time^2),
                subject='ID',
                ng=1, data=data)
    # summary(mj1)
    
    pi_init <- matrix(1/k_tofit, nrow = length(x), ncol = k_tofit)
    vcovmj1 <- Matrix::nearPD(vcov(mj1))$mat
    init <- estimates(mj1, cholesky=F)
    # sample from the asymptotic distribution
    sample_init <- mvnfast::rmvn(k_tofit, init, sigma = vcovmj1)
    colnames(sample_init) <- names(init)
    
    D_est <- sample_init[1,which(grepl('varcov', colnames(sample_init)))]
    
    cholesky <- metafor::vec2mat(D_est, diag=T)
    
    cholesky[3,1] <-cholesky[2,2] 
    cholesky[2,2] <- cholesky[1,3]
    cholesky[1,3] <- cholesky[3,1]
    
    # idWeib_shape <- which(grepl('Weibull2', colnames(sample_init)))
    # idWeib_scale <- which(grepl('Weibull1', colnames(sample_init)))
    
    # (scales <- (exp(sample_init[,idWeib_scale]))^(-1/exp(sample_init[,idWeib_shape])))
    # (shapes <- exp(sample_init[,idWeib_shape]))
    # scales_est <- init[c(10,8,6,4)] 
    
    # scales_est <- (init[,1]^2)^(-1/(init[,2]^2))
    # shapes_est <- exp(sample_init[,2])
    
    # gammas_matrix <- sample_init[1,(max(idWeib_shape)+1):(max(idWeib_shape)+2)]
    
    # betas_matrix <- sample_init[,2:8]
    
    idxs_slope <- which(grepl(c('intercept', "measure_time", "measure_time2")[2], colnames(sample_init)))
    
    b1 <- sample_init[,(idxs_slope[1]-1)]
    b2 <- sample_init[,idxs_slope[1]]
    b3 <- sample_init[(idxs_slope[2])]
    
    df <- df2 <- (data.frame(int=b1,mt=b2, mt2=b3))
    # df[2,] <- df2[1,]
    # df[1,] <- df2[3,]
    # df[3,] <- df2[4,]
    # df[4,] <- df[2,]
    betas_matrix <- as.matrix(df)
    # betas_matrix
    # 
    library(fastDummies)
    true_resp <- dummy_cols(component_assignments)[,-1]
    var.time <- "measure_time"
    break.times <- NULL
    # obtain the measure times and compute break times using quantiles of the measure times
    timeInterv <- datalong[,var.time]
    if(is.null(break.times)){
      break.times <- quantile(timeInterv,prob=seq(0,1,length.out=10))
    }else if(length(break.times)==1){
      break.times <- quantile(timeInterv,prob=seq(0,1,length.out=break.times))
    }
    data <- datalong
    
    ng <- k_tofit
    # nea <- sum(x$idea)
    
    
    ntps <- length(break.times)-1 # number of time intervals for which predictions will be computed
    # C <- apply(P, 1, which.max)
    C <- component_assignments
    # table(C)
    
    pprob <- data.frame(ID=data_wide$ID,true_resp)
    colnames(pprob) <- c("ID", paste0("prob",1:ng))
    
    numsg <- lapply(1:ng,function(g) data_wide$ID[which(C==g)])
    times <- data[,var.time]
    maxT <- sapply(numsg,function(y) max(times[which(data[,"ID"] %in% y)]))
    
    pred_complete <- data.frame(ID=datalong$ID, y=datalong$y, XB_true)
    colnames(pred_complete) <- c("ID", "y", paste0("pred_m",1:ng))
    pred0 <- pred_complete #ID,obs,pred_m1...ng i.e. longitudinal recordings and marginal predictions for each class
    
    ppi <- pprob #ID,prob1...ng i.e. posterior probabilities of class membership
    name.id <- colnames(data)[1]
    colnames(ppi)[1] <- name.id
    
    pred0 <- cbind(pred0,times)
    pred_complet0 <- merge(pred0,ppi,by=name.id)
    
    # pred_complet0
    
    pred_times0 <- matrix(NA,length(break.times)-1,4*ng+1)
    pred_times0[,1] <- (break.times[-length(break.times)] + break.times[-1])/2
    colnames(pred_times0) <- c(var.time,paste("pred.class",1:ng,sep=""),paste("obs.class",1:ng,sep=""),
                               paste("lower.obs.class",1:ng,sep=""),paste("upper.obs.class",1:ng,sep=""))
    
    for(i in 1:ntps){
      indg <- lapply(1:ng,function(g) which( (pred_complet0[,ng+3]>=break.times[i]) & (pred_complet0[,ng+3]<min(break.times[i+1],maxT[g])) ) ) #which sets of predictions correspond to times in the interval [break.times[i],break.times[i+1]) for class g
      if(i==(length(break.times)-1)) indg <- lapply(1:ng, function(g){ if(maxT[g]>break.times[i]) c(indg[[g]],which( (pred_complet0[,ng+3]==min(break.times[i+1],maxT[g]))  )) else indg[[g]]})
      
      if(any(sapply(indg,length)>0)){
        ##predictions ponderees
        pred_times0[i,1+1:ng] <- sapply(1:ng,function(g) sum(pred_complet0[indg[[g]],2+g] * pred_complet0[indg[[g]],ng+3+g])/sum(pred_complet0[indg[[g]],ng+3+g]))
        
        ##observations ponderees
        pred_times0[i,1+ng+1:ng] <- sapply(1:ng,function(g) sum(pred_complet0[indg[[g]],2] * pred_complet0[indg[[g]],ng+3+g])/sum(pred_complet0[indg[[g]],ng+3+g]))
        
        ##calcul des ecarts-types des moyennes ponderes
        sd_obs <- sapply(1:ng,function(g) sd(pred_complet0[indg[[g]],2])*sqrt(sum(pred_complet0[indg[[g]],ng+3+g]**2)/sum(pred_complet0[indg[[g]],ng+3+g])**2)  )
        ##ic des obs moyennes
        pred_times0[i,1+2*ng+1:ng] <- pred_times0[i,1+ng+1:ng] - 1.96 * sd_obs
        pred_times0[i,1+3*ng+1:ng] <- pred_times0[i,1+ng+1:ng] + 1.96 * sd_obs
      }
      
    }
    
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
    
    p <- ggplot(pred_long, aes(x = measure_time, colour = class, fill = class)) +
      geom_ribbon(aes(ymin = lower.obs, ymax = upper.obs),
                  alpha = 0.2, colour = NA) +
      geom_point(aes(y = obs)) +
      geom_line(aes(y = pred)) +
      facet_grid(~class)+
      labs(x = "Measurement time", y = "Longitudinal responses", title = "Weighted Marginal Predictions vs Weighted Observed Trajectories", subtitle = "Weighted by posterior probabilities of class membership",
           colour = "Class", fill = "Class",
           caption = "Solid = predicted, dashed = observed, ribbon = obs 95% CI") +
      theme_minimal()
    
    print(p)
    # [1,] 68.2015 64.5285 75.9308 59.3352
    # [2,] 24.2364 57.8511 22.0628 48.6829
    # [3,] 13.5360 56.0213 23.8152 39.4700
    # [4,] 68.2015 64.5285 75.9308 59.3352
    # [5,] 59.5253 64.8401 30.0949 64.5957
    # [6,] 10.0633 55.5094 25.2694 36.9508
    # [7,] 68.2015 64.5285 75.9308 59.3352
    # [8,] 69.1451 65.5881 58.0926 61.1216
    # [9,] 66.7514 65.8853 42.8476 66.2642
    # [10,] 65.7692 65.9348 39.0409 67.1376
    # 
    # intercept class1 intercept class2 intercept class3 intercept class4
    # 1         67.69441          59.4171          65.2603          77.0778
    # 2         25.57216          47.1815          60.5106          22.6542
    # 3         12.67933          39.9419          58.2702          25.2003
    # 4         67.69441          59.4171          65.2603          77.0778
    # 5         60.02183          63.9490          65.9170          30.0107
    # 6          8.82468          37.7300          57.5897          26.2225
    # > 
    
    sd_e <- residual_err <-abs(sample_init[1,"stderr"])
    #  matrix(D[,1], nrow=3, ncol=3, byrow = F)
    # print(pi_init)
    
    
    library(eha)
    # this weibull model from the EHA package is parameterized EXACTLY in terms of my Weibull model. 
    prg <- phreg(Surv(event_time, delta) ~ 1, data = data_wide, dist = "weibull")
    scale_shape <- mvnfast::rmvn(k_tofit, prg$coefficients, sigma = prg$var)
    
    (scales <- exp(scale_shape[,1]))
    (shapes <- exp(scale_shape[,2]))
    
    # library(lme4)
    # lme_model <- lmer(y ~ measure_time + I(measure_time^2) + (measure_time  | ID), data = datalong)
    # merDeriv::vcov.lmerMod(lme_model, full=TRUE)
    # summary(lme_model)
    
    # lme_model
    cholesky <- as.matrix(t(chol(Matrix::nearPD(cholesky)$mat)))
    
    
    theta <- list(pj=pi_init, scale=scales, shape=shapes, 
                  betas=betas_matrix,  
                  cholesky=cholesky, D=D, sigma_e=sd_e)
    print(theta)
    XB <- X %*% t(betas_matrix)
    # XB1 <- nnet(y~X[,2], linout=T, size=50, weights = theta$pj_init[datalong$ID,1], maxit=2000, decay=0.1)
    # XB2 <- nnet(y~X[,2], linout=T, size=50, weights = theta$pj_init[datalong$ID,2], maxit=2000, decay=0.1)
    # XB3 <- nnet(y~X[,2], linout=T, size=50, weights = theta$pj_init[datalong$ID,3], maxit=2000, decay=0.1)
    # XB4 <- nnet(y~X[,2], linout=T, size=50, weights = theta$pj_init[datalong$ID,4], maxit=2000, decay=0.1)
    
    # XB1 <- fitted(XB1)
    # XB2 <- fitted(XB2)
    # XB3 <- fitted(XB3)
    # XB4 <- fitted(XB4)
    
    # XB <- cbind(XB1, XB2, XB3, XB4)
    
    
    return(list(theta=theta, datalong=datalong, XB=XB, pi_init=pi_init, data_wide=data_wide, true=trueXB, truepi=true_pi, true_C=component_assignments, theta_orig=theta_orig))}
  else{
    return(list(theta=NULL, datalong=datalong, XB=NULL, pi_init=pi_init, data_wide=data_wide, true=trueXB, truepi=true_pi, true_C=component_assignments, theta_orig=theta_orig))
  }
}


train <- nl_sim2(seed=3, n=250, censor_rate = 0.05, k_tofit=4)
test <- nl_sim2(seed=4, n=1000,censor_rate = 0.05, k_tofit=4)
# Run the main model

result <- nnem(k=4, long_form_g =  y~1+ measure_time + I(measure_time^2), ranform =~1+ measure_time + I(measure_time^2), data_long = train$datalong, data_wide = train$data_wide, id_var = "ID", gform=~x, y_var = "y", event_var = "delta", measure_var = "measure_time", event_time_var = "event_time", learning_rate_lme = 0.001, maxit =20)

result <- nnem_nnlme(k = 4, long_form_g = y~1+ measure_time, ran_form =~1+ measure_time + I(measure_time^2), g_form = ~x, subject_id = "ID", datalong = train$datalong, data_wide=train$data_wide, XB_init = train$XB, theta = train$theta, print_surv = FALSE, max.iter = 30, prog_plot = F, NWG = FALSE, nproc=-1, plotSurvScores = FALSE, optimizer = "mla", learning_rate = 0.001, batch_size = nrow(train$datalong), weight_decay_moe = 0.01)

run_simulation_nnlme(train, test, 3, 1000, getwd(), 4)