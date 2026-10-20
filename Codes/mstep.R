#! M - step
mstep <- function(postprob, theta, data_wide, data_long, ranform, long_form_g,  gform, EM_iter, learning_rate_lme, inheritedModels, batch_size){

  n <- nrow(data_wide)
  m <- ncol(postprob)
  x <- model.frame(gform, data_wide)
# mixing proportions using nnet

    # postprob <- e_step_result$postprob
# CM Step I
    nn_model <- nnet(
      x = x,          # Input features
      y = postprob,          # Target is component probabilities
      size = 10,      # Hidden layer size
      decay = 0.01,   #! Weight decay: this is something which can be tweaked
      maxit = 5000,   # Max iterations
      trace = FALSE,  # No output during training
      softmax = TRUE,  # Use softmax for multinomial
      entropy=TRUE # maximum conditional likelihood fitting
    )
    
    # Predict probabilities (this code doesn't work for some stupid reason, 
  # extract directly from the nnet object)
    # prob <- predict(nn_model, data_wide$X3)
    new_prior <- fitted(nn_model)


print("New priors")
print(head(new_prior))
  
random_effects <- theta$ranef

init <- c(theta$cholesky[lower.tri(theta$cholesky, diag = T)],log(theta$sigma_e)) 

k_tofit <- ncol(theta$pj)
# IDs <- datalong$ID
z_cols <- c('z0', 't')
# max(IDs)
gammas <- postprob
gammas <- as.data.frame(gammas)
IDs <- as.character(data_wide$ID)          # attach the key as a column
rownames(gammas) <- IDs                     # now rownames are clean "1","2",...
IDs <- data_long$ID
# 2. Index by CHARACTER, not integer, so it matches by name not position:
IDs <- as.character(data_long$ID)                           # c("1","2","3","5",...)
gammas <- gammas[IDs, ]   
gammas <- as.matrix(gammas)
# gammas <- as.numeric(gammas)
storage.mode(gammas) <- "double"
IDs <- as.numeric(IDs)


# Initialize number of variance components
X <- model.matrix(long_form_g, data = data_long)
Z <- model.matrix(long_form_g, data = data_long)
y <- model.frame(long_form_g, data=data_long)
y <- model.response(y)

# Create input layers
Z_input <- layer_input(shape = c(1), dtype = "int64", name = "Z_input")
t_input <- layer_input(shape = c(1), name = "t_input")

# List of Z inputs
Z_inputs <- list(Z_input, t_input)

# colnames(data) <- c('z0', "t", "event_time", "delta", "y", "g")
# X_train <- matrix(data$t)
# x_cols <- setdiff(colnames(X_train), c('z0'))
# (input_dim <- dim(X_train)[2])
# z_train <- data[c(z_cols)]
# input <- layer_input(shape = 1L)
X_input <- layer_input(shape=c(1L), name="X_input")
y_true_input <- layer_input(shape=c(1L), name="y_true_input")
# gammas_inputs_k1 <- layer_input(shape = (1L), name="gammas_k1")
# gammas_inputs_k2 <- layer_input(shape=(1L), name="gammas_k2")
# x_cols <- setdiff(colnames(X), c('z0'))
trainedModel <- list()
G <- as.integer(ncol(gammas))
print(G)
chols <- theta$cholesky
sigma_e <- theta$sigma_e
k_tofit <- ncol(postprob)
print("EM_iter"); print(EM_iter)
if(EM_iter==1){
  
  # this uses the functional api
  hidden_out <- X_input %>%
    layer_dense(units = 256, name = "net1_dense2") %>%
    layer_activation_leaky_relu() 
  
  
  y_pred_output <- layer_dense(hidden_out, units = k_tofit, name = "net1_output")
  
  # linear activation
  # New input for the per-observation weights (posteriors)
  w_input <- layer_input(shape = shape(k_tofit), name = "w_input", dtype = "float32")
  nll_layer <- NLL2(q = ncol(chols), G=G, chol = as.numeric(chols[lower.tri(chols, diag = T)]), sigma_e = as.matrix(sigma_e^2))
  nll <- nll_layer(y_true_input, y_pred_output, Z_inputs, w_input)
  model <- keras_model(
    inputs = c(list(X_input, y_true_input), Z_inputs, list(w_input)),
    outputs = nll
  )
  
  model |> 
    compile(
      optimizer = optimizer_lion(learning_rate = learning_rate_lme), 
      jit_compile = FALSE
    )
  # print(learning_rate_lme)

  train_data <- list(
    X_input = as.matrix(X[,2], ncol=1),                          # 1st
    y_true_input = matrix(y, ncol = 1),  # 2nd
    Z_input = matrix(IDs, ncol = 1),    # 3rd
    t_input = matrix(X[,2], ncol = 1),
    w_input = as.matrix((gammas[,1:k_tofit]))
  )
  
  
  # print(1)
  es <- callback_early_stopping(
    monitor = "loss",          
    patience = 20,
    min_delta = 1e-6,
    restore_best_weights = TRUE
  )
  history <- model |> fit(
    x = train_data,
    epochs = 250, #an arbitrarily large amount of epochs to run. 
    batch_size = batch_size,
    verbose = 0,
    shuffle = FALSE,
    callbacks = list(es)
  )
  
  
  y_pred <- predict(model, train_data, batch_size = batch_size, verbose = 0)
  # cat("Initializing keras1 mean-response at k =",iter, ":",  "\n")
  print(head(y_pred,30))
  
  trainedModel <- list(XB=y_pred,  nnlme=model,  nn_weights = lapply(model$get_weights(), as.array))
} else {
  
  
  train_data <- list(
    X_input = as.matrix(X[,2], ncol=1),                          # 1st
    y_true_input = matrix(y, ncol = 1),  # 2nd
    Z_input = matrix(IDs, ncol = 1),    # 3rd
    t_input = matrix(X[,2], ncol = 1),
    w_input = as.matrix((gammas[,1:k_tofit]))
  )
  
  # print(inheritedModels)
  
  model <- inheritedModels$nnlme; print(model)
  y_pred <- predict(model, train_data, batch_size = batch_size, verbose = 0); print(y_pred)
  layers <- sapply(model$layers, \(l) l$name)
  
  nll <- model$layers[[length(layers)]]; #print(nll)
  nll$sigma_e
  nll$chol
  #   tensorflow is weird. assign dynamically updates the weight without physically assigning using an operator i.e. <- or =.
  len_chol <- length(chols[lower.tri(chols, diag = T)]);# print(len_chol)
  nll$chol$assign(init[1:len_chol]); #print(init[1:len_chol])
  nll$sigma_e$assign(exp(init[(len_chol+1)])^2); #print((exp(init[(len_chol+1)])^2));
  # print(is.nan(gammas))
  # model
  es <- callback_early_stopping(
    monitor = "loss",           
    patience = 20,
    min_delta = 1e-6,
    restore_best_weights = TRUE
  )
  history <- model |> fit(
    x = train_data,
    epochs = 125, #an arbitrarily large amount of epochs to run. 
    batch_size = batch_size,
    verbose = 0,
    shuffle = FALSE,
    callbacks = list(es)
  )
  
  
  
  y_pred <- predict(model, train_data, batch_size = batch_size, verbose = 0)
  
  #  model |> 
  #   compile(
  #     optimizer = optimizer_sgd(learning_rate = 0.01, momentum = 0.9),
  
  trainedModel <- list(XB=y_pred,nnlme=model,  nn_weights = lapply(model$get_weights(), as.array))
  
}


fitted_nnlme_matrix <-trainedModel$XB

# marginal likelihood of the longitudinal submodel function to be fed into mla
nll_marg <- function(parms,XB_in, data, posterior, long_form_g, ranform, theta_true){
  
l2 <<- sum(lower.tri(theta_true$cholesky, diag = T))
log_marginal_density_yis <- function(data_long, XB_in, theta, long_form_g,ranform, m=ncol(posterior)) {
  # parms prep
  sigmas <- theta$sigma_e
  # D_matrices <- theta$D

  # data prep
  y <- model.frame(long_form_g,data_long)[,1]
  Xg <- model.matrix(long_form_g, data_long)
  Zg <- model.matrix(ranform, data_long)
  data_long$y <- y

  N <- length(unique(data_long$ID));#print(N) # number of unique subjects under study.
  m <- nrow(betas_matrix);#print(m);#print(betas_matrix) # number of unique groups for mixture model
  f_yi <- matrix(0,nrow = N, ncol = m)
  nis<-table(data_long$ID)
  IDs<-unique(data_long$ID)
  # length(nis)

  D_matrices <- tcrossprod(theta$cholesky);# print(D_matrices)
  for (k in 1:m) {
    sigma_g <- sigmas # error variances 
    D_g <- D_matrices # random effect variance-covariance matrices
    
    for(i in seq_len(length(IDs))){

      # i<-2
      # n_i <- nis[i]
      # k <- 1
      ids <- IDs[i] 
      
      # data_long <- data;i <- 1
      subset <- data_long[data_long$ID == ids, ] #! ensure that IDs are consistently named across datasets
            # subset<-sim[sim$ID==i,]
      Z_ig <- Zg[data_long$ID == ids, ,drop = FALSE]    # Random Effects Design Matrix
      X_ig <- Xg[data_long$ID==ids,,drop=FALSE ]# Cluster-specific Fixed Effects Design Matrix

      p <- ncol(X_ig)
      Xi_beta_k <- matrix(XB_in[data_long$ID==ids,k,drop=FALSE]) # Fixed Effects Linear Predictor
   
      q <- ncol(Z_ig)
      Y_i <- subset$y # Longitudinal recordings

      #? calculating the distributional parameters for the group specific random effects
      # D_g_inv <- solve(D_g)
      # 3-by-2 dimensions
      V_ig <- as.matrix(Z_ig %*% (D_g) %*% t(Z_ig) + sigma_g^2 * diag(nrow(X_ig)))
      # simgas[i] <- list(Sigma_i)
      f_yi[i, k] <- mvnfast::dmvn(t(Y_i), Xi_beta_k, V_ig, log = T)
    }
    }
  return(f_yi)
}


  init <- parms
  cholesky <- metafor::vec2mat(init[(1):(l2)], diag = T)
  # cholesky[3,1] <- cholesky[2,2]
  # cholesky[2,2] <- cholesky[1,3]

  cholesky[1,2] <- cholesky[1,3] <- cholesky[2,3]<- 0
  #print(tcrossprod(cholesky))
  sd_e <- exp(init[length(init)])
  #print(sd_e)
    tmp_theta <- list( sigma_e=sd_e,cholesky=cholesky)
ysu <- log_marginal_density_yis(
  data_long = data, 
  XB_in = XB_in,
theta=tmp_theta,
long_form_g = long_form_g,
ranform = ranform)
  
  llk <- sum(rowSums(posterior*(ysu)))
  return(-llk)
}


library(parallel)
library(marqLevAlg)
n_cores <- 1
opt_result <- marqLevAlg::mla(init, fn=nll_marg, XB_in=fitted_nnlme_matrix,data=data_long, posterior=postprob,nproc=n_cores, print.info=T, long_form_g=long_form_g, ranform=ranform, theta_true=theta)

cat("============================================")
cat("\nNN-JLCM Longitudinal Parameter Estimates:\n")
cat("============================================")
  sd_info_fixef <- summary(opt_result, loglik=T,, digits = 10)

  # l1 <- nrow(theta$betas)*ncol(theta$betas)
l2 <- sum(lower.tri(theta$cholesky, diag = T))
# obj <- opt_result
# (betas_matrix <- matrix(opt_result$b[1:l1], nrow = m, ncol=3))
  cholesky <- metafor::vec2mat(opt_result$b[(1):(l2)], diag=T)
  cholesky[1,2] <- cholesky[1,3] <- cholesky[2,3]<- 0
   (D <- tcrossprod(cholesky))
  (sd_e <- exp(opt_result$b[length(opt_result$b)]))

# theta <- theta_orig
shapes <- as.matrix(theta$shape)
scales <- as.matrix(theta$scale)

D_matrices <- D
# cholesky <- cholesky
sigmas <- sd_e
pi_init <- new_prior

parms <- c(log(shapes), log(scales), theta$gamma_matrix)
neg_loglik_wb <- function(parms_j, data_wide, posterior, ng){
  # parms_j <- parms
  log_shapes <- matrix(parms_j[1:m],ncol = 1); log_shapes # first m parameters are the shapes
  log_scales <- matrix(parms_j[(m+1):(2*m)]); log_scales # next m parameters are the scales
  gammas_matrix <- as.matrix(parms_j[(2*m+1):(length(parms_j))])

  f_td <-  function(shape, scale, t, d) {
    (d*(dweibull(t,shape=shape, scale = scale,log = T))) + (1-d)*(pweibull(t, shape = shape, scale =scale,lower.tail = F,log.p = T)) 
  }

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
      delta <- delta
      # fixed effects in survival component
      # wxi <- as.matrix(w)
      ti <- data_wide$event_time
      for(i in 1:N){
      f_tds[i, j] <- f_td(shape, scale, ti[i], delta[i])
      }
      }
    return(f_tds)
  }
  scales <- exp(log_scales)
  shapes <- exp(log_shapes)
  lwb <-  density_weibull(shape_params = shapes, scale_params = scales, data_wide=data_wide)
  component_likelihood <- sum(rowSums(posterior * (lwb)))
  return(-component_likelihood)
}
neg_loglik_wb(parms, data_wide = data_wide, posterior = postprob, ng=m)
opt_result_wb <- tryCatch({
  marqLevAlg::mla(b=parms, fn=neg_loglik_wb, data_wide=data_wide, posterior=postprob, ng=4, nproc=n_cores, print.info=T)}, error = function(e) {
    optim(par=parms, fn=neg_loglik_wb, data_wide=data_wide, posterior=postprob)})

cat("==========================================")
cat("\nNN-JLCM Survival Parameter Estimates:\n")
cat("==========================================")
sd_res_Td <- summary(opt_result_wb, loglik=T)
# sd_res_Td[,1:3] <- as.numeric(unlist(sd_res_Td[,1:3]))
# sd_res_Td[1:(2*m),c(1,5:6)] <- exp(sd_res_Td[1:(2*m),c(1,5:6)])
shapes <- exp(matrix(opt_result_wb$b[1:m],ncol = 1)) # first m bameters are the shapes
# sd_res_Td
scales <- exp(matrix(opt_result_wb$b[(m+1):(2*m)])) # next m bameters are the scales

   shapes_sd <- NULL
   scales_sd <- NULL
#  compute the random effects
  bi <- list()
  Xg <- model.matrix(long_form_g, data=data_long)
Zg <-  model.matrix(ranform, data=data_long)
yhat_mat <- matrix(0, nrow=nrow(data_long), ncol=ncol(gammas))
for(j in 1:m){
# j <- 1
    sigma_g <- sd_e # error variances 

    D_g <- D
     # random effect variance-covariance matrices
    # beta_g <- betas_matrix[j,] # fixed effects

  # i<-1
  bigs<-matrix(nrow = nrow(data_wide), ncol = ncol(D))
  seq <- data_wide$ID
  sim <- data_long
  yhat <- c()
    for(i in seq_len(length(seq))){

      ids <- seq[i]
      subset<-sim[sim$ID==ids,]
      Z_ig <- Z[sim$ID==ids, , drop = FALSE]
      # X_ig <- X[sim$ID==ids, ,drop=FALSE]
      Xi_beta <- fitted_nnlme_matrix[sim$ID==ids, j ,drop = FALSE]

      # p<-ncol(X_ig)
      # Xi_beta <- X_ig%*%beta_g 
      # xbs[i]<-list(Xi_beta)
      q<-ncol(Z_ig)
      Y_i <- subset$y
      # ys<-list(Y_i)
      sigma<-sigma_g
      n_i <- length(subset$ID)
      # 3by2 %*% 
      Vi <- as.matrix(Z_ig%*% D_g%*% t(Z_ig) + sigma_g^2 * diag(n_i))
      big <- D_g %*% t(Z_ig) %*% matrixcalc::matrix.inverse(Vi) %*% (Y_i - Xi_beta)
      bigs[i,] <- big
      
      yhats <- Xi_beta + Z_ig %*% big
      yhat <- rbind(yhat, yhats)
      
    }
  yhat_mat[,j] <- yhat
  bi[j] <- list(bigs)
}


# standard_errors <- list()
# update parameters
 theta <- list(pj=new_prior, scale=scales, XB=fitted_nnlme_matrix,shape=shapes, betas=betas_matrix, ranef=bi, D=D, cholesky=cholesky, sigma_e=sd_e, trainedModel=trainedModel, yhat=yhat_mat)
  results_standard_errors <- list("Longitudinal Submodel"=sd_info_fixef, "Survival Submodel"=sd_res_Td)
 ret <- list("Parameters"=theta, "Model Estimates & S.E."=results_standard_errors)
return(ret)
}

