library(reticulate)
use_virtualenv("r-keras", required = TRUE)
tf <- import("tensorflow")
library(keras3)
NLL2 <<- new_layer_class(
  classname = "NLL",
  initialize = function( chol, q,  G, sigma_e=0.5, ...) {
    super$initialize(...)
    # initialize sigma
    self$G <- G
    self$q <- q
    self$n_params <- self$q * (self$q + 1) / 2
    self$sigma <- sigma_e
    self$chol <- chol
    self$sigma_e <- self$add_weight(
      name="sigma_e",
      shape=shape(),
      initializer =tf$constant_initializer(sqrt(self$sigma)),
      trainable = F
    )
    self$chol <-  self$add_weight(
      shape=list(as.integer(c(self$n_params))),
      initializer = initializer_constant(self$chol), 
      trainable=F,
      name="Cholesky"
    )
  },
  get_vars = function(){
    return(
      list(
        sig2e=as.numeric(self$sigma_e),
        D = as.array(self$chol)))
  },
  # setting up the random effects matrix
  get_indices = function(N, Z_idx, min_Z) {
    cbind(0:(N-1), Z_idx - min_Z)
  },
  getZ = function(N, Z_idx, min_Z, max_Z) {
    # Use list() not c() for the shape parameter
    Z_idx_flat <- tf$reshape(Z_idx, shape = list(-1L))
    
    # Create indices
    row_indices <- tf$range(N, dtype = tf$int64)
    col_indices <- tf$cast(Z_idx_flat - min_Z, dtype = tf$int64)
    indices <- tf$stack(list(row_indices, col_indices), axis = 1L)
    
    # Create shape
    n_rows <- tf$cast(N, dtype = tf$int64)
    n_cols <- tf$cast(max_Z - min_Z + 1L, dtype = tf$int64)
    output_shape <- tf$stack(list(n_rows, n_cols))
    
    # Scatter
    Z0 <- tf$scatter_nd(
      indices = indices,
      updates = tf$ones(tf$cast(N, tf$int32), dtype = tf$float32),
      shape = output_shape
    )
    
    return(Z0)
  },
  kronecker_product_tf = function(W_small, n_subjects) {
    # Create identity operator
    I_n <- tf$eye((n_subjects))
    
    # Use TensorFlow's LinearOperator for efficient Kronecker
    W_op <- tf$linalg$LinearOperatorFullMatrix(W_small)
    I_op <- tf$linalg$LinearOperatorFullMatrix(I_n)
    
    # Kronecker product operator
    kron_op <- tf$linalg$LinearOperatorKronecker(list(W_op, I_op))
    
    # Convert to dense matrix
    W_full <- kron_op$to_dense()
    
    return(W_full)
  },
  single_k_loss_fn = function(y_true, y_pred, Z_idxs, posterior_probs) {
    # loop
    
    ll_g <- vector("list", self$G)
    for(ki in seq_len(self$G)){
      mu <- y_pred[, ki, drop = FALSE]; #tf$print(mu); tf$print(y_pred)
      residuals_flat <- tf$reshape(y_true - mu, list(-1L)); #tf$print("residuals_flat", residuals_flat)   # note: using mu, not y_pred
      
      Z_idx <- tf$reshape(Z_idxs[[1]], list(-1L)); #tf$print("Z_idx", Z_idx)
      t     <- tf$cast(tf$reshape(Z_idxs[[2]], list(-1L)), tf$float32);# tf$print("time:", t); tf$print("unmod_post", posterior_probs)
      w     <- tf$cast(tf$reshape(posterior_probs[,ki, drop=FALSE], list(-1L)), tf$float32); #tf$print("posteriors:", w)
      w_safe <- tf$maximum(tf$ones_like(w), tf$constant(1e-8, dtype = tf$float32)); #tf$print(w_safe)
      
      var <- tf$exp(self$sigma_e); #tf$print(var)
      
      # --- Reconstruct D (q x q) from Cholesky vector ---
      row_idx <- unlist(lapply(1:self$q, function(i) rep(i, i)))
      col_idx <- unlist(lapply(1:self$q, function(i) 1:i))
      tri_indices <- tf$constant(cbind(row_idx - 1L, col_idx - 1L), dtype = "int32")
      L_D <- tf$scatter_nd(tri_indices, self$chol,
                           shape = as.integer(c(self$q, self$q)))
      D <- tf$matmul(L_D, L_D, transpose_b = TRUE)
      self$est_D_full <- D
      
      # --- Group observations by subject (requires sorted rowids for RaggedTensor) ---
      min_Z <- tf$reduce_min(Z_idx); #tf$print(min_Z)
      subj0 <- tf$cast(Z_idx - min_Z, tf$int64); #tf$print("subj0",subj0)
      sort_order <- tf$argsort(subj0); #tf$print("ordered index", sort_order)
      subj0_s <- tf$gather(subj0, sort_order)
      t_s     <- tf$gather(t,     sort_order)
      w_s     <- tf$gather(w_safe, sort_order)
      res_s   <- tf$gather(residuals_flat, sort_order)
      
      t_padded   <- tf$RaggedTensor$from_value_rowids(t_s,   subj0_s)$to_tensor(default_value = 0.0)
      w_padded   <- tf$ones_like(t_padded)
      res_padded <- tf$RaggedTensor$from_value_rowids(res_s, subj0_s)$to_tensor(default_value = 0.0)
      n_i_vec    <- tf$RaggedTensor$from_value_rowids(t_s,   subj0_s)$row_lengths()  # (S,)
      
      S     <- tf$shape(t_padded)[1]; #tf$print("Size of training", S)
      n_max <- tf$shape(t_padded)[2]; #tf$print("n_max", n_max) # size of padding, not size of training set
      mask  <- tf$cast(tf$sequence_mask(n_i_vec, maxlen = n_max), tf$float32)  # (S, n_max)
      
      # --- Build T_padded: (S, n_max, q), columns are [1, t, t^2, ...] ---
      T_cols <- list(tf$ones_like(t_padded))
      if (self$q > 1) {
        for (k in 1:(self$q - 1)) T_cols[[k + 1]] <- t_padded ^ k
      }
      T_padded <- tf$stack(T_cols, axis = -1L)                                  # (S, n_max, q)
      
      # --- Per-subject V blocks in one batched op ---
      TD  <- tf$matmul(T_padded, D)                                             # (S, n_max, q)
      TDT <- tf$matmul(TD, T_padded, transpose_b = TRUE)                        # (S, n_max, n_max)
      # tf$print("Variance-Covariance Matrix of Batch (TDT):",TDT, summarize=200)
      # tf$print("Diagonal of Variance-Covariance Matrix of Batch (R_diag):",var / w_padded, summarize=200)
      # tf$print("Weights (w_padded):",w_padded, summarize=200)
      # tf$print("Variance (var):",var, summarize=200)
      R_diag <- var/w_padded                                                  # (S, n_max)
      V_batch <- TDT + tf$linalg$diag(R_diag) 
      # tf$print("Variance-Covariance Matrix of Batch:",V_batch, summarize=200)
      
      
      # Make padded rows/cols behave as identity so Cholesky is well-defined
      mask_outer <- tf$expand_dims(mask, 2L) * tf$expand_dims(mask, 1L)
      I_batch <- tf$eye(n_max, batch_shape = list(S))
      V_safe  <- V_batch * mask_outer + I_batch * (1.0 - mask_outer)
      #print(V_safe)
      # --- One batched Cholesky + triangular solve for ALL subjects ---
      L_batch <- tf$linalg$cholesky(V_safe)                                     # (S, n_max, n_max)
      alpha   <- tf$linalg$triangular_solve(
        L_batch, tf$expand_dims(res_padded, -1L), lower = TRUE
      )                                                              # (S, n_max, 1)
      # chk <- function(x, nm) cat(sprintf("%-12s nan=%s\n", nm, as.logical(tf$reduce_any(tf$math$is_nan(tf$cast(x, tf$float32))))))
      # chk(res_padded, "res_padded")
      alpha   <- alpha * tf$expand_dims(mask, -1L)                              # zero out padding
      # tf$print("Cholesky Factor of Batch:",L_batch, summarize=200)
      # tf$print("Residuals (Padded):",res_padded, summarize=200)
      
      # if singular, L_batch will have NaNs, and alpha will have NaNs. This will propagate to the loss.
      # tf$print("Alpha:",alpha, summarize=200)
      # find approximate inverse of V_batch using the Cholesky factor
      L_inv_batch <- tf$linalg$triangular_solve(L_batch, tf$eye(n_max, batch_shape = list(S)), lower = TRUE)  # (S, n_max, n_max)
      V_inv_batch <- tf$matmul(L_inv_batch, L_inv_batch, transpose_a = TRUE)  # (S, n_max, n_max)
      # tf$print("Inverse of Variance-Covariance Matrix of Batch:",V_inv_batch, summarize=200)
      
      quad_vec    <- tf$reduce_sum(alpha * alpha, axis = c(1L, 2L))             # (S,)
      diag_L      <- tf$linalg$diag_part(L_batch)                               # (S, n_max)
      log_det_vec <- 2.0 * tf$reduce_sum(tf$math$log(diag_L) * mask, axis = 1L) # (S,)
      n_i_float   <- tf$cast(n_i_vec, tf$float32)
      
      #   u        <- tf$unique(Z_idx)            # $y = sorted unique ids, $idx = dense 0-based labels
      # seg_dense <- u$idx                       # [0 0 0 ... 127 127 127], length N, contiguous
      # S         <- tf$size(u$y)                # actual number of subjects (128, not 250)
      
      # per_subject <- tf$math$segment_mean(weights, seg_dense)   # (S,) — no phantom rows
      # seg <- tf$cast(Z_idx, tf$int32) - 1L; tf$print("segments", seg)                 # 0-based: [0 0 0 ... 0 1 1]
      # subject 0 is the reindexed series of subject IDs used for the padding, just repurposed here to compute the segment mean, it's perfect!
      per_subject <- tf$math$segment_mean(w, subj0) ; #tf$print("weights", w); tf$print("weights for k1", per_subject);
      nll_per_subj <- 0.5 * (quad_vec + log_det_vec + n_i_float * 1.8378770664);# tf$print("nll_per_subject_k1", nll_per_subj);
      
      ll_g[[ki]] <- per_subject * nll_per_subj ; #tf$print("weighted_nll_per_subjk1", weighted_nll_per_subjk1)
      # 
    }
    
    ll <- op_stack(ll_g, axis = 1L)
    # return the average loss over the batch size, to ensure the gradients are computed properly
    loss <- tf$reduce_mean(tf$reduce_sum(tf$cast(ll, tf$float64), axis = 1L))
    return(loss)
  }, 
  call = function(y_true, y_pred, Z_idxs, posterior_probs) {
    loss <- self$single_k_loss_fn(y_true, y_pred, Z_idxs, posterior_probs)
    self$add_loss(loss)
    return(y_pred)
  }
)

# library(keras3)

# reticulate::py_config()                       
# reticulate::virtualenv_python("r-keras")     
# 
# tf <- reticulate::import("tensorflow")
# k  <- reticulate::import("keras")
# np <- reticulate::import("numpy")
# 
# c(tensorflow = tf$`__version__`,
#   keras      = k$`__version__`,
#   numpy      = np$`__version__`,
#   python     = reticulate::py_config()$version)
# 
# tf$config$list_physical_devices()           
# for reproducability
# pkgs <- reticulate::py_list_packages(envname = "r-keras")
# pkgs[, c("package", "version")]

# Save a pip-style requirements file
# writeLines(pkgs$requirement, "requirements-r-keras.txt")

# reticulate::virtualenv_create("r-keras", python = "3.10")
# reticulate::virtualenv_install("r-keras", requirements = "requirements-r-keras.txt")

