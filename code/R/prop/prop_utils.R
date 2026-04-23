prop_symmetrize <- function(mat) {
  0.5 * (mat + t(mat))
}

prop_as_matrix <- function(x) {
  if (is.null(dim(x))) {
    return(matrix(x, ncol = 1))
  }
  x
}

prop_extract_xt <- function(x, t_index) {
  xt <- x[, t_index, , drop = FALSE]
  dim(xt) <- c(dim(x)[1], dim(x)[3])
  xt
}

prop_expand_scalar <- function(value, size, name) {
  if (length(value) == 1L) {
    return(rep(as.numeric(value), size))
  }
  if (length(value) != size) {
    stop(sprintf("%s must have length 1 or %d", name, size), call. = FALSE)
  }
  as.numeric(value)
}

prop_stable_chol <- function(mat, base_jitter = 1e-8, max_tries = 8L) {
  mat <- prop_symmetrize(mat)
  n <- nrow(mat)
  if (n == 0L) {
    return(matrix(0, 0, 0))
  }

  for (ii in seq_len(max_tries)) {
    jitter <- base_jitter * (10 ^ (ii - 1L))
    chol_try <- tryCatch(
      chol(mat + diag(jitter, n)),
      error = function(e) NULL
    )
    if (!is.null(chol_try)) {
      return(chol_try)
    }
  }

  eig <- eigen(mat, symmetric = TRUE)
  vals <- pmax(eig$values, base_jitter)
  eig$vectors %*% diag(sqrt(vals), nrow = length(vals))
}

prop_mvnorm_draw <- function(mean_vec, cov_mat) {
  chol_cov <- prop_stable_chol(cov_mat)
  as.vector(mean_vec + t(chol_cov) %*% rnorm(length(mean_vec)))
}

prop_build_design_stack <- function(x) {
  nt <- dim(x)[2]
  do.call(rbind, lapply(seq_len(nt), function(tt) prop_extract_xt(x, tt)))
}
