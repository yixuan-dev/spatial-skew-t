prop_symmetrize <- function(mat) {
  0.5 * (mat + t(mat))
}

prop_as_matrix <- function(x) {
  if (is.null(dim(x))) {
    return(matrix(x, ncol = 1))
  }
  x
}

prop_check_scalar_numeric <- function(value, name,
                                      lower = -Inf, upper = Inf,
                                      strict_lower = FALSE,
                                      strict_upper = FALSE,
                                      integerish = FALSE) {
  if (length(value) != 1L || !is.numeric(value) || is.na(value) || !is.finite(value)) {
    stop(sprintf("%s must be a single finite numeric value.", name), call. = FALSE)
  }

  value_num <- as.numeric(value)
  if (integerish && !isTRUE(all.equal(value_num, round(value_num), tolerance = 1e-8))) {
    stop(sprintf("%s must be an integer-valued scalar.", name), call. = FALSE)
  }
  if ((strict_lower && value_num <= lower) || (!strict_lower && value_num < lower)) {
    comparator <- if (strict_lower) ">" else ">="
    stop(sprintf("%s must be %s %s.", name, comparator, format(lower)), call. = FALSE)
  }
  if ((strict_upper && value_num >= upper) || (!strict_upper && value_num > upper)) {
    comparator <- if (strict_upper) "<" else "<="
    stop(sprintf("%s must be %s %s.", name, comparator, format(upper)), call. = FALSE)
  }

  if (integerish) {
    return(as.integer(round(value_num)))
  }
  value_num
}

prop_check_flag <- function(value, name) {
  if (length(value) != 1L || !is.logical(value) || is.na(value)) {
    stop(sprintf("%s must be a single TRUE/FALSE value.", name), call. = FALSE)
  }
  as.logical(value)
}

prop_check_numeric_vector <- function(x, name, length_expected = NULL) {
  if (is.null(x) || !is.atomic(x)) {
    stop(sprintf("%s must be a numeric vector.", name), call. = FALSE)
  }
  x_num <- as.numeric(x)
  if (length(x_num) < 1L || anyNA(x_num) || any(!is.finite(x_num))) {
    stop(sprintf("%s must contain only finite numeric values.", name), call. = FALSE)
  }
  if (!is.null(length_expected) && length(x_num) != length_expected) {
    stop(sprintf("%s must have length %d.", name, length_expected), call. = FALSE)
  }
  x_num
}

prop_check_numeric_matrix <- function(x, name,
                                      nrow_expected = NULL,
                                      ncol_expected = NULL,
                                      positive = FALSE,
                                      nonnegative = FALSE) {
  x_mat <- as.matrix(x)
  if (length(dim(x_mat)) != 2L) {
    stop(sprintf("%s must be a two-dimensional matrix-like object.", name), call. = FALSE)
  }
  if (!is.numeric(x_mat) || anyNA(x_mat) || any(!is.finite(x_mat))) {
    stop(sprintf("%s must contain only finite numeric values.", name), call. = FALSE)
  }
  if (!is.null(nrow_expected) && nrow(x_mat) != nrow_expected) {
    stop(sprintf("%s must have %d rows.", name, nrow_expected), call. = FALSE)
  }
  if (!is.null(ncol_expected) && ncol(x_mat) != ncol_expected) {
    stop(sprintf("%s must have %d columns.", name, ncol_expected), call. = FALSE)
  }
  if (positive && any(x_mat <= 0)) {
    stop(sprintf("%s must contain only positive values.", name), call. = FALSE)
  }
  if (nonnegative && any(x_mat < 0)) {
    stop(sprintf("%s must contain only nonnegative values.", name), call. = FALSE)
  }
  x_mat
}

prop_check_numeric_array3 <- function(x, name,
                                      dim1_expected = NULL,
                                      dim2_expected = NULL,
                                      dim3_expected = NULL) {
  dims <- dim(x)
  if (is.null(dims) || length(dims) != 3L) {
    stop(sprintf("%s must be a three-dimensional numeric array.", name), call. = FALSE)
  }
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("%s must contain only finite numeric values.", name), call. = FALSE)
  }
  if (!is.null(dim1_expected) && dims[1] != dim1_expected) {
    stop(sprintf("%s must have first dimension %d.", name, dim1_expected), call. = FALSE)
  }
  if (!is.null(dim2_expected) && dims[2] != dim2_expected) {
    stop(sprintf("%s must have second dimension %d.", name, dim2_expected), call. = FALSE)
  }
  if (!is.null(dim3_expected) && dims[3] != dim3_expected) {
    stop(sprintf("%s must have third dimension %d.", name, dim3_expected), call. = FALSE)
  }
  x
}

prop_require_names <- function(x, name, required) {
  if (!is.list(x)) {
    stop(sprintf("%s must be a list.", name), call. = FALSE)
  }
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(
      sprintf("%s is missing required entries: %s.", name, paste(missing, collapse = ", ")),
      call. = FALSE
    )
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
