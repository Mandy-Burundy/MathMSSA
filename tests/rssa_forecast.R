#!/usr/bin/env Rscript

# Stage A reference: Rssa 1.1, kind=mssa, svd.method=svd,
# rforecast(direction="row", base="original", len=1).
# Package installation stays outside this script. The Rssa version must be exactly 1.1.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  message("usage: rssa_forecast.R fixtures.json rssa_output.json rawdir")
  quit(status = 1L)
}

fixtures_path <- args[[1L]]
output_path <- args[[2L]]
raw_dir <- args[[3L]]

suppressPackageStartupMessages(library(Rssa))
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  message("jsonlite is not installed")
  quit(status = 1L)
}
suppressPackageStartupMessages(library(jsonlite))

rssa_version <- as.character(packageVersion("Rssa"))
if (!identical(rssa_version, "1.1")) {
  message(sprintf("Rssa version %s, required 1.1", rssa_version))
  quit(status = 2L)
}

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

write_f64 <- function(path, values) {
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.double(values), con, endian = "little")
}

row_major <- function(mat) as.double(as.vector(t(mat)))

embed_hankel <- function(series, L) {
  series <- as.numeric(series)
  k <- length(series) - L + 1L
  idx <- outer(seq_len(L), seq_len(k), function(i, j) i + j - 1L)
  matrix(series[idx], nrow = L, ncol = k)
}

as_json_number <- function(x) {
  if (length(x) != 1L || is.na(x)) return(NULL)
  if (is.infinite(x)) return("Infinity")
  as.numeric(x)
}

forecast_rank <- function(s, r, V, K) {
  fc <- rforecast(
    s,
    groups = list(seq_len(r)),
    len = 1L,
    direction = "row",
    base = "original",
    only.new = TRUE,
    drop = FALSE
  )
  forecast <- as.numeric(fc[[1L]][[1L]])
  Vr <- V[, seq_len(r), drop = FALSE]
  ends <- cumsum(K)
  W <- Vr[ends, , drop = FALSE]
  Q <- Vr[-ends, , drop = FALSE]
  A <- diag(length(K)) - tcrossprod(W)
  Cmat <- tcrossprod(W, Q)
  ev <- tryCatch(svd(A, nu = 0L, nv = 0L)$d, error = function(e) NA_real_)
  cond <- if (any(!is.finite(ev)) || min(ev) == 0) Inf else max(ev) / min(ev)
  fit <- tryCatch(qr.coef(qr(A), Cmat), error = function(e) NULL)
  residual <- if (is.null(fit)) {
    NA_real_
  } else {
    numerator <- norm(A %*% fit - Cmat, type = "2")
    denominator <- norm(A, type = "2") * norm(fit, type = "2") + norm(Cmat, type = "2")
    if (denominator == 0) 0 else numerator / denominator
  }
  rc <- reconstruct(s, groups = list(seq_len(r)))
  list(
    forecast = forecast,
    condition = cond,
    residual = residual,
    recon = lapply(rc[[1L]], as.numeric),
    A = A,
    C = Cmat
  )
}

run_point <- function(point) {
  id <- point$id
  channels <- lapply(point$channels, as.numeric)
  L <- as.integer(point$l)
  ranks <- as.integer(unlist(point$ranks))
  lens <- vapply(channels, length, 1L)
  K <- lens - L + 1L
  neig <- min(max(ranks) + 1L, L, sum(K))
  fitted <- tryCatch(
    ssa(channels, L = L, kind = "mssa", svd.method = "svd", neig = neig),
    error = function(e) e
  )
  if (inherits(fitted, "error")) {
    return(list(id = id, ok = FALSE, message = conditionMessage(fitted)))
  }
  ns <- asNamespace("Rssa")
  tm <- do.call(cbind, lapply(channels, embed_hankel, L = L))
  stored <- as.matrix(get(".get.or.create.trajmat.mssa", ns)(fitted))
  extmat_diff <- max(abs(stored - tm))
  ff <- get(".F", ns)(fitted)
  z <- unlist(lapply(seq_along(ff), function(idx) {
    ff[[idx]][seq(to = fitted$length[[idx]], length.out = K[[idx]] - 1L)]
  }), use.names = FALSE)
  sigma <- as.numeric(fitted$sigma)
  sigma <- sigma[seq_len(min(neig, length(sigma)))]
  U <- fitted$U[, seq_len(min(neig, ncol(fitted$U))), drop = FALSE]
  V <- calc.v(fitted, seq_len(ncol(U)))
  write_f64(file.path(raw_dir, paste0(id, ".traj.f64")), row_major(tm))
  write_f64(file.path(raw_dir, paste0(id, ".z.f64")), z)
  write_f64(file.path(raw_dir, paste0(id, ".U.f64")), row_major(U))
  write_f64(file.path(raw_dir, paste0(id, ".V.f64")), row_major(V))
  ranks_out <- list()
  for (r in ranks) {
    piece <- tryCatch(forecast_rank(fitted, r, V, K), error = function(e) e)
    if (inherits(piece, "error")) {
      ranks_out[[as.character(r)]] <- list(ok = FALSE, message = conditionMessage(piece))
      next
    }
    write_f64(file.path(raw_dir, sprintf("%s.r%d.A.f64", id, r)), row_major(piece$A))
    write_f64(file.path(raw_dir, sprintf("%s.r%d.C.f64", id, r)), row_major(piece$C))
    ranks_out[[as.character(r)]] <- list(
      ok = TRUE,
      forecast = piece$forecast,
      condition = as_json_number(piece$condition),
      residual = as_json_number(piece$residual),
      recon = piece$recon,
      aDim = dim(piece$A),
      cDim = dim(piece$C)
    )
  }
  list(
    id = id,
    ok = TRUE,
    message = NULL,
    lengths = as.integer(lens),
    blockWidths = as.integer(K),
    endpoints = as.integer(cumsum(K)),
    sigma = sigma,
    frobenius = sqrt(sum(tm * tm)),
    extmatDiff = extmat_diff,
    dims = list(traj = dim(tm), z = length(z), U = dim(U), V = dim(V)),
    ranks = ranks_out
  )
}

fixtures <- fromJSON(fixtures_path, simplifyVector = FALSE)
points <- lapply(fixtures$points, run_point)
payload <- list(
  rssa = rssa_version,
  rVersion = paste(R.version$major, R.version$minor, sep = "."),
  fixture = fixtures_path,
  points = points
)
writeLines(
  toJSON(payload, digits = 17L, auto_unbox = TRUE, pretty = TRUE, null = "null"),
  output_path,
  useBytes = TRUE
)
quit(status = 0L)
