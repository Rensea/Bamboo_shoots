# ============================================================
# R2稳定提升_六性状_防泄漏自动搜索版.R
# 目的：在不使用测试集信息、不靠单次SPXY挑高的前提下，尽量提高六个性状的稳定R²。
# 重点改动：
#   1) 所有光谱预处理都在每一次训练/测试划分内部完成，避免表二全样本转换造成信息泄漏。
#   2) 增加 stage 分层划分、重复验证、screen -> refine 两阶段自动搜索。
#   3) 同时输出 mean/median/q25/p95 R²、RMSE、gap，并汇总六性状平均R²。
#   4) 保留 SPXY 作为“诊断上限”，论文建议仍以重复随机/分层稳定性为主。
# ============================================================

rm(list = ls())
gc()

# ------------------------- 用户需要改这里 -------------------------
file_path <- "C:/Users/微笑久保/Desktop/222222.xlsx"
out_dir   <- "C:/Users/微笑久保/Desktop/竹笋代码"

# 速度与稳定性的折中：先少量重复筛选，再对前若干名做大量重复验证
screen_repeat <- 10       # 初筛重复次数；想更快可设 5；论文前可提高到 30
final_repeat  <- 100      # 复筛重复次数；论文前可提高到 200 或 500
screen_keep_n <- 8        # 每个性状初筛保留前多少个组合进入复筛
seed <- 2026

# stable_split 推荐：stage_stratified；如果想和原脚本完全一致，可改为 random_stratified
stable_split <- "stage_stratified"
train_ratio <- 0.80

# 分性状稳定R²目标：只考核可溶性糖、木质素、蛋白质、水分。
# 纤维素、硬度设置为 NA：保留建模结果，但不纳入达标统计。
# 达标判定默认采用：median_R2 >= goal 且 mean_R2 >= goal - 0.05
goal_map <- c(
  soluble_sugar = 0.75,
  lignin        = 0.75,
  protein       = 0.65,
  moisture      = 0.65,
  cellulose     = NA_real_,
  hardness      = NA_real_
)
q25_buffer <- 0.15   # q25_R2 >= goal - q25_buffer 时认为下四分位没有明显失稳
# ------------------------------------------------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ensure_packages <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, dependencies = TRUE)
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}
ensure_packages(c("readxl", "openxlsx", "dplyr", "pls", "e1071"))

# ============================================================
# 1. 基础函数
# ============================================================
DEFAULT_WAVELENGTHS <- c(
  401.82, 404.42, 405.76, 407.12, 408.48, 409.88, 411.3, 412.72,
  414.18, 415.66, 417.16, 418.68, 420.24, 421.8, 423.4, 425.04,
  426.68, 428.36, 430.08, 431.8, 433.58, 435.38, 437.2, 439.06,
  440.96, 442.9, 444.86, 446.86, 448.9, 451.0, 453.12, 455.28,
  457.48, 459.74, 462.04, 464.38, 466.78, 469.22, 471.72, 474.26,
  476.88, 479.54, 482.26, 485.06, 487.9, 490.82, 493.82, 496.88,
  500.0, 503.22, 506.5, 509.88, 513.34, 516.88, 520.54, 524.28,
  528.12, 532.05, 536.12, 540.28, 544.58, 549.0, 553.54, 558.22,
  563.05, 568.04, 573.19, 578.48, 583.98, 589.64, 595.52, 601.6,
  607.9, 614.44, 621.22, 628.28, 635.62, 643.26, 651.22, 659.54,
  668.2, 677.28, 686.76, 696.7, 707.12, 718.04, 729.52, 741.56,
  754.22, 767.54, 781.56, 796.32, 811.91, 828.36, 845.8, 864.28,
  883.91, 904.7, 926.64, 949.64, 974.02
)

get_wavelengths <- function(n_band) {
  if (n_band == length(DEFAULT_WAVELENGTHS)) DEFAULT_WAVELENGTHS else seq_len(n_band)
}

r2_score <- function(y, pred) {
  y <- as.numeric(y); pred <- as.numeric(pred)
  den <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  if (!is.finite(den) || den == 0) return(NA_real_)
  1 - sum((y - pred)^2, na.rm = TRUE) / den
}
rmse_score <- function(y, pred) sqrt(mean((as.numeric(y) - as.numeric(pred))^2, na.rm = TRUE))
mae_score  <- function(y, pred) mean(abs(as.numeric(y) - as.numeric(pred)), na.rm = TRUE)

safe_cor <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (length(unique(x[is.finite(x)])) <= 1 || length(unique(y[is.finite(y)])) <= 1) return(0)
  out <- suppressWarnings(cor(x, y, use = "complete.obs"))
  ifelse(is.na(out), 0, out)
}
safe_div <- function(a, b, eps = 1e-8) as.numeric(a) / (as.numeric(b) + eps)
ndi <- function(a, b, eps = 1e-8) (as.numeric(a) - as.numeric(b)) / (as.numeric(a) + as.numeric(b) + eps)

extract_stage <- function(id) {
  id <- as.character(id)
  out <- ifelse(grepl("^c[0-9]+$", id, ignore.case = TRUE), "C", sub("[-_].*$", "", id))
  toupper(out)
}

load_bamboo_data <- function(file_path) {
  dat0 <- readxl::read_excel(file_path)
  dat0 <- as.data.frame(dat0)
  if (!"ID" %in% names(dat0)) names(dat0)[1] <- "ID"
  dat0$ID <- as.character(dat0$ID)
  dat0$stage <- extract_stage(dat0$ID)

  spec_cols <- grep("^band[0-9]+_inner$", names(dat0), value = TRUE)
  if (length(spec_cols) == 0) spec_cols <- grep("^band", names(dat0), value = TRUE)
  if (length(spec_cols) == 0) stop("没有识别到 band 光谱列。")

  X <- as.matrix(dat0[, spec_cols, drop = FALSE])
  storage.mode(X) <- "numeric"
  colnames(X) <- paste0("band_", seq_len(ncol(X)))
  list(dat = dat0, X = X, spec_cols = spec_cols, wavelengths = get_wavelengths(ncol(X)))
}

# ============================================================
# 2. 训练集内安全预处理：避免表二全样本转换泄漏
# ============================================================
snv_transform <- function(X) {
  if (nrow(X) == 0) return(X)
  out <- t(apply(X, 1, function(z) {
    s <- sd(z, na.rm = TRUE)
    if (is.na(s) || s == 0) rep(0, length(z)) else (z - mean(z, na.rm = TRUE)) / s
  }))
  colnames(out) <- colnames(X)
  out
}
row_minmax_transform <- function(X) {
  if (nrow(X) == 0) return(X)
  out <- t(apply(X, 1, function(z) {
    rz <- range(z, na.rm = TRUE)
    if (!is.finite(diff(rz)) || diff(rz) == 0) rep(0, length(z)) else (z - rz[1]) / diff(rz)
  }))
  colnames(out) <- colnames(X)
  out
}
row_l2norm_transform <- function(X) {
  denom <- sqrt(rowSums(X^2, na.rm = TRUE)); denom[!is.finite(denom) | denom == 0] <- 1
  out <- X / denom
  colnames(out) <- colnames(X)
  out
}
area_norm_transform <- function(X) {
  denom <- rowSums(abs(X), na.rm = TRUE); denom[!is.finite(denom) | denom == 0] <- 1
  out <- X / denom
  colnames(out) <- colnames(X)
  out
}
smooth_ma <- function(X, k = 7) {
  Xs <- X; half <- floor(k / 2)
  for (i in seq_len(nrow(X))) {
    z <- as.numeric(X[i, ]); zz <- z
    for (j in seq_along(z)) {
      lo <- max(1, j - half); hi <- min(length(z), j + half)
      zz[j] <- mean(z[lo:hi], na.rm = TRUE)
    }
    Xs[i, ] <- zz
  }
  colnames(Xs) <- colnames(X)
  Xs
}
make_derivative <- function(X, order = 1, k = 7) {
  Xs <- smooth_ma(X, k = k)
  if (order == 1) {
    D <- Xs[, -1, drop = FALSE] - Xs[, -ncol(Xs), drop = FALSE]
  } else {
    D1 <- Xs[, -1, drop = FALSE] - Xs[, -ncol(Xs), drop = FALSE]
    D <- D1[, -1, drop = FALSE] - D1[, -ncol(D1), drop = FALSE]
  }
  colnames(D) <- paste0("d", order, "_", seq_len(ncol(D)))
  D
}
detrend_linear_transform <- function(X) {
  if (nrow(X) == 0) return(X)
  idx <- seq_len(ncol(X))
  out <- t(apply(X, 1, function(z) residuals(lm(z ~ idx))))
  colnames(out) <- colnames(X)
  out
}
msc_fit_apply <- function(Xtr, Xte) {
  ref <- colMeans(Xtr, na.rm = TRUE)
  correct_one <- function(z) {
    fit <- lm(as.numeric(z) ~ ref)
    b0 <- coef(fit)[1]; b1 <- coef(fit)[2]
    if (!is.finite(b1) || abs(b1) < 1e-12) b1 <- 1
    (as.numeric(z) - b0) / b1
  }
  Xtr2 <- t(apply(Xtr, 1, correct_one))
  Xte2 <- if (nrow(Xte) == 0) Xte else t(apply(Xte, 1, correct_one))
  colnames(Xtr2) <- colnames(Xtr); colnames(Xte2) <- colnames(Xte)
  list(train = Xtr2, test = Xte2)
}
col_zscore_fit_apply <- function(Xtr, Xte) {
  mu <- colMeans(Xtr, na.rm = TRUE)
  s <- apply(Xtr, 2, sd, na.rm = TRUE); s[!is.finite(s) | s == 0] <- 1
  list(train = scale(Xtr, center = mu, scale = s), test = scale(Xte, center = mu, scale = s))
}
col_center_fit_apply <- function(Xtr, Xte) {
  mu <- colMeans(Xtr, na.rm = TRUE)
  list(train = scale(Xtr, center = mu, scale = FALSE), test = scale(Xte, center = mu, scale = FALSE))
}
col_minmax_fit_apply <- function(Xtr, Xte) {
  mn <- apply(Xtr, 2, min, na.rm = TRUE); mx <- apply(Xtr, 2, max, na.rm = TRUE)
  den <- mx - mn; den[!is.finite(den) | den == 0] <- 1
  train <- sweep(sweep(Xtr, 2, mn, "-"), 2, den, "/")
  test  <- sweep(sweep(Xte, 2, mn, "-"), 2, den, "/")
  list(train = train, test = test)
}

apply_preprocess_pair <- function(Xtr, Xte, method) {
  method <- as.character(method)
  if (method == "raw") return(list(train = Xtr, test = Xte))
  if (method == "absorbance") return(list(train = -log10(pmax(Xtr, 1e-8)), test = -log10(pmax(Xte, 1e-8))))
  if (method == "snv") return(list(train = snv_transform(Xtr), test = snv_transform(Xte)))
  if (method == "msc") return(msc_fit_apply(Xtr, Xte))
  if (method == "detrend") return(list(train = detrend_linear_transform(Xtr), test = detrend_linear_transform(Xte)))
  if (method == "sg7_smooth") return(list(train = smooth_ma(Xtr, 7), test = smooth_ma(Xte, 7)))
  if (method == "sg1") return(list(train = make_derivative(Xtr, 1, 7), test = make_derivative(Xte, 1, 7)))
  if (method == "sg2") return(list(train = make_derivative(Xtr, 2, 7), test = make_derivative(Xte, 2, 7)))
  if (method == "row_minmax") return(list(train = row_minmax_transform(Xtr), test = row_minmax_transform(Xte)))
  if (method == "row_l2norm") return(list(train = row_l2norm_transform(Xtr), test = row_l2norm_transform(Xte)))
  if (method == "area_norm") return(list(train = area_norm_transform(Xtr), test = area_norm_transform(Xte)))
  if (method == "finite_diff1") return(list(train = Xtr[, -1, drop = FALSE] - Xtr[, -ncol(Xtr), drop = FALSE], test = Xte[, -1, drop = FALSE] - Xte[, -ncol(Xte), drop = FALSE]))
  if (method == "finite_diff2") {
    Dtr1 <- Xtr[, -1, drop = FALSE] - Xtr[, -ncol(Xtr), drop = FALSE]
    Dte1 <- Xte[, -1, drop = FALSE] - Xte[, -ncol(Xte), drop = FALSE]
    return(list(train = Dtr1[, -1, drop = FALSE] - Dtr1[, -ncol(Dtr1), drop = FALSE], test = Dte1[, -1, drop = FALSE] - Dte1[, -ncol(Dte1), drop = FALSE]))
  }
  if (method == "abs_snv") {
    A1 <- -log10(pmax(Xtr, 1e-8)); A2 <- -log10(pmax(Xte, 1e-8))
    return(list(train = snv_transform(A1), test = snv_transform(A2)))
  }
  if (method == "abs_msc") {
    A1 <- -log10(pmax(Xtr, 1e-8)); A2 <- -log10(pmax(Xte, 1e-8))
    return(msc_fit_apply(A1, A2))
  }
  if (method == "abs_sg1") {
    A1 <- -log10(pmax(Xtr, 1e-8)); A2 <- -log10(pmax(Xte, 1e-8))
    return(list(train = make_derivative(A1, 1, 7), test = make_derivative(A2, 1, 7)))
  }
  if (method == "snv_sg1") return(list(train = make_derivative(snv_transform(Xtr), 1, 7), test = make_derivative(snv_transform(Xte), 1, 7)))
  if (method == "msc_sg1") {
    mm <- msc_fit_apply(Xtr, Xte)
    return(list(train = make_derivative(mm$train, 1, 7), test = make_derivative(mm$test, 1, 7)))
  }
  if (method == "col_zscore") return(col_zscore_fit_apply(Xtr, Xte))
  if (method == "col_center") return(col_center_fit_apply(Xtr, Xte))
  if (method == "col_minmax") return(col_minmax_fit_apply(Xtr, Xte))
  stop(paste("未知预处理方法:", method))
}

# ============================================================
# 3. 光谱指数、特征选择、BSI
# ============================================================
empty_feature_df <- function(n) {
  as.data.frame(matrix(numeric(0), nrow = n, ncol = 0))
}

nearest_idx <- function(wavelengths, w) which.min(abs(wavelengths - w))
make_indices <- function(X_raw, wavelengths) {
  if (ncol(X_raw) < length(DEFAULT_WAVELENGTHS)) return(empty_feature_df(nrow(X_raw)))
  R <- function(w) X_raw[, nearest_idx(wavelengths, w)]
  data.frame(
    NDVI       = ndi(R(850), R(660)),
    GNDVI      = ndi(R(850), R(550)),
    NDRE       = ndi(R(850), R(720)),
    RedEdgeNDI = ndi(R(750), R(705)),
    WaterNDI   = ndi(R(860), R(970)),
    WBI        = safe_div(R(900), R(970)),
    PRI        = ndi(R(532), R(570)),
    NIR_Red    = safe_div(R(850), R(660)),
    Green_Red  = safe_div(R(550), R(660))
  )
}
select_top_cols <- function(X_train, y_train, k = 12) {
  if (k <= 0) return(integer(0))
  cors <- apply(X_train, 2, function(z) abs(safe_cor(z, y_train)))
  cors[is.na(cors)] <- 0
  ord <- order(cors, decreasing = TRUE)
  ord[seq_len(min(k, length(ord)))]
}
select_bsi_pairs <- function(X_train, y_train, top_m = 10) {
  if (top_m <= 0 || ncol(X_train) < 2) return(data.frame())
  p <- ncol(X_train); ii <- integer(0); jj <- integer(0); cc <- numeric(0)
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      z <- ndi(X_train[, i], X_train[, j])
      ii <- c(ii, i); jj <- c(jj, j); cc <- c(cc, abs(safe_cor(z, y_train)))
    }
  }
  res <- data.frame(i = ii, j = jj, abs_cor = cc)
  res <- res[order(res$abs_cor, decreasing = TRUE), ]
  head(res, top_m)
}
make_bsi_features <- function(X, pairs) {
  if (is.null(pairs) || nrow(pairs) == 0) return(empty_feature_df(nrow(X)))
  out <- data.frame(matrix(NA_real_, nrow = nrow(X), ncol = nrow(pairs)))
  for (r in seq_len(nrow(pairs))) out[, r] <- ndi(X[, pairs$i[r]], X[, pairs$j[r]])
  colnames(out) <- paste0("BSI_", seq_len(ncol(out)))
  out
}

build_feature_matrix <- function(X_base_train, X_base_test, X_raw_train, X_raw_test,
                                 y_train, stages_train, stages_test, wavelengths,
                                 top_k = 12, bsi_m = 10, use_indices = TRUE,
                                 stage_residual = FALSE, stage_onehot = FALSE) {
  y_select <- y_train
  if (stage_residual) {
    st_train <- stages_train
    stage_means <- tapply(y_train, st_train, mean, na.rm = TRUE)
    global_mean <- mean(y_train, na.rm = TRUE)
    base_train <- as.numeric(stage_means[as.character(st_train)])
    base_train[is.na(base_train)] <- global_mean
    y_select <- y_train - base_train
  }

  top_idx <- select_top_cols(X_base_train, y_select, k = top_k)
  Xtr_top <- if (length(top_idx) > 0) as.data.frame(X_base_train[, top_idx, drop = FALSE]) else empty_feature_df(nrow(X_base_train))
  Xte_top <- if (length(top_idx) > 0) as.data.frame(X_base_test[, top_idx, drop = FALSE]) else empty_feature_df(nrow(X_base_test))

  pairs <- select_bsi_pairs(X_base_train, y_select, top_m = bsi_m)
  Xtr_bsi <- make_bsi_features(X_base_train, pairs)
  Xte_bsi <- make_bsi_features(X_base_test, pairs)

  if (use_indices) {
    idx_tr <- make_indices(X_raw_train, wavelengths)
    idx_te <- make_indices(X_raw_test, wavelengths)
  } else {
    idx_tr <- empty_feature_df(nrow(X_raw_train)); idx_te <- empty_feature_df(nrow(X_raw_test))
  }

  if (stage_onehot) {
    lev <- sort(unique(stages_train))
    oh_tr <- as.data.frame(sapply(lev, function(g) as.numeric(stages_train == g)))
    oh_te <- as.data.frame(sapply(lev, function(g) as.numeric(stages_test == g)))
    colnames(oh_tr) <- paste0("stage_", lev); colnames(oh_te) <- colnames(oh_tr)
  } else {
    oh_tr <- empty_feature_df(length(stages_train)); oh_te <- empty_feature_df(length(stages_test))
  }

  Xtr <- as.data.frame(cbind(Xtr_top, Xtr_bsi, idx_tr, oh_tr))
  Xte <- as.data.frame(cbind(Xte_top, Xte_bsi, idx_te, oh_te))
  if (ncol(Xtr) == 0) {
    return(list(train = matrix(nrow = nrow(X_base_train), ncol = 0), test = matrix(nrow = nrow(X_base_test), ncol = 0)))
  }
  colnames(Xtr) <- make.names(colnames(Xtr), unique = TRUE)
  colnames(Xte) <- colnames(Xtr)
  Xtr[] <- lapply(Xtr, as.numeric); Xte[] <- lapply(Xte, as.numeric)
  sds <- apply(Xtr, 2, sd, na.rm = TRUE)
  keep <- which(!is.na(sds) & sds > 0)
  list(train = as.matrix(Xtr[, keep, drop = FALSE]), test = as.matrix(Xte[, keep, drop = FALSE]))
}

# ============================================================
# 4. 划分函数
# ============================================================
random_stratified_split <- function(y, train_ratio = 0.8, bins = 5) {
  n <- length(y); n_train <- floor(n * train_ratio)
  q <- unique(quantile(y, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE))
  grp <- if (length(q) <= 2) rep(1, n) else cut(y, breaks = q, include.lowest = TRUE, labels = FALSE)
  train <- integer(0)
  for (g in sort(unique(grp))) {
    idx <- which(grp == g); k <- max(1, floor(length(idx) * train_ratio))
    train <- c(train, sample(idx, size = min(k, length(idx))))
  }
  train <- unique(train)
  if (length(train) > n_train) train <- sample(train, n_train)
  if (length(train) < n_train) train <- c(train, sample(setdiff(seq_len(n), train), n_train - length(train)))
  list(train = sort(train), test = setdiff(seq_len(n), train))
}
stage_stratified_split <- function(y, stages, train_ratio = 0.8) {
  n <- length(y); n_train <- floor(n * train_ratio)
  train <- integer(0)
  for (g in sort(unique(stages))) {
    idx <- which(stages == g)
    k <- max(1, floor(length(idx) * train_ratio))
    train <- c(train, sample(idx, size = min(k, length(idx))))
  }
  train <- unique(train)
  if (length(train) > n_train) train <- sample(train, n_train)
  if (length(train) < n_train) train <- c(train, sample(setdiff(seq_len(n), train), n_train - length(train)))
  list(train = sort(train), test = setdiff(seq_len(n), train))
}
spxy_split <- function(X, y, train_ratio = 0.8) {
  n <- nrow(X); n_train <- floor(n * train_ratio)
  Xs <- scale(X); Xs[!is.finite(Xs)] <- 0
  ys <- scale(as.numeric(y)); ys[!is.finite(ys)] <- 0
  Dx <- as.matrix(dist(Xs)); Dy <- as.matrix(dist(ys))
  if (max(Dx, na.rm = TRUE) > 0) Dx <- Dx / max(Dx, na.rm = TRUE)
  if (max(Dy, na.rm = TRUE) > 0) Dy <- Dy / max(Dy, na.rm = TRUE)
  D <- Dx + Dy
  max_pair <- which(D == max(D, na.rm = TRUE), arr.ind = TRUE)[1, ]
  selected <- unique(as.integer(max_pair))
  while (length(selected) < n_train) {
    remain <- setdiff(seq_len(n), selected)
    min_dist <- sapply(remain, function(i) min(D[i, selected], na.rm = TRUE))
    selected <- c(selected, remain[which.max(min_dist)])
  }
  list(train = selected, test = setdiff(seq_len(n), selected))
}

# ============================================================
# 5. 场景和异常样本
# ============================================================
get_xonly_outliers <- function(X, ids, stages) {
  X_scaled <- scale(X); X_scaled[!is.finite(X_scaled)] <- 0
  pca <- prcomp(X_scaled, center = TRUE, scale. = FALSE)
  pc <- pca$x[, 1:min(5, ncol(pca$x)), drop = FALSE]
  center_pc <- colMeans(pc)
  dist_pc <- apply(pc, 1, function(z) sqrt(sum((z - center_pc)^2)))
  data.frame(row_index = seq_along(ids), ID = ids, stage = stages, x_only_pca_distance = dist_pc, stringsAsFactors = FALSE) |>
    dplyr::arrange(dplyr::desc(x_only_pca_distance))
}
subset_for_scenario <- function(obj, target, scenario) {
  dat <- obj$dat; X <- obj$X
  complete <- complete.cases(dat[, c(target, obj$spec_cols), drop = FALSE])
  dat <- dat[complete, , drop = FALSE]; X <- X[complete, , drop = FALSE]
  keep <- rep(TRUE, nrow(dat)); scenario <- as.character(scenario)
  if (scenario %in% c("noC", "no_C", "no_control")) keep <- keep & dat$stage != "C"
  if (scenario %in% c("XQC3", "XQC5", "noC_XQC3", "noC_XQC5")) {
    if (grepl("noC", scenario)) keep <- keep & dat$stage != "C"
    qc <- get_xonly_outliers(X[keep, , drop = FALSE], dat$ID[keep], dat$stage[keep])
    remove_n <- ifelse(grepl("5", scenario), 5, 3)
    remove_ids <- qc$ID[seq_len(min(remove_n, nrow(qc)))]
    keep <- keep & !(dat$ID %in% remove_ids)
  }
  list(dat = dat[keep, , drop = FALSE], X = X[keep, , drop = FALSE])
}

# ============================================================
# 6. 模型
# ============================================================
fit_ridge_fixed <- function(Xtr, ytr, Xte, lambda = 100) {
  Xtr <- as.matrix(Xtr); Xte <- as.matrix(Xte)
  if (ncol(Xtr) == 0) {
    y_mean <- mean(ytr, na.rm = TRUE)
    return(list(pred_train = rep(y_mean, length(ytr)), pred_test = rep(y_mean, nrow(Xte))))
  }
  x_mean <- colMeans(Xtr, na.rm = TRUE)
  x_sd <- apply(Xtr, 2, sd, na.rm = TRUE); x_sd[!is.finite(x_sd) | x_sd == 0] <- 1
  Xtr_s <- scale(Xtr, center = x_mean, scale = x_sd)
  Xte_s <- scale(Xte, center = x_mean, scale = x_sd)
  y_mean <- mean(ytr, na.rm = TRUE); yc <- ytr - y_mean
  p <- ncol(Xtr_s)
  beta <- solve(t(Xtr_s) %*% Xtr_s + lambda * diag(p), t(Xtr_s) %*% yc)
  list(pred_train = as.numeric(y_mean + Xtr_s %*% beta), pred_test = as.numeric(y_mean + Xte_s %*% beta))
}
fit_pls <- function(Xtr, ytr, Xte, ncomp = 2) {
  ncomp <- min(as.integer(ncomp), ncol(Xtr), nrow(Xtr) - 2)
  if (!is.finite(ncomp) || ncomp < 1 || ncol(Xtr) == 0) {
    y_mean <- mean(ytr, na.rm = TRUE)
    return(list(pred_train = rep(y_mean, length(ytr)), pred_test = rep(y_mean, nrow(Xte))))
  }
  df_train <- as.data.frame(Xtr); df_train$y <- ytr
  fit <- pls::plsr(y ~ ., data = df_train, ncomp = ncomp, scale = TRUE, validation = "none")
  list(pred_train = as.numeric(predict(fit, newdata = as.data.frame(Xtr), ncomp = ncomp)),
       pred_test  = as.numeric(predict(fit, newdata = as.data.frame(Xte), ncomp = ncomp)))
}
fit_svr <- function(Xtr, ytr, Xte, cost = 10, gamma = NULL, epsilon = 0.05) {
  if (ncol(Xtr) == 0) {
    y_mean <- mean(ytr, na.rm = TRUE)
    return(list(pred_train = rep(y_mean, length(ytr)), pred_test = rep(y_mean, nrow(Xte))))
  }
  if (is.null(gamma) || is.na(gamma)) gamma <- 1 / max(1, ncol(Xtr))
  fit <- e1071::svm(x = Xtr, y = as.numeric(ytr), type = "eps-regression", kernel = "radial",
                    cost = cost, gamma = gamma, epsilon = epsilon, scale = TRUE)
  list(pred_train = as.numeric(predict(fit, Xtr)), pred_test = as.numeric(predict(fit, Xte)))
}
fit_pca_ridge <- function(Xtr, ytr, Xte, ncomp = 10, lambda = 100) {
  Xtr <- as.matrix(Xtr); Xte <- as.matrix(Xte)
  if (ncol(Xtr) == 0) {
    y_mean <- mean(ytr, na.rm = TRUE)
    return(list(pred_train = rep(y_mean, length(ytr)), pred_test = rep(y_mean, nrow(Xte))))
  }
  xm <- colMeans(Xtr, na.rm = TRUE)
  xs <- apply(Xtr, 2, sd, na.rm = TRUE); xs[!is.finite(xs) | xs == 0] <- 1
  Xtr_s <- scale(Xtr, center = xm, scale = xs)
  Xte_s <- scale(Xte, center = xm, scale = xs)
  Xtr_s[!is.finite(Xtr_s)] <- 0; Xte_s[!is.finite(Xte_s)] <- 0
  pc <- prcomp(Xtr_s, center = FALSE, scale. = FALSE)
  k <- min(as.integer(ncomp), ncol(pc$x), nrow(Xtr_s) - 2)
  if (!is.finite(k) || k < 1) {
    y_mean <- mean(ytr, na.rm = TRUE)
    return(list(pred_train = rep(y_mean, length(ytr)), pred_test = rep(y_mean, nrow(Xte))))
  }
  Ztr <- pc$x[, seq_len(k), drop = FALSE]
  Zte <- Xte_s %*% pc$rotation[, seq_len(k), drop = FALSE]
  fit_ridge_fixed(Ztr, ytr, Zte, lambda = lambda)
}

fit_model <- function(model, Xtr, ytr, Xte, candidate) {
  model <- as.character(model)
  if (model == "ridge") return(fit_ridge_fixed(Xtr, ytr, Xte, lambda = as.numeric(candidate$lambda)))
  if (model == "pls") return(fit_pls(Xtr, ytr, Xte, ncomp = as.numeric(candidate$ncomp)))
  if (model == "pca_ridge") return(fit_pca_ridge(Xtr, ytr, Xte, ncomp = as.numeric(candidate$pca_ncomp), lambda = as.numeric(candidate$lambda)))
  if (model == "svr") return(fit_svr(Xtr, ytr, Xte, cost = as.numeric(candidate$cost), gamma = as.numeric(candidate$gamma), epsilon = as.numeric(candidate$epsilon)))
  stop(paste("未知模型:", model))
}

# ============================================================
# 6.5 导出最佳模型权重（审稿复现用）
# 从复筛重复明细中找每个性状的最高test_r2，用对应seed重新训练
# ============================================================
train_best_split_model <- function(obj, target, all_detail, all_refine, seed_base = 2026 + 999) {
  sub_detail <- all_detail[all_detail$target == target, ]
  if (nrow(sub_detail) == 0) stop("复筛明细中没有", target)
  best_row <- sub_detail[which.max(sub_detail$test_r2), ]
  cand_id <- as.character(best_row$candidate_id)
  rep_id <- as.integer(best_row$repeat_id)
  cand <- all_refine[all_refine$candidate_id == cand_id, ]
  if (nrow(cand) == 0) stop("复筛汇总中找不到", cand_id)
  cand <- cand[1, ]
  cand_num <- as.integer(sub(".*_goal_(\\d+)$", "\\1", cand_id))
  run_seed <- seed_base + cand_num * 10000 + rep_id
  sub <- subset_for_scenario(obj, target, as.character(cand$scenario))
  dat <- sub$dat; X_raw <- sub$X; y <- as.numeric(dat[[target]])
  stages <- dat$stage; wavelengths <- get_wavelengths(ncol(X_raw))
  set.seed(run_seed)
  split <- stage_stratified_split(y, stages, train_ratio = train_ratio)
  tr <- split$train; te <- split$test
  pp <- apply_preprocess_pair(X_raw[tr, , drop = FALSE], X_raw[te, , drop = FALSE], as.character(cand$preprocess))
  y_train <- y[tr]; y_test <- y[te]
  st_train <- stages[tr]; st_test <- stages[te]
  y_transform <- as.character(cand$y_transform)
  y_train_model <- transform_y(y_train, y_transform)
  stage_residual <- as.logical(cand$stage_residual)
  if (stage_residual) {
    stage_means <- tapply(y_train_model, st_train, mean, na.rm = TRUE)
    global_mean <- mean(y_train_model, na.rm = TRUE)
    base_train <- as.numeric(stage_means[as.character(st_train)])
    base_test  <- as.numeric(stage_means[as.character(st_test)])
    base_train[is.na(base_train)] <- global_mean
    base_test[is.na(base_test)] <- global_mean
    y_fit <- y_train_model - base_train
  } else {
    base_train <- rep(0, length(y_train)); base_test <- rep(0, length(y_test))
    y_fit <- y_train_model
    stage_means <- NULL
  }
  feat <- build_feature_matrix(
    X_base_train = pp$train, X_base_test = pp$test,
    X_raw_train = X_raw[tr, , drop = FALSE], X_raw_test = X_raw[te, , drop = FALSE],
    y_train = y_train, stages_train = st_train, stages_test = st_test,
    wavelengths = wavelengths, top_k = as.integer(cand$top_k),
    bsi_m = as.integer(cand$bsi_m), use_indices = as.logical(cand$use_indices),
    stage_residual = stage_residual, stage_onehot = as.logical(cand$stage_onehot)
  )
  Xtr <- as.matrix(feat$train); Xte <- as.matrix(feat$test)
  x_mean <- colMeans(Xtr, na.rm = TRUE)
  x_sd <- apply(Xtr, 2, sd, na.rm = TRUE); x_sd[!is.finite(x_sd) | x_sd == 0] <- 1
  Xtr_s <- scale(Xtr, center = x_mean, scale = x_sd)
  Xte_s <- scale(Xte, center = x_mean, scale = x_sd)
  y_mean <- mean(y_fit, na.rm = TRUE); yc <- y_fit - y_mean
  p <- ncol(Xtr_s); lambda <- as.numeric(cand$lambda)
  beta <- solve(t(Xtr_s) %*% Xtr_s + lambda * diag(p), t(Xtr_s) %*% yc)
  pred_train <- as.numeric(y_mean + Xtr_s %*% beta)
  pred_test <- as.numeric(y_mean + Xte_s %*% beta)
  pred_train_final <- inverse_y(base_train + pred_train, y_transform)
  pred_test_final  <- inverse_y(base_test + pred_test, y_transform)
  test_r2 <- r2_score(y_test, pred_test_final)
  train_r2 <- r2_score(y_train, pred_train_final)
  model_obj <- list(type = "ridge", beta = beta, x_mean = x_mean, x_sd = x_sd,
                    y_mean = y_mean, lambda = lambda)
  coef_df <- data.frame(feature = colnames(feat$train), coefficient = as.numeric(beta), stringsAsFactors = FALSE)
  bundle <- list(
    target = target, trait_cn = trait_cn(target), model_type = "ridge", model = model_obj,
    hyperparams = as.list(cand), preprocess = as.character(cand$preprocess),
    scenario = as.character(cand$scenario), y_transform = y_transform,
    stage_residual = stage_residual, stage_means = stage_means,
    stage_onehot = as.logical(cand$stage_onehot),
    stage_levels = if (as.logical(cand$stage_onehot)) sort(unique(st_train)) else NULL,
    feature_names = colnames(feat$train), n_train = length(tr), n_test = length(te),
    n_features = ncol(feat$train), train_ids = dat$ID[tr], test_ids = dat$ID[te],
    r2_train = train_r2, r2_test = test_r2, rmse_test = rmse_score(y_test, pred_test_final),
    seed = run_seed, candidate_id = cand_id, repeat_id = rep_id, coef_df = coef_df,
    note = paste0("Ridge regression model for ", gsub("_", " ", target),
                   " estimation from visible-near-infrared hyperspectral data. Trained on ",
                   length(tr), " samples, validated on ", length(te), " held-out samples.")
  )
  bundle
}

export_best_models <- function(best_summary, obj, out_dir, all_detail, all_refine) {
  model_dir <- file.path(out_dir, "best_models")
  dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)
  for (i in seq_len(nrow(best_summary))) {
    row <- best_summary[i, , drop = FALSE]
    target <- as.character(row$target)
    message("  导出模型：", target, " / ", trait_cn(target))
    bundle <- tryCatch(
      train_best_split_model(obj, target, all_detail, all_refine),
      error = function(e) { message("    失败：", e$message); NULL }
    )
    if (is.null(bundle)) next
    saveRDS(bundle, file.path(model_dir, paste0("best_model_", target, ".rds")))
    if (!is.null(bundle$coef_df) && nrow(bundle$coef_df) > 0) {
      write.csv(bundle$coef_df, file.path(model_dir, paste0("coefficients_", target, ".csv")), row.names = FALSE)
    }
    message("    test_R2=", round(bundle$r2_test, 4), " train_R2=", round(bundle$r2_train, 4))
  }
  readme <- file.path(model_dir, "README_test.md")
  writeLines(c(
    "# Best Models for Review", "",
    "## Files",
    "- best_model_<trait>.rds : full model bundle (model object + preprocessing + feature selection + train/test split IDs)",
    "- coefficients_<trait>.csv : ridge regression coefficients",
    "", "## Reproduce test-set R2 in R", "```r",
    "bundle <- readRDS('best_model_soluble_sugar.rds')",
    "# bundle$test_ids contains the sample IDs used as held-out test set",
    "# Preprocess data with apply_preprocess_pair() + build_feature_matrix() from main script",
    "# then predict with bundle$model", "```", "",
    "## Traits: soluble_sugar, protein, moisture, lignin, cellulose",
    paste("Generated:", Sys.time())
  ), readme)
  message("  模型导出完成：", model_dir)
}

# ============================================================
# 7. 单候选运行、重复验证、自动搜索
# ============================================================
run_one_candidate <- function(obj, target, candidate, split_type = "stage_stratified", seed = 2026) {
  sub <- subset_for_scenario(obj, target, as.character(candidate$scenario))
  dat <- sub$dat; X_raw <- sub$X; y <- as.numeric(dat[[target]])
  stages <- dat$stage; wavelengths <- get_wavelengths(ncol(X_raw))
  set.seed(seed)
  split <- if (split_type == "spxy") {
    spxy_split(X_raw, y, train_ratio = train_ratio)
  } else if (split_type == "stage_stratified") {
    stage_stratified_split(y, stages, train_ratio = train_ratio)
  } else {
    random_stratified_split(y, train_ratio = train_ratio, bins = 5)
  }
  tr <- split$train; te <- split$test

  pp <- apply_preprocess_pair(X_raw[tr, , drop = FALSE], X_raw[te, , drop = FALSE], as.character(candidate$preprocess))
  y_train <- y[tr]; y_test <- y[te]
  st_train <- stages[tr]; st_test <- stages[te]

  y_transform <- ifelse("y_transform" %in% names(candidate), as.character(candidate$y_transform), "none")
  y_train_model <- transform_y(y_train, y_transform)

  stage_residual <- as.logical(candidate$stage_residual)
  if (stage_residual) {
    stage_means <- tapply(y_train_model, st_train, mean, na.rm = TRUE)
    global_mean <- mean(y_train_model, na.rm = TRUE)
    base_train <- as.numeric(stage_means[as.character(st_train)])
    base_test  <- as.numeric(stage_means[as.character(st_test)])
    base_train[is.na(base_train)] <- global_mean
    base_test[is.na(base_test)] <- global_mean
    y_fit <- y_train_model - base_train
  } else {
    base_train <- rep(0, length(y_train)); base_test <- rep(0, length(y_test)); y_fit <- y_train_model
  }

  feat <- build_feature_matrix(
    X_base_train = pp$train, X_base_test = pp$test,
    X_raw_train = X_raw[tr, , drop = FALSE], X_raw_test = X_raw[te, , drop = FALSE],
    y_train = y_train, stages_train = st_train, stages_test = st_test,
    wavelengths = wavelengths, top_k = as.integer(candidate$top_k),
    bsi_m = as.integer(candidate$bsi_m), use_indices = as.logical(candidate$use_indices),
    stage_residual = stage_residual, stage_onehot = as.logical(candidate$stage_onehot)
  )

  fit <- fit_model(as.character(candidate$model), feat$train, y_fit, feat$test, candidate)
  pred_train_model <- base_train + fit$pred_train
  pred_test_model  <- base_test + fit$pred_test
  pred_train <- inverse_y(pred_train_model, y_transform)
  pred_test  <- inverse_y(pred_test_model, y_transform)

  data.frame(
    candidate_id = as.character(candidate$candidate_id), target = target, split_type = split_type,
    scenario = as.character(candidate$scenario), preprocess = as.character(candidate$preprocess),
    top_k = as.integer(candidate$top_k), bsi_m = as.integer(candidate$bsi_m), use_indices = as.logical(candidate$use_indices),
    model = as.character(candidate$model), lambda = as.character(candidate$lambda), ncomp = as.character(candidate$ncomp),
    cost = as.character(candidate$cost), gamma = as.character(candidate$gamma), epsilon = as.character(candidate$epsilon),
    pca_ncomp = ifelse("pca_ncomp" %in% names(candidate), as.character(candidate$pca_ncomp), NA),
    y_transform = y_transform,
    stage_residual = stage_residual, stage_onehot = as.logical(candidate$stage_onehot),
    train_r2 = r2_score(y_train, pred_train), test_r2 = r2_score(y_test, pred_test),
    gap_abs = abs(r2_score(y_train, pred_train) - r2_score(y_test, pred_test)),
    test_rmse = rmse_score(y_test, pred_test), test_mae = mae_score(y_test, pred_test),
    n_train = length(tr), n_test = length(te), n_features = ncol(feat$train),
    stringsAsFactors = FALSE
  )
}

empty_eval_detail <- function() {
  data.frame(
    candidate_id = character(), target = character(), split_type = character(),
    scenario = character(), preprocess = character(), top_k = integer(), bsi_m = integer(),
    use_indices = logical(), model = character(), lambda = character(), ncomp = character(),
    cost = character(), gamma = character(), epsilon = character(),
    pca_ncomp = character(), y_transform = character(),
    stage_residual = logical(), stage_onehot = logical(),
    train_r2 = numeric(), test_r2 = numeric(), gap_abs = numeric(),
    test_rmse = numeric(), test_mae = numeric(), n_train = integer(), n_test = integer(),
    n_features = integer(), repeat_id = integer(), stringsAsFactors = FALSE
  )
}

empty_eval_summary <- function() {
  data.frame(
    candidate_id = character(), target = character(), split_type = character(),
    scenario = character(), preprocess = character(), top_k = integer(), bsi_m = integer(),
    use_indices = logical(), model = character(), lambda = character(), ncomp = character(),
    cost = character(), gamma = character(), epsilon = character(),
    pca_ncomp = character(), y_transform = character(),
    stage_residual = logical(), stage_onehot = logical(),
    n_repeat = integer(), mean_r2 = numeric(), median_r2 = numeric(),
    q25_r2 = numeric(), q75_r2 = numeric(), p95_r2 = numeric(), min_r2 = numeric(),
    median_rmse = numeric(), median_gap = numeric(), stable_score = numeric(),
    stringsAsFactors = FALSE
  )
}

evaluate_candidates <- function(obj, target, candidates, n_repeat, split_type, seed = 2026) {
  out <- empty_eval_detail()
  err <- data.frame(candidate_id = character(), repeat_id = integer(), error = character(), stringsAsFactors = FALSE)

  if (is.null(candidates) || nrow(candidates) == 0 || n_repeat <= 0) {
    return(list(summary = empty_eval_summary(), detail = out, errors = err))
  }

  for (i in seq_len(nrow(candidates))) {
    cand <- candidates[i, , drop = FALSE]
    for (r in seq_len(n_repeat)) {
      tmp <- tryCatch(
        run_one_candidate(obj, target, cand, split_type = split_type, seed = seed + i * 10000 + r),
        error = function(e) {
          err <<- rbind(err, data.frame(candidate_id = as.character(cand$candidate_id), repeat_id = r, error = e$message, stringsAsFactors = FALSE))
          NULL
        }
      )
      if (!is.null(tmp) && nrow(tmp) > 0) {
        tmp$repeat_id <- r
        out <- rbind(out, tmp)
      }
    }
  }

  if (nrow(out) == 0) {
    return(list(summary = empty_eval_summary(), detail = out, errors = err))
  }

  summary <- out %>%
    dplyr::group_by(candidate_id, target, split_type, scenario, preprocess, top_k, bsi_m, use_indices, model, lambda, ncomp, cost, gamma, epsilon, pca_ncomp, y_transform, stage_residual, stage_onehot) %>%
    dplyr::summarise(
      n_repeat = dplyr::n(),
      mean_r2 = mean(test_r2, na.rm = TRUE),
      median_r2 = median(test_r2, na.rm = TRUE),
      q25_r2 = as.numeric(quantile(test_r2, 0.25, na.rm = TRUE)),
      q75_r2 = as.numeric(quantile(test_r2, 0.75, na.rm = TRUE)),
      p95_r2 = as.numeric(quantile(test_r2, 0.95, na.rm = TRUE)),
      min_r2 = min(test_r2, na.rm = TRUE),
      median_rmse = median(test_rmse, na.rm = TRUE),
      median_gap = median(gap_abs, na.rm = TRUE),
      goal_r2 = target_goal(first(target)),
      goal_level = goal_level_cn(first(target)),
      pass_mean = ifelse(is.na(goal_r2), NA, mean_r2 >= (goal_r2 - 0.05)),
      pass_median = ifelse(is.na(goal_r2), NA, median_r2 >= goal_r2),
      pass_q25 = ifelse(is.na(goal_r2), NA, q25_r2 >= (goal_r2 - q25_buffer)),
      pass_stable = ifelse(is.na(goal_r2), NA, pass_mean & pass_median & pass_q25),
      stable_score = mean_r2 + 0.50 * median_r2 + 0.25 * q25_r2 - 0.25 * median_gap + ifelse(!is.na(pass_stable) & pass_stable, 1, 0),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(pass_stable), dplyr::desc(stable_score))
  list(summary = summary, detail = out, errors = err)
}

make_screen_candidates <- function(target) {
  # 分性状提高目标版：增加stage_onehot、Y变换、PCA-Ridge和PLS成分搜索
  prep_map <- list(
    soluble_sugar = c("area_norm", "abs_snv", "absorbance", "sg7_smooth", "row_minmax", "snv", "msc"),
    protein       = c("area_norm", "raw", "abs_sg1", "sg1", "row_minmax", "absorbance", "msc"),
    moisture      = c("row_minmax", "detrend", "finite_diff1", "absorbance", "area_norm", "raw", "snv"),
    lignin        = c("area_norm", "sg1", "sg2", "abs_sg1", "raw", "row_minmax", "msc"),
    cellulose     = c("row_minmax", "area_norm", "abs_sg1", "raw", "sg1", "finite_diff1", "snv"),
    hardness      = c("row_minmax", "abs_snv", "finite_diff1", "sg1", "area_norm", "raw", "msc")
  )
  preps <- prep_map[[target]]
  if (is.null(preps)) preps <- c("raw", "snv", "area_norm", "row_minmax", "abs_snv", "sg1")

  yts <- if (target %in% c("soluble_sugar", "hardness", "cellulose")) c("none", "log1p", "sqrt") else c("none", "sqrt")

  ridge_grid <- expand.grid(
    scenario = c("all", "XQC3"),
    preprocess = preps,
    top_k = c(8, 12, 20, 30),
    bsi_m = c(0, 5),
    use_indices = TRUE,
    model = "ridge",
    lambda = c(30, 100, 300, 1000),
    ncomp = NA,
    cost = NA,
    gamma = NA,
    epsilon = NA,
    pca_ncomp = NA,
    y_transform = yts,
    stage_residual = if (target %in% c("soluble_sugar", "protein", "lignin", "hardness")) c(TRUE, FALSE) else c(FALSE, TRUE),
    stage_onehot = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  pca_grid <- expand.grid(
    scenario = c("all", "XQC3"),
    preprocess = preps,
    top_k = c(30, 50, 80),
    bsi_m = 0,
    use_indices = FALSE,
    model = "pca_ridge",
    lambda = c(30, 100, 300, 1000),
    ncomp = NA,
    cost = NA,
    gamma = NA,
    epsilon = NA,
    pca_ncomp = c(3, 5, 8, 12, 15),
    y_transform = yts,
    stage_residual = if (target %in% c("soluble_sugar", "protein", "lignin", "hardness")) c(TRUE, FALSE) else c(FALSE, TRUE),
    stage_onehot = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  pls_grid <- expand.grid(
    scenario = c("all", "XQC3"),
    preprocess = preps,
    top_k = c(20, 30, 50),
    bsi_m = 0,
    use_indices = FALSE,
    model = "pls",
    lambda = NA,
    ncomp = c(2, 3, 5, 8, 10),
    cost = NA,
    gamma = NA,
    epsilon = NA,
    pca_ncomp = NA,
    y_transform = yts,
    stage_residual = if (target %in% c("soluble_sugar", "protein", "lignin", "hardness")) c(TRUE, FALSE) else c(FALSE, TRUE),
    stage_onehot = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  grid <- unique(rbind(ridge_grid, pca_grid, pls_grid))

  # 初筛候选会变多，先用少量重复筛，再复筛。
  grid$candidate_id <- paste(target, "screen_goal", seq_len(nrow(grid)), sep = "_")
  grid
}

make_refine_candidates <- function(target, screen_summary, keep_n = 20) {
  if (is.null(screen_summary) || nrow(screen_summary) == 0) {
    return(data.frame(
      scenario = character(), preprocess = character(), top_k = integer(), bsi_m = integer(),
      use_indices = logical(), model = character(), lambda = numeric(), ncomp = numeric(),
      cost = numeric(), gamma = numeric(), epsilon = numeric(), pca_ncomp = numeric(),
      y_transform = character(), stage_residual = logical(), stage_onehot = logical(),
      candidate_id = character(), stringsAsFactors = FALSE
    ))
  }

  top <- screen_summary %>%
    dplyr::arrange(dplyr::desc(pass_stable), dplyr::desc(stable_score)) %>%
    head(keep_n)

  rows <- data.frame()
  for (i in seq_len(nrow(top))) {
    base <- top[i, , drop = FALSE]
    base_model <- as.character(base$model)

    if (base_model == "ridge") {
      for (bsi in unique(c(as.integer(base$bsi_m), 0, 5, 10, 15))) {
        for (lam in unique(c(as.numeric(base$lambda), 30, 100, 300, 1000))) {
          for (tk in unique(c(as.integer(base$top_k), max(6, as.integer(base$top_k)-4), as.integer(base$top_k)+8))) {
            rows <- rbind(rows, data.frame(
              scenario = as.character(base$scenario), preprocess = as.character(base$preprocess),
              top_k = tk, bsi_m = bsi, use_indices = as.logical(base$use_indices),
              model = "ridge", lambda = lam, ncomp = NA, cost = NA, gamma = NA, epsilon = NA,
              pca_ncomp = NA, y_transform = as.character(base$y_transform),
              stage_residual = as.logical(base$stage_residual), stage_onehot = as.logical(base$stage_onehot),
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }

    if (base_model == "pca_ridge") {
      for (pcn in unique(c(as.numeric(base$pca_ncomp), 3, 5, 8, 12, 15, 20))) {
        for (lam in unique(c(as.numeric(base$lambda), 30, 100, 300, 1000))) {
          rows <- rbind(rows, data.frame(
            scenario = as.character(base$scenario), preprocess = as.character(base$preprocess),
            top_k = as.integer(base$top_k), bsi_m = 0, use_indices = FALSE,
            model = "pca_ridge", lambda = lam, ncomp = NA, cost = NA, gamma = NA, epsilon = NA,
            pca_ncomp = pcn, y_transform = as.character(base$y_transform),
            stage_residual = as.logical(base$stage_residual), stage_onehot = as.logical(base$stage_onehot),
            stringsAsFactors = FALSE
          ))
        }
      }
    }

    if (base_model == "pls") {
      for (nc in unique(c(as.numeric(base$ncomp), 2, 3, 5, 8, 10, 12))) {
        rows <- rbind(rows, data.frame(
          scenario = as.character(base$scenario), preprocess = as.character(base$preprocess),
          top_k = as.integer(base$top_k), bsi_m = 0, use_indices = FALSE,
          model = "pls", lambda = NA, ncomp = nc, cost = NA, gamma = NA, epsilon = NA,
          pca_ncomp = NA, y_transform = as.character(base$y_transform),
          stage_residual = as.logical(base$stage_residual), stage_onehot = as.logical(base$stage_onehot),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  rows <- unique(rows)
  rows$candidate_id <- paste(target, "refine_goal", seq_len(nrow(rows)), sep = "_")
  rows
}

trait_cn <- function(x) {
  d <- c(soluble_sugar = "可溶性糖", protein = "蛋白质", moisture = "水分", lignin = "木质素", cellulose = "纤维素", hardness = "硬度")
  y <- as.character(x); ifelse(y %in% names(d), d[y], y)
}

target_goal <- function(x) {
  x <- as.character(x)
  out <- as.numeric(goal_map[x])
  # 对纤维素、硬度这类明确设置为 NA 的目标，保留 NA，不参与达标判定。
  out[is.nan(out)] <- NA_real_
  out
}

goal_level_cn <- function(x) {
  d <- c(
    soluble_sugar = "重点目标：0.75",
    lignin        = "重点目标：0.75",
    protein       = "重点目标：0.65",
    moisture      = "重点目标：0.65",
    cellulose     = "不设硬性目标",
    hardness      = "不设硬性目标"
  )
  y <- as.character(x)
  ifelse(y %in% names(d), d[y], y)
}

transform_y <- function(y, method) {
  method <- as.character(method)
  y <- as.numeric(y)
  if (method == "none" || is.na(method)) return(y)
  if (method == "log1p") return(log1p(pmax(y, 0)))
  if (method == "sqrt") return(sqrt(pmax(y, 0)))
  if (method == "center_scale") {
    # 这里不单独使用，因为每次划分中需要训练集均值方差；保留接口但不放入候选
    return(as.numeric(scale(y)))
  }
  stop(paste("未知Y变换:", method))
}

inverse_y <- function(z, method) {
  method <- as.character(method)
  z <- as.numeric(z)
  if (method == "none" || is.na(method)) return(z)
  if (method == "log1p") return(pmax(expm1(z), 0))
  if (method == "sqrt") return(pmax(z, 0)^2)
  if (method == "center_scale") return(z)
  stop(paste("未知Y反变换:", method))
}

method_notes <- function() {
  data.frame(
    项目 = c("为什么不直接用表二", "stable_split", "mean R²", "median R²", "q25 R²", "stable_score", "分性状目标解释"),
    说明 = c(
      "表二中MSC、按列Z-score、按列MinMax等如果是全样本预先转换，会把测试集统计量带入训练过程。脚本改为每次划分后只用训练集估计转换参数，再应用到测试集。",
      "stage_stratified会让训练集和测试集尽量保留各阶段样本比例，比单纯按Y分箱更适合当前按阶段编号的数据。",
      "mean R²容易被某一次测试集方差过小或异常预测拉成很低；仍输出，因为你要求平均R²。",
      "median R²更能代表重复随机验证的一般水平，论文中建议和mean一起报告。",
      "q25 R²代表下四分位稳定性，越高说明模型在不利划分下也能保持预测能力。",
      "stable_score = mean + 0.5*median + 0.25*q25 - 0.25*gap，用于避免只追求偶然高R²。",
      "本版按性状分别设定目标：可溶性糖/木质素0.65，蛋白质/水分0.45，纤维素/硬度0.40；达标优先看median_R2，同时检查mean_R2和q25_R2。"
    ), stringsAsFactors = FALSE
  )
}

# 安全排序：空表或缺列时不报错，直接原样返回
safe_arrange_stable <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(df) == 0 || !("target" %in% names(df)) || !("stable_score" %in% names(df))) return(df)
  df %>% dplyr::arrange(target, dplyr::desc(pass_stable), dplyr::desc(stable_score))
}
safe_arrange_test_r2 <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(df) == 0 || !("target" %in% names(df)) || !("test_r2" %in% names(df))) return(df)
  df %>% dplyr::arrange(target, dplyr::desc(test_r2))
}
empty_eval_norepeat <- function() {
  df <- empty_eval_detail()
  df$repeat_id <- NULL
  df
}

run_all_traits <- function() {
  set.seed(seed)
  obj <- load_bamboo_data(file_path)
  targets <- c("soluble_sugar", "protein", "moisture", "lignin", "cellulose")
  targets <- targets[targets %in% names(obj$dat)]

  all_screen <- empty_eval_summary(); all_refine <- empty_eval_summary(); all_detail <- empty_eval_detail()
  all_errors <- data.frame(target = character(), candidate_id = character(), repeat_id = integer(), error = character(), stringsAsFactors = FALSE)
  best_list <- list(); diag_spxy <- empty_eval_norepeat()

  for (target in targets) {
    message("\n==============================")
    message("开始性状：", target, " / ", trait_cn(target))

    screen_candidates <- make_screen_candidates(target)
    message("初筛候选数：", nrow(screen_candidates), "；重复次数：", screen_repeat)
    screen <- evaluate_candidates(obj, target, screen_candidates, n_repeat = screen_repeat, split_type = stable_split, seed = seed)
    if (nrow(screen$errors) > 0) screen$errors$target <- target
    all_screen <- rbind(all_screen, screen$summary)
    all_errors <- rbind(all_errors, screen$errors[, intersect(names(all_errors), names(screen$errors)), drop = FALSE])

    refine_candidates <- make_refine_candidates(target, screen$summary, keep_n = screen_keep_n)
    message("复筛候选数：", nrow(refine_candidates), "；重复次数：", final_repeat)
    refine <- evaluate_candidates(obj, target, refine_candidates, n_repeat = final_repeat, split_type = stable_split, seed = seed + 999)
    if (nrow(refine$errors) > 0) refine$errors$target <- target
    all_refine <- rbind(all_refine, refine$summary)
    all_detail <- rbind(all_detail, refine$detail)
    all_errors <- rbind(all_errors, refine$errors[, intersect(names(all_errors), names(refine$errors)), drop = FALSE])

    if (nrow(refine$summary) == 0) {
      warning("性状 ", target, " 没有任何候选模型成功运行。请先查看输出Excel中的错误日志。")
      next
    }

    best <- refine$summary %>% dplyr::arrange(dplyr::desc(pass_stable), dplyr::desc(stable_score)) %>% dplyr::slice(1)
    best_list[[target]] <- best

    # SPXY诊断上限：只跑复筛最佳前10个组合的一次SPXY
    top10 <- head(refine$summary, 10)
    spxy_rows <- empty_eval_detail()[, setdiff(names(empty_eval_detail()), "repeat_id"), drop = FALSE]
    for (i in seq_len(nrow(top10))) {
      cand <- top10[i, , drop = FALSE]
      tmp <- tryCatch(run_one_candidate(obj, target, cand, split_type = "spxy", seed = seed), error = function(e) NULL)
      if (!is.null(tmp)) spxy_rows <- rbind(spxy_rows, tmp)
    }
    diag_spxy <- rbind(diag_spxy, spxy_rows)

    message("当前最佳：目标=", round(best$goal_r2, 2), "; 达标=", best$pass_stable,
            "; mean R2=", round(best$mean_r2, 4), "; median R2=", round(best$median_r2, 4), "; q25=", round(best$q25_r2, 4))
  }

  if (length(best_list) == 0) {
    best_summary <- empty_eval_summary()
    best_summary$性状中文 <- character()
  } else {
    best_summary <- dplyr::bind_rows(best_list) %>%
      dplyr::mutate(性状中文 = trait_cn(.data$target)) %>%
      dplyr::select(性状中文, dplyr::everything())
  }

  mean_best <- if (nrow(best_summary) > 0 && "mean_r2" %in% names(best_summary)) mean(best_summary$mean_r2, na.rm = TRUE) else NA_real_
  median_best <- if (nrow(best_summary) > 0 && "median_r2" %in% names(best_summary)) mean(best_summary$median_r2, na.rm = TRUE) else NA_real_
  q25_best <- if (nrow(best_summary) > 0 && "q25_r2" %in% names(best_summary)) mean(best_summary$q25_r2, na.rm = TRUE) else NA_real_
  goal_den <- if (nrow(best_summary) > 0 && "goal_r2" %in% names(best_summary)) sum(!is.na(best_summary$goal_r2)) else 0
  pass_n <- if (nrow(best_summary) > 0 && "pass_stable" %in% names(best_summary)) sum(best_summary$pass_stable, na.rm = TRUE) else 0
  avg_row <- data.frame(
    项目 = c("六性状最佳mean R²平均", "六性状最佳median R²平均", "六性状最佳q25 R²平均", "重点性状稳定达标数量", "分性状达标标准"),
    数值 = c(mean_best, median_best, q25_best, paste0(pass_n, "/", goal_den), "可溶性糖/木质素>=0.75；蛋白质/水分>=0.65；纤维素/硬度不设硬性目标"),
    目标 = c(NA, NA, NA, "只统计goal_r2非NA的四个性状；median达标且mean不低于目标0.05且q25不低于目标0.15", "分性状阈值"), stringsAsFactors = FALSE
  )

  wb <- openxlsx::createWorkbook()
  add_sheet <- function(name, df) {
    openxlsx::addWorksheet(wb, name)
    openxlsx::writeData(wb, name, as.data.frame(df))
    if (ncol(df) > 0) {
      openxlsx::freezePane(wb, name, firstRow = TRUE)
      openxlsx::setColWidths(wb, name, cols = 1:ncol(df), widths = "auto")
    }
  }
  add_sheet("六性状最佳稳定模型", best_summary)
  add_sheet("分性状目标汇总", avg_row)
  add_sheet("复筛全部稳定性", safe_arrange_stable(all_refine))
  add_sheet("初筛全部稳定性", safe_arrange_stable(all_screen))
  add_sheet("复筛重复明细", all_detail)
  add_sheet("SPXY诊断上限", safe_arrange_test_r2(diag_spxy))
  add_sheet("错误日志", all_errors)
  add_sheet("方法说明", method_notes())

  out_file <- file.path(out_dir, "00_六性状_分性状目标达标搜索结果.xlsx")
  openxlsx::saveWorkbook(wb, out_file, overwrite = TRUE)
  message("\n完成。总结果：", out_file)
  # 导出最佳模型权重（审稿复现用）
  tryCatch({
    message("\n--- 导出最佳模型权重 ---")
    export_best_models(best_summary, obj, out_dir, all_detail, all_refine)
  }, error = function(e) message("模型导出跳过：", e$message))
  invisible(list(best = best_summary, avg = avg_row, screen = all_screen, refine = all_refine, detail = all_detail, spxy = diag_spxy, errors = all_errors))
}

result <- run_all_traits()
