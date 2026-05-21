ECEA <- function(
    # --- 兼容性输入 ---
  gene_set_ind = NULL,
  gene_expression_matrix = NULL,
  num_diff_samples = NULL,

  # --- 高效输入 (用于 run_ECEA 包装器) ---
  gene_set = NULL,
  precomputed_ranks = NULL,
  pre_sorted_data = NULL
) {

  # ============================================================================
  # 内部数据准备：兼容所有调用方式
  # ============================================================================

  if (!is.null(pre_sorted_data) && !is.null(gene_set)) {
    # ---- 路径 A (最高效) ----
    sorted_data <- pre_sorted_data
    sorted_data$is_in_set <- sorted_data$gene_id %in% gene_set

  } else {
    # ---- 路径 B/C ----
    if (!is.null(gene_set_ind)) {
      if (is.null(gene_expression_matrix)) stop("错误: 当使用 'gene_set_ind' 时, 必须提供 'gene_expression_matrix'")
      final_gene_set <- rownames(gene_expression_matrix)[gene_set_ind]
    } else if (!is.null(gene_set)) {
      final_gene_set <- gene_set
    } else {
      stop("错误: 必须提供 'gene_set' 或 'gene_set_ind'")
    }

    if (is.null(precomputed_ranks)) {
      # ---- 步骤 1: 对每个基因进行t检验，并作为排序指标 ----
      if (is.null(gene_expression_matrix)) stop("错误: 当需要计算t-test时, 必须提供 'gene_expression_matrix'")
      t_test_results_res <- matrixTests::row_t_welch(
        gene_expression_matrix[, 1:num_diff_samples],
        gene_expression_matrix[, (num_diff_samples + 1):ncol(gene_expression_matrix)]
      )
      t_test_results <- t_test_results_res$statistic
      names(t_test_results) <- rownames(gene_expression_matrix)
    } else {
      t_test_results <- precomputed_ranks
    }

    # ---- 步骤 2: 清理排序指标中的无效值 (NaN, Inf) ----
    if (any(is.nan(t_test_results))) { t_test_results[is.nan(t_test_results)] <- 0 }
    if (any(is.infinite(t_test_results))) {
      finite_vals <- t_test_results[is.finite(t_test_results)]
      max_finite_val <- if(length(finite_vals) > 0) max(finite_vals, na.rm = TRUE) else 1
      min_finite_val <- if(length(finite_vals) > 0) min(finite_vals, na.rm = TRUE) else -1
      t_test_results[t_test_results == Inf] <- max_finite_val + 1; t_test_results[t_test_results == -Inf] <- min_finite_val - 1
    }

    # ---- 步骤 3: 构建包含基因信息和排序指标的数据框 ----
    if(is.null(names(t_test_results))) stop("排序指标 (ranks) 必须是一个命名的向量 (基因ID作为名字).")
    gene_ids <- names(t_test_results)
    is_in_set <- gene_ids %in% final_gene_set
    gene_data_df <- data.frame(gene_id = gene_ids, rank_metric = t_test_results, is_in_set = is_in_set)

    # ---- 步骤 4: 统一按指标降序排列基因列表 ----
    sorted_data <- gene_data_df[order(gene_data_df$rank_metric, decreasing = TRUE), ]
  }

  # ========================================================================
  # GSEA 核心计算
  # ========================================================================

  # ---- 步骤 5: 计算GSEA随机游走路径 ----
  N <- nrow(sorted_data)
  N_H <- sum(sorted_data$is_in_set)

  # 边缘情况处理：如果基因集为空或包含所有基因，则无富集意义
  if (N_H == 0 || N_H == N) {
    return(list(ES = 0, tau = NA, W = 0, pvalue = 1.0, direction = "NEUTRAL", rs_values = if(N > 0) rep(0, N) else numeric(0), sorted_data = sorted_data))
  }

  val_hit <- sqrt((N - N_H) / N_H); val_miss <- -sqrt(N_H / (N - N_H))
  rs_values <- cumsum(ifelse(sorted_data$is_in_set, val_hit, val_miss))

  # 步骤 6: 同时找到路径的最大值(波峰)和最小值(波谷)
  max_es <- max(rs_values); min_es <- min(rs_values)
  max_tau <- which.max(rs_values); min_tau <- which.min(rs_values)

  # 步骤 7: 比较绝对值，确定最终富集得分(ES)、方向(direction)和位置(tau)
  if (max_es > abs(min_es)) {
    ES <- max_es; direction <- "UP"; tau <- max_tau
  } else {
    ES <- min_es; direction <- "DOWN"; tau <- min_tau
  }

  M_effective_for_W <- abs(ES)

  if (M_effective_for_W <= 1e-9) {
    return(list(ES = 0, tau = tau, W = 0, pvalue = 1.0, direction = "NEUTRAL", rs_values = rs_values, sorted_data = sorted_data))
  }

  # ---- 步骤 8: 计算近似卡方统计量 W 并得出 p值 ----
  denominator <- tau * (N - tau)
  denominator_max <- max_tau * (N - max_tau)
  denominator_min <- min_tau * (N - min_tau)
  W_val <- if(denominator == 0) 0 else { (N * M_effective_for_W^2) / denominator }
  W_max <- (N * max_es^2) / denominator_max
  W_min <- (N * min_es^2) / denominator_min
  pvalue <- pchisq(W_val*qchisq(0.95, df = 3)/qchisq(0.975, df = 3), df = 3, lower.tail = FALSE)
  log_pvalue <- pchisq(W_val, df = 3, lower.tail = FALSE, log.p = TRUE)

  # ---- 步骤 9: 返回包含所有关键信息的完整结果 ----
  return(list(
    ES = ES,
    pvalue = pvalue,
    log_pvalue = log_pvalue,
    direction = direction,
    tau = tau,
    W = W_val,
    rs_values = rs_values,
    sorted_data = sorted_data,
    max_es = max_es,
    min_es = min_es,
    max_tau = max_tau,
    min_tau = min_tau,
    W_max = W_max,
    W_min = W_min
  ))
}

ECEA_VIF <- function(
    gene_set_ind = NULL,
    gene_expression_matrix = NULL,
    num_diff_samples = NULL,
    gene_set = NULL,
    precomputed_ranks = NULL,
    pre_sorted_data = NULL,
    use_correction = TRUE
) {
  if (!is.null(pre_sorted_data) && !is.null(gene_set)) {
    sorted_data <- pre_sorted_data
    sorted_data$is_in_set <- sorted_data$gene_id %in% gene_set
    final_gene_set <- gene_set
  } else {
    if (!is.null(gene_set_ind)) {
      if (is.null(gene_expression_matrix)) stop("错误: 当使用 'gene_set_ind' 时, 必须提供 'gene_expression_matrix'")
      final_gene_set <- rownames(gene_expression_matrix)[gene_set_ind]
    } else if (!is.null(gene_set)) {
      final_gene_set <- gene_set
    } else {
      stop("错误: 必须提供 'gene_set' 或 'gene_set_ind'")
    }
    if (is.null(precomputed_ranks)) {
      if (is.null(gene_expression_matrix)) stop("错误: 当需要计算 t-test 时, 必须提供 'gene_expression_matrix'")
      t_test_results_res <- matrixTests::row_t_welch(
        gene_expression_matrix[, 1:num_diff_samples],
        gene_expression_matrix[, (num_diff_samples + 1):ncol(gene_expression_matrix)]
      )
      t_test_results <- t_test_results_res$statistic
      names(t_test_results) <- rownames(gene_expression_matrix)
    } else {
      t_test_results <- precomputed_ranks
    }
    if (any(is.nan(t_test_results))) { t_test_results[is.nan(t_test_results)] <- 0 }
    if (any(is.infinite(t_test_results))) {
      finite_vals <- t_test_results[is.finite(t_test_results)]
      max_finite_val <- if(length(finite_vals) > 0) max(finite_vals, na.rm = TRUE) else 1
      min_finite_val <- if(length(finite_vals) > 0) min(finite_vals, na.rm = TRUE) else -1
      t_test_results[t_test_results == Inf] <- max_finite_val + 1
      t_test_results[t_test_results == -Inf] <- min_finite_val - 1
    }
    if(is.null(names(t_test_results))) stop("排序指标 (ranks) 必须是一个命名的向量 (基因ID作为名字).")
    gene_ids <- names(t_test_results)
    is_in_set <- gene_ids %in% final_gene_set
    gene_data_df <- data.frame(gene_id = gene_ids, rank_metric = t_test_results, is_in_set = is_in_set)
    sorted_data <- gene_data_df[order(gene_data_df$rank_metric, decreasing = TRUE), ]
  }
  N <- nrow(sorted_data)
  N_H <- sum(sorted_data$is_in_set)
  mean_rho <- 0
  VIF <- 1.0
  if (N_H == 0 || N_H == N) {
    return(list(ES = 0, pvalue = 1.0, log_pvalue = 0, direction = "NEUTRAL", tau = NA, W = 0, W_adjusted = 0, VIF = 1.0, mean_rho = 0, N_H = N_H, rs_values = if(N > 0) rep(0, N) else numeric(0), sorted_data = sorted_data, max_es = 0, min_es = 0, max_tau = NA, min_tau = NA, W_max = 0, W_min = 0))
  }
  val_hit <- sqrt((N - N_H) / N_H)
  val_miss <- -sqrt(N_H / (N - N_H))
  rs_values <- cumsum(ifelse(sorted_data$is_in_set, val_hit, val_miss))
  max_es <- max(rs_values); min_es <- min(rs_values)
  max_tau <- which.max(rs_values); min_tau <- which.min(rs_values)
  if (max_es > abs(min_es)) { ES <- max_es; direction <- "UP"; tau <- max_tau
  } else { ES <- min_es; direction <- "DOWN"; tau <- min_tau }
  M_effective_for_W <- abs(ES)
  if (M_effective_for_W <= 1e-9) {
    return(list(ES = 0, pvalue = 1.0, log_pvalue = 0, direction = "NEUTRAL", tau = tau, W = 0, W_adjusted = 0, VIF = 1.0, mean_rho = 0, N_H = N_H, rs_values = rs_values, sorted_data = sorted_data, max_es = max_es, min_es = min_es, max_tau = max_tau, min_tau = min_tau, W_max = 0, W_min = 0))
  }
  denominator <- tau * (N - tau); denominator_max <- max_tau * (N - max_tau); denominator_min <- min_tau * (N - min_tau)
  W_val <- if(denominator == 0) 0 else { (N * M_effective_for_W^2) / denominator }
  W_max <- if(denominator_max == 0) 0 else { (N * max_es^2) / denominator_max }
  W_min <- if(denominator_min == 0) 0 else { (N * min_es^2) / denominator_min }

  # ========================================================================
  # VIF 校正
  # ========================================================================
  if (use_correction && N_H > 1 && !is.null(gene_expression_matrix)) {
    current_set_indices <- match(final_gene_set, rownames(gene_expression_matrix))
    current_set_indices <- current_set_indices[!is.na(current_set_indices)]

    if(length(current_set_indices) > 1) {
      sub_expr <- gene_expression_matrix[current_set_indices, , drop=FALSE]

      # 子采样控制计算复杂度
      if(nrow(sub_expr) > 100) {
        sub_expr <- sub_expr[sample(1:nrow(sub_expr), 100), ]
      }

      # 计算 pathway 内相关性
      cor_mat <- cor(t(sub_expr))
      rho_vals <- cor_mat[upper.tri(cor_mat)]

      # 恢复整体平均值计算，允许正负噪音相互抵消
      mean_rho <- mean(rho_vals, na.rm = TRUE)

      # 只有当整体表现为正相关（存在共表达结构）时，才触发惩罚
      mean_rho <- max(0, mean_rho)


      VIF <- 1 + (N_H - 1) * mean_rho
    }
  }


  W_adjusted <- W_val / (VIF ^ 0.8)
  if (W_adjusted < 0) W_adjusted <- 0

  pvalue <- pchisq(W_adjusted * qchisq(0.95, df = 3) / qchisq(0.975, df = 3), df = 3, lower.tail = FALSE)
  log_pvalue <- pchisq(W_adjusted, df = 3, lower.tail = FALSE, log.p = TRUE)
  return(list(ES = ES, pvalue = pvalue, log_pvalue = log_pvalue, direction = direction, tau = tau, W = W_val, W_adjusted = W_adjusted, VIF = VIF, mean_rho = mean_rho, N_H = N_H, rs_values = rs_values, sorted_data = sorted_data, max_es = max_es, min_es = min_es, max_tau = max_tau, min_tau = min_tau, W_max = W_max, W_min = W_min))
}
