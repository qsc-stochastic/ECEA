WDECEA <- function(
    # --- 兼容性输入 ---
  gene_set_ind = NULL,
  gene_expression_matrix = NULL,
  num_diff_samples = NULL,

  # --- 高效输入 (用于 run_ECEA 包装器) ---
  gene_set = NULL,
  precomputed_ranks = NULL,
  pre_sorted_data = NULL,

  # --- 新增：统计学校正参数 ---
  use_correlation_correction = TRUE  # 是否开启基于长期方差的校正
) {

  # ============================================================================
  # 内部数据准备：兼容所有调用方式 (保持原样)
  # ============================================================================

  if (!is.null(pre_sorted_data) && !is.null(gene_set)) {
    # ---- 路径 A (最高效) ----
    sorted_data <- pre_sorted_data
    sorted_data$is_in_set <- sorted_data$gene_id %in% gene_set

  } else {
    # ---- 路径 B/C (常规处理) ----
    if (!is.null(gene_set_ind)) {
      if (is.null(gene_expression_matrix)) stop("错误: 当使用 'gene_set_ind' 时, 必须提供 'gene_expression_matrix'")
      final_gene_set <- rownames(gene_expression_matrix)[gene_set_ind]
    } else if (!is.null(gene_set)) {
      final_gene_set <- gene_set
    } else {
      stop("错误: 必须提供 'gene_set' 或 'gene_set_ind'")
    }

    if (is.null(precomputed_ranks)) {
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

    # 清理无效值
    if (any(is.nan(t_test_results))) { t_test_results[is.nan(t_test_results)] <- 0 }
    if (any(is.infinite(t_test_results))) {
      finite_vals <- t_test_results[is.finite(t_test_results)]
      max_finite_val <- if(length(finite_vals) > 0) max(finite_vals, na.rm = TRUE) else 1
      min_finite_val <- if(length(finite_vals) > 0) min(finite_vals, na.rm = TRUE) else -1
      t_test_results[t_test_results == Inf] <- max_finite_val + 1; t_test_results[t_test_results == -Inf] <- min_finite_val - 1
    }

    if(is.null(names(t_test_results))) stop("排序指标 (ranks) 必须是一个命名的向量.")
    gene_ids <- names(t_test_results)
    is_in_set <- gene_ids %in% final_gene_set
    gene_data_df <- data.frame(gene_id = gene_ids, rank_metric = t_test_results, is_in_set = is_in_set)

    sorted_data <- gene_data_df[order(gene_data_df$rank_metric, decreasing = TRUE), ]
  }

  # ========================================================================
  # GSEA 核心计算 & 理论修正
  # ========================================================================

  N <- nrow(sorted_data)
  N_H <- sum(sorted_data$is_in_set)

  if (N_H == 0 || N_H == N) {
    return(list(ES = 0, tau = NA, W = 0, pvalue = 1.0, direction = "NEUTRAL",
                rs_values = if(N > 0) rep(0, N) else numeric(0), sigma_sq = 1))
  }

  # 定义每一步的增量 (Increments)
  val_hit <- sqrt((N - N_H) / N_H)
  val_miss <- -sqrt(N_H / (N - N_H))
  increments <- ifelse(sorted_data$is_in_set, val_hit, val_miss)

  # 计算累积和 (随机游走路径)
  rs_values <- cumsum(increments)

  # ---- 新增核心逻辑：长期方差 (Long-run Variance) 估计 ----
  # 理论基础：sigma^2 = gamma(0) + 2 * sum(gamma(k))
  # 目的：校正基因排序中的局部相关性

  sigma_sq <- 1.0 # 默认值 (独立假设)

  if (use_correlation_correction && N > 10) {
    # 1. 确定滞后阶数 (Lag truncation parameter)
    # 经验法则：m 约等于 N^(1/4) 或 N^(1/3)。这里取 N^(1/3) 以捕获稍微长一点的相关性
    max_lag <- floor(N^(1/3))

    # 2. 计算自协方差 (Autocovariance)
    # type="covariance" 返回未经标准化的协方差 gamma(k)
    acf_res <- acf(increments, lag.max = max_lag, plot = FALSE, type = "covariance")
    gammas <- acf_res$acf # gammas[1] 是 gamma(0), gammas[k+1] 是 gamma(k)

    # 3. 应用 Bartlett Kernel (Newey-West estimator) 计算加权和
    # 这种加权方式 (1 - k/(m+1)) 保证了估计出的 sigma_sq 必定为正，避免负方差的数学尴尬
    weights <- 1 - (1:max_lag) / (max_lag + 1)

    # 公式: gamma(0) + 2 * sum( weight_k * gamma(k) )
    gamma_0 <- gammas[1]
    weighted_sum_cov <- sum(weights * gammas[2:(max_lag + 1)])

    estimated_sigma_sq <- gamma_0 + 2 * weighted_sum_cov

    # 4. 边界保护：防止过度校正或数值错误
    # 如果数据极其反常导致负值(虽然Bartlett极少出现)或接近0，强制设为下限
    if (!is.na(estimated_sigma_sq) && estimated_sigma_sq > 1e-4) {
      sigma_sq <- estimated_sigma_sq
    }
  }

  # ========================================================================
  # 统计量计算
  # ========================================================================

  max_es <- max(rs_values); min_es <- min(rs_values)
  max_tau <- which.max(rs_values); min_tau <- which.min(rs_values)

  if (max_es > abs(min_es)) {
    ES <- max_es; direction <- "UP"; tau <- max_tau
  } else {
    ES <- min_es; direction <- "DOWN"; tau <- min_tau
  }

  M_effective_for_W <- abs(ES)

  if (M_effective_for_W <= 1e-9) {
    return(list(ES = 0, tau = tau, W = 0, pvalue = 1.0, direction = "NEUTRAL", sigma_sq = sigma_sq))
  }

  # ---- 步骤 8 (修正): 计算 W 并应用 sigma_sq 校正 ----
  denominator <- tau * (N - tau)
  denominator_max <- max_tau * (N - max_tau)
  denominator_min <- min_tau * (N - min_tau)

  # 原始 W (假设独立)
  W_raw <- if(denominator == 0) 0 else { (N * M_effective_for_W^2) / denominator }
  W_max_raw <- (N * max_es^2) / denominator_max
  W_min_raw <- (N * min_es^2) / denominator_min

  # 【关键修正】: 将 W 除以长期方差 sigma^2
  # 如果序列正相关 (sigma_sq > 1)，W 会变小，P值变大 (更保守)
  # 如果序列负相关 (sigma_sq < 1)，W 会变大，P值变小 (更敏感)
  W_val <- W_raw / sigma_sq
  W_max <- W_max_raw / sigma_sq
  W_min <- W_min_raw / sigma_sq

  # 计算 P 值 (自由度 df=3 的卡方近似)
  # 注意：这里我们假设校正后的 W_val 恢复到了标准布朗桥的尺度
  pvalue <- pchisq(W_val*qchisq(0.95, df = 3)/qchisq(0.975, df = 3), df = 3, lower.tail = FALSE)
  log_pvalue <- pchisq(W_val, df = 3, lower.tail = FALSE, log.p = TRUE)

  # ---- 步骤 9: 返回结果 (包含校正因子 sigma_sq 供诊断) ----
  return(list(
    ES = ES,
    pvalue = pvalue,
    log_pvalue = log_pvalue,
    direction = direction,
    tau = tau,
    W = W_val,         # 校正后的统计量
    W_raw = W_raw,     # 原始统计量 (用于对比)
    sigma_sq = sigma_sq, # 长期方差校正因子 (重要诊断指标)
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

weighted_ECEA <- function(
    # --- 兼容性输入 (与原函数保持一致) ---
  gene_set_ind = NULL,
  gene_expression_matrix = NULL,
  num_diff_samples = NULL,

  # --- 高效输入 (用于 run_ECEA 包装器) ---
  gene_set = NULL,
  precomputed_ranks = NULL,
  pre_sorted_data = NULL,

  # --- 新增参数: 权重指数 ---
  weight_exponent = 1.0  # 1.0 = 标准GSEA加权; 0 = 经典KS统计量(无权)
) {

  # ============================================================================
  # 1. 内部数据准备 (保持原逻辑不变)
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
      if (is.null(gene_expression_matrix)) stop("错误: 当需要计算t-test时, 必须提供 'gene_expression_matrix'")
      # 使用 matrixTests 快速计算
      t_test_results_res <- matrixTests::row_t_welch(
        gene_expression_matrix[, 1:num_diff_samples],
        gene_expression_matrix[, (num_diff_samples + 1):ncol(gene_expression_matrix)]
      )
      t_test_results <- t_test_results_res$statistic
      names(t_test_results) <- rownames(gene_expression_matrix)
    } else {
      t_test_results <- precomputed_ranks
    }

    # 清理无效值
    if (any(is.nan(t_test_results))) { t_test_results[is.nan(t_test_results)] <- 0 }
    if (any(is.infinite(t_test_results))) {
      finite_vals <- t_test_results[is.finite(t_test_results)]
      max_finite_val <- if(length(finite_vals) > 0) max(finite_vals, na.rm = TRUE) else 1
      min_finite_val <- if(length(finite_vals) > 0) min(finite_vals, na.rm = TRUE) else -1
      t_test_results[t_test_results == Inf] <- max_finite_val + 1
      t_test_results[t_test_results == -Inf] <- min_finite_val - 1
    }

    # 构建并排序
    if(is.null(names(t_test_results))) stop("排序指标必须命名")
    gene_ids <- names(t_test_results)
    is_in_set <- gene_ids %in% final_gene_set
    gene_data_df <- data.frame(gene_id = gene_ids, rank_metric = t_test_results, is_in_set = is_in_set)
    sorted_data <- gene_data_df[order(gene_data_df$rank_metric, decreasing = TRUE), ]
  }

  # ============================================================================
  # 2. GSEA 核心计算 (已针对加权和方差修正进行升级)
  # ============================================================================

  N <- nrow(sorted_data)
  N_H <- sum(sorted_data$is_in_set)

  # 边缘情况：空集或全集
  if (N_H == 0 || N_H == N) {
    return(list(ES = 0, tau = NA, W = 0, pvalue = 1.0, direction = "NEUTRAL",
                rs_values = if(N > 0) rep(0, N) else numeric(0), sorted_data = sorted_data))
  }

  # ---- 步骤 A: 计算加权步长 (Weighted Steps) ----
  # 这是一个关键改变。如果是 Hit，步长与 |metric|^p 成正比。
  # 如果是 Miss，步长保持均匀 (标准GSEA做法)，或者也可以加权 (视具体定义)。
  # 这里采用标准 GSEA (Subramanian et al.) 逻辑：Hit 加权，Miss 均匀。

  abs_metric <- abs(sorted_data$rank_metric)
  # 应用权重指数
  adj_weights <- abs_metric ^ weight_exponent

  # 计算 Hit 的总权重，用于归一化
  sum_weights_hit <- sum(adj_weights[sorted_data$is_in_set])

  # 防止分母为0 (极少数情况)
  if (sum_weights_hit == 0) sum_weights_hit <- 1

  # 初始化步长向量
  steps <- numeric(N)

  # Hit 的步长：权重 / 总权重 (正向贡献)
  steps[sorted_data$is_in_set] <- adj_weights[sorted_data$is_in_set] / sum_weights_hit

  # Miss 的步长：-1 / (N - N_H) (负向惩罚，确保总和为0)
  # 注意：Miss 通常不加权，因为它们代表背景噪音
  steps[!sorted_data$is_in_set] <- -1 / (N - N_H)

  # ---- 步骤 B: 计算随机游走路径 (Running Sum) ----
  rs_values <- cumsum(steps)

  # 找到波峰和波谷
  max_es <- max(rs_values)
  min_es <- min(rs_values)
  max_tau_idx <- which.max(rs_values) # 物理位置索引
  min_tau_idx <- which.min(rs_values)

  # 确定 ES, 方向, 和对应的物理位置索引
  if (max_es > abs(min_es)) {
    ES <- max_es; direction <- "UP"; tau_idx <- max_tau_idx
  } else {
    ES <- min_es; direction <- "DOWN"; tau_idx <- min_tau_idx
  }

  # ---- 步骤 C: 基于"内在时间"(Intrinsic Time) 计算统计量 ----
  # 这是解决"加权导致P值失效"的核心步骤。
  # 我们需要计算在 tau_idx 这一刻，方差累积了多少。

  # 1. 计算每一步的方差贡献 (即步长的平方)
  step_variances <- steps^2

  # 2. 计算累积方差向量
  cum_variances <- cumsum(step_variances)
  total_variance <- cum_variances[N] # 总方差 V

  # 3. 找到峰值位置对应的归一化累积方差 (即 G(tau))
  # 这就是我们在分布图中对应的 "x轴位置"
  g_star <- cum_variances[tau_idx] / total_variance

  # ---- 步骤 D: 计算 W 统计量和 P 值 ----

  # 边缘处理：如果峰值出现在方差累积的极早期或极晚期 (例如前1%或后1%)
  # 分母会趋近于0，导致 W 爆炸。通常这时认为是噪音。
  if (g_star < 0.001 || g_star > 0.999) {
    W_val <- 0; pvalue <- 1.0
  } else {
    # 核心公式: W = ES^2 / (V * g * (1-g))
    # 这个统计量服从 Chi-sq (df=3) 的分布规律 (基于布朗桥最大值理论)
    W_val <- (ES^2) / (total_variance * g_star * (1 - g_star))

    # 计算近似 P 值
    pvalue <- pchisq(W_val, df = 3, lower.tail = FALSE)
  }

  # Log P-value 用于高精度需求
  log_pvalue <- if (W_val > 0) pchisq(W_val, df = 3, lower.tail = FALSE, log.p = TRUE) else 0

  # ---- 步骤 E: 辅助统计量 (可选) ----
  # 为了兼容性，也可以计算一下 max/min 的 W 值
  # 注意：这里需要对所有位置计算 g_vec
  g_vec <- cum_variances / total_variance
  # 避免除以0
  valid_mask <- (g_vec > 0.001 & g_vec < 0.999)

  W_max <- 0; W_min <- 0
  if (any(valid_mask)) {
    # 仅在有效范围内估算潜在的最大W
    denom_vec <- total_variance * g_vec * (1 - g_vec)
    # 这里的计算仅作参考，不影响主 ES 的 P 值
    # W_vec <- (rs_values^2) / denom_vec
  }

  # ============================================================================
  # 3. 返回结果
  # ============================================================================
  return(list(
    ES = ES,
    pvalue = pvalue,
    log_pvalue = log_pvalue,
    direction = direction,
    tau = tau_idx,             # 物理位置 (第几个基因)
    g_star = g_star,           # 内在时间位置 (累积方差占比，0~1)
    W = W_val,
    rs_values = rs_values,
    sorted_data = sorted_data,
    max_es = max_es,
    min_es = min_es,
    total_variance = total_variance # 返回总方差供参考
  ))
}

ECEA_perm <- function(
    # --- 兼容性输入 ---
  gene_set_ind = NULL,
  gene_expression_matrix = NULL,
  num_diff_samples = NULL,

  # --- 高效输入 (用于 run_ECEA 包装器) ---
  gene_set = NULL,
  precomputed_ranks = NULL,
  pre_sorted_data = NULL,

  # --- 新增：置换检验参数 ---
  n_perm = 1000,
  seed = 123,
  epsilon = 1e-9 # 防止分母为0的微小量
) {

  # 设置随机种子以保证结果可重复
  if (!is.null(seed)) set.seed(seed)

  # ============================================================================
  # 1. 内部辅助函数：计算核心统计量 W
  # ============================================================================
  # 该函数接收一个逻辑向量（表示排序后的基因是否在集合中），返回 W 统计量及其他详情
  calc_GSEA_stat <- function(is_in_set_vec) {
    N <- length(is_in_set_vec)
    N_H <- sum(is_in_set_vec)

    if (N_H == 0 | N_H == N) {
      return(list(W = 0, ES = 0, tau = NA, direction = "NEUTRAL", rs_values = rep(0, N)))
    }

    val_hit <- sqrt((N - N_H) / N_H)
    val_miss <- -sqrt(N_H / (N - N_H))

    # 随机游走路径
    rs_values <- cumsum(ifelse(is_in_set_vec, val_hit, val_miss))

    # 寻找波峰波谷
    max_es <- max(rs_values)
    min_es <- min(rs_values)
    max_tau <- which.max(rs_values)
    min_tau <- which.min(rs_values)

    if (max_es > abs(min_es)) {
      ES <- max_es; direction <- "UP"; tau <- max_tau
    } else {
      ES <- min_es; direction <- "DOWN"; tau <- min_tau
    }

    M_effective <- abs(ES)

    if (M_effective <= epsilon) {
      return(list(W = 0, ES = 0, tau = tau, direction = "NEUTRAL", rs_values = rs_values))
    }

    # 计算 W 统计量 (M^2 / Var_bridge)
    # 引入 epsilon 防止 tau 在边界时分母为 0 导致数值爆炸
    denominator <- tau * (N - tau)
    W_val <- if (denominator == 0) 0 else { (N * M_effective^2) / denominator }

    return(list(
      W = W_val,
      ES = ES,
      tau = tau,
      direction = direction,
      rs_values = rs_values,
      max_es = max_es, min_es = min_es,
      max_tau = max_tau, min_tau = min_tau
    ))
  }

  # ============================================================================
  # 2. 数据准备与观测值计算
  # ============================================================================

  # --- 路径 A: 使用预排序数据 (只能做基因集置换) ---
  if (!is.null(pre_sorted_data) &&!is.null(gene_set)) {
    sorted_data <- pre_sorted_data
    sorted_data$is_in_set <- sorted_data$gene_id %in% gene_set
    perm_mode <- "gene_set"

  } else {
    # --- 路径 B: 使用原始表达矩阵 (推荐：可以做表型置换) ---
    perm_mode <- "phenotype"

    if (!is.null(gene_set_ind)) {
      if (is.null(gene_expression_matrix)) stop("错误: 使用 'gene_set_ind' 需提供表达矩阵")
      final_gene_set <- rownames(gene_expression_matrix)[gene_set_ind]
    } else if (!is.null(gene_set)) {
      final_gene_set <- gene_set
    } else {
      stop("错误: 必须提供 'gene_set' 或 'gene_set_ind'")
    }

    # 计算观测值的 t-test 排序
    if (is.null(precomputed_ranks)) {
      if (is.null(gene_expression_matrix)) stop("错误: 计算 t-test 需要表达矩阵")
      obs_t_res <- matrixTests::row_t_welch(
        gene_expression_matrix[, 1:num_diff_samples],
        gene_expression_matrix[, (num_diff_samples + 1):ncol(gene_expression_matrix)]
      )
      t_test_results <- obs_t_res$statistic
      names(t_test_results) <- rownames(gene_expression_matrix)
    } else {
      t_test_results <- precomputed_ranks
    }

    # 清理无效值
    t_test_results[is.nan(t_test_results)] <- 0
    if (any(is.infinite(t_test_results))) {
      finite_vals <- t_test_results[is.finite(t_test_results)]
      max_val <- if(length(finite_vals)) max(finite_vals) else 1
      min_val <- if(length(finite_vals)) min(finite_vals) else -1
      t_test_results[t_test_results == Inf] <- max_val + 1
      t_test_results[t_test_results == -Inf] <- min_val - 1
    }

    # 构建排序列表
    gene_ids <- names(t_test_results)
    is_in_set <- gene_ids %in% final_gene_set
    gene_data_df <- data.frame(gene_id = gene_ids, rank_metric = t_test_results, is_in_set = is_in_set)
    sorted_data <- gene_data_df
  }

  # --- 计算观测到的统计量 ---
  obs_res <- calc_GSEA_stat(sorted_data$is_in_set)
  W_obs <- obs_res$W

  # ============================================================================
  # 3. 置换检验 (Permutation Test)
  # ============================================================================

  perm_Ws <- numeric(n_perm)

  if (perm_mode == "phenotype" &&!is.null(gene_expression_matrix) && is.null(precomputed_ranks)) {
    # 【表型置换】：打乱样本标签 -> 重新计算t值 -> 重新排序 -> 计算W
    # 这是处理基因相关性的正确方法

    total_samples <- ncol(gene_expression_matrix)
    all_indices <- 1:total_samples

    for (i in 1:n_perm) {
      # 1. 随机打乱样本索引
      shuffled_indices <- sample(all_indices)

      # 2. 重新划分组别
      # 注意：为了速度，这里不完全重新构建矩阵，而是直接通过索引传给 row_t_welch
      # matrixTests 允许直接传矩阵，我们这里手动切分索引可能会慢，
      # 但这是为了模拟“样本标签打乱”。

      # 模拟：Group A 是打乱后的前 num_diff_samples 个，Group B 是剩下的
      # 实际上 matrixTests::row_t_welch(x, y)
      x_idx <- shuffled_indices[1:num_diff_samples]
      y_idx <- shuffled_indices[(num_diff_samples + 1):total_samples]

      # 3. 快速计算 t 统计量
      # 注意：如果矩阵很大，这一步在R中循环1000次可能会慢。
      # 生产环境中通常会用 C++ (Rcpp) 或矩阵乘法优化这一步。
      perm_res <- matrixTests::row_t_welch(
        gene_expression_matrix,
        gene_expression_matrix
      )
      perm_stats <- perm_res$statistic

      # 处理 Inf/NaN (简化版，为了速度)
      perm_stats[is.na(perm_stats)] <- 0

      # 4. 重新获取 is_in_set 的排序
      # 我们不需要完整的dataframe，只需要知道根据 perm_stats 排序后，is_in_set 变成了什么样
      # 利用 order 获取排序索引
      rank_idx <- order(perm_stats, decreasing = TRUE)

      # 原始的 is_in_set 向量 (对应行号顺序)
      # 也就是上面构建 gene_data_df 之前的 is_in_set
      # 注意：sorted_data 是排过序的，不能直接用。我们需要原始顺序的 logical 向量。
      raw_is_in_set <- if(exists("gene_data_df")) gene_data_df$is_in_set else (rownames(gene_expression_matrix) %in% final_gene_set)

      perm_sorted_is_in_set <- raw_is_in_set[rank_idx]

      # 5. 计算 W
      perm_Ws[i] <- calc_GSEA_stat(perm_sorted_is_in_set)$W
    }

  } else {
    # 【基因集置换】：仅打乱排序列表中的基因标签
    # 警告：这种方法比表型置换快，但不能校正基因间相关性，会导致 P 值偏小（假阳性）。
    # 当没有表达矩阵时（路径 A），只能用这个。

    current_bool_vec <- sorted_data$is_in_set
    for (i in 1:n_perm) {
      # 简单地打乱 TRUE/FALSE 的位置
      perm_bool_vec <- sample(current_bool_vec)
      perm_Ws[i] <- calc_GSEA_stat(perm_bool_vec)$W
    }
  }

  # ============================================================================
  # 4. 计算 P 值与结果汇总
  # ============================================================================

  # P值计算公式：(超过观测值的次数 + 1) / (置换总次数 + 1)
  # 这里的 +1 是为了防止 P值为0 (pseudocount)
  pvalue <- (sum(perm_Ws >= W_obs) + 1) / (n_perm + 1)

  # 返回结果
  return(list(
    ES = obs_res$ES,
    pvalue = pvalue,
    log_pvalue = log(pvalue), # 自然对数
    W = W_obs,
    perm_Ws = perm_Ws, # 返回置换分布供画图用
    perm_mode = perm_mode,
    direction = obs_res$direction,
    tau = obs_res$tau,
    rs_values = obs_res$rs_values,
    sorted_data = sorted_data,
    max_es = obs_res$max_es,
    min_es = obs_res$min_es,
    max_tau = obs_res$max_tau,
    min_tau = obs_res$min_tau
  ))
}

ECEA_weighted_perm <- function(
    # --- 兼容性输入 ---
  gene_set_ind = NULL,
  gene_expression_matrix = NULL,
  num_diff_samples = NULL,

  # --- 高效输入 ---
  gene_set = NULL,
  precomputed_ranks = NULL,
  pre_sorted_data = NULL,

  # --- 置换检验参数 ---
  n_perm = 1000,
  seed = 123,
  epsilon = 1e-9,
  weight_p = 1.0 # 新增：权重指数，标准GSEA通常为1，设为0则退回无权
) {

  if (!is.null(seed)) set.seed(seed)

  # ============================================================================
  # 1. 内部核心函数：计算加权 ES 和 W 统计量
  # ============================================================================
  # 修改：现在需要接收 metric_vec (排序指标值) 来计算权重
  calc_weighted_GSEA_stat <- function(is_in_set_vec, metric_vec) {
    N <- length(is_in_set_vec)
    N_H <- sum(is_in_set_vec)

    if (N_H == 0 | N_H == N) {
      return(list(W = 0, ES = 0, tau = NA, direction = "NEUTRAL", rs_values = rep(0, N)))
    }

    # --- 加权逻辑开始 ---
    # 1. 取指标的绝对值并应用权重指数 p (通常 p=1)
    abs_metric <- abs(metric_vec)
    if (weight_p != 1) abs_metric <- abs_metric ^ weight_p

    # 2. 计算归一化因子 N_R (所有在集合内基因的权重之和)
    # 仅累加 is_in_set_vec 为 TRUE 的基因
    N_R <- sum(abs_metric[is_in_set_vec])

    # 3. 定义步长向量
    # Hit: |r|^p / N_R
    # Miss: -1 / (N - N_H)
    if (N_R == 0) {
      # 极端情况：集合内基因指标全为0
      step_vec <- ifelse(is_in_set_vec, 0, -1/(N - N_H))
    } else {
      step_vec <- ifelse(is_in_set_vec,
                         abs_metric / N_R,       # Hit step
                         -1 / (N - N_H))         # Miss step
    }
    # --- 加权逻辑结束 ---

    # 随机游走路径
    rs_values <- cumsum(step_vec)

    # 寻找波峰波谷
    max_es <- max(rs_values)
    min_es <- min(rs_values)
    max_tau <- which.max(rs_values)
    min_tau <- which.min(rs_values)

    if (max_es > abs(min_es)) {
      ES <- max_es; direction <- "UP"; tau <- max_tau
    } else {
      ES <- min_es; direction <- "DOWN"; tau <- min_tau
    }

    M_effective <- abs(ES)

    if (M_effective <= epsilon) {
      return(list(W = 0, ES = 0, tau = tau, direction = "NEUTRAL", rs_values = rs_values))
    }

    # 计算 W 统计量 (公式保持不变)
    denominator <- tau * (N - tau)
    W_val <- if (denominator == 0) 0 else { (N * M_effective^2) / denominator }

    return(list(
      W = W_val,
      ES = ES,
      tau = tau,
      direction = direction,
      rs_values = rs_values,
      max_es = max_es, min_es = min_es,
      max_tau = max_tau, min_tau = min_tau
    ))
  }

  # ============================================================================
  # 2. 数据准备
  # ============================================================================

  # --- 路径 A: 使用预排序数据 ---
  if (!is.null(pre_sorted_data) && !is.null(gene_set)) {
    sorted_data <- pre_sorted_data
    sorted_data$is_in_set <- sorted_data$gene_id %in% gene_set
    perm_mode <- "gene_set"

  } else {
    # --- 路径 B: 使用原始表达矩阵 ---
    perm_mode <- "phenotype"

    if (!is.null(gene_set_ind)) {
      if (is.null(gene_expression_matrix)) stop("错误: 需提供表达矩阵")
      final_gene_set <- rownames(gene_expression_matrix)[gene_set_ind]
    } else if (!is.null(gene_set)) {
      final_gene_set <- gene_set
    } else {
      stop("错误: 必须提供 'gene_set'")
    }

    # 计算 t-test 排序
    if (is.null(precomputed_ranks)) {
      if (is.null(gene_expression_matrix)) stop("错误: 计算 t-test 需要表达矩阵")
      obs_t_res <- matrixTests::row_t_welch(
        gene_expression_matrix[, 1:num_diff_samples],
        gene_expression_matrix[, (num_diff_samples + 1):ncol(gene_expression_matrix)]
      )
      t_test_results <- obs_t_res$statistic
      names(t_test_results) <- rownames(gene_expression_matrix)
    } else {
      t_test_results <- precomputed_ranks
    }

    # 清理无效值
    t_test_results[is.nan(t_test_results)] <- 0
    if (any(is.infinite(t_test_results))) {
      finite_vals <- t_test_results[is.finite(t_test_results)]
      max_val <- if(length(finite_vals)) max(finite_vals) else 1
      min_val <- if(length(finite_vals)) min(finite_vals) else -1
      t_test_results[t_test_results == Inf] <- max_val + 1
      t_test_results[t_test_results == -Inf] <- min_val - 1
    }

    # 构建排序列表 (必须包含 rank_metric 用于加权)
    gene_ids <- names(t_test_results)
    is_in_set <- gene_ids %in% final_gene_set
    gene_data_df <- data.frame(gene_id = gene_ids, rank_metric = t_test_results, is_in_set = is_in_set)
    # 按指标降序排列
    sorted_data <- gene_data_df[order(gene_data_df$rank_metric, decreasing = TRUE), ]
  }

  # --- 计算观测到的统计量 ---
  # 注意：这里传入了 metric 向量
  obs_res <- calc_weighted_GSEA_stat(sorted_data$is_in_set, sorted_data$rank_metric)
  W_obs <- obs_res$W

  # ============================================================================
  # 3. 置换检验 (Weighted)
  # ============================================================================

  perm_Ws <- numeric(n_perm)

  if (perm_mode == "phenotype" && !is.null(gene_expression_matrix) && is.null(precomputed_ranks)) {
    # 【表型置换】
    total_samples <- ncol(gene_expression_matrix)
    all_indices <- 1:total_samples

    # 获取原始的 logical 向量顺序 (对应矩阵行顺序)
    raw_is_in_set <- rownames(gene_expression_matrix) %in% final_gene_set

    for (i in 1:n_perm) {
      shuffled_indices <- sample(all_indices)

      # 重新计算 t 值
      perm_res <- matrixTests::row_t_welch(
        gene_expression_matrix[, shuffled_indices[1:num_diff_samples]],
        gene_expression_matrix[, shuffled_indices[(num_diff_samples + 1):total_samples]]
      )
      perm_stats <- perm_res$statistic
      perm_stats[is.na(perm_stats)] <- 0

      # 获取新排序
      rank_idx <- order(perm_stats, decreasing = TRUE)

      # 获取排序后的 logical 向量
      perm_sorted_is_in_set <- raw_is_in_set[rank_idx]
      # 获取排序后的 metric 向量 (用于加权)
      perm_sorted_metric <- perm_stats[rank_idx]

      # 计算加权 W
      perm_Ws[i] <- calc_weighted_GSEA_stat(perm_sorted_is_in_set, perm_sorted_metric)$W
    }

  } else {
    # 【基因集置换】
    # 固定 metric 的值，只打乱 is_in_set 标签
    # 这是加权 GSEA 中基因集置换的标准做法

    current_bool_vec <- sorted_data$is_in_set
    fixed_metric_vec <- sorted_data$rank_metric # 指标值保持不变

    for (i in 1:n_perm) {
      perm_bool_vec <- sample(current_bool_vec)
      # 计算加权 W (传入打乱的标签 + 固定的权重)
      perm_Ws[i] <- calc_weighted_GSEA_stat(perm_bool_vec, fixed_metric_vec)$W
    }
  }

  # ============================================================================
  # 4. 结果汇总
  # ============================================================================

  pvalue <- (sum(perm_Ws >= W_obs) + 1) / (n_perm + 1)

  return(list(
    ES = obs_res$ES,
    pvalue = pvalue,
    log_pvalue = log(pvalue),
    W = W_obs,
    perm_Ws = perm_Ws,
    perm_mode = perm_mode,
    direction = obs_res$direction,
    tau = obs_res$tau,
    rs_values = obs_res$rs_values,
    sorted_data = sorted_data
  ))
}


# ==============================================================================
# === 2. 预处理与参数估计 (Preprocessing & Parameter Estimation) ===
# ==============================================================================

#' Estimate ZTNB parameters (lambda, phi) and dropout probabilities p_ij
#'
#' @param counts matrix of raw counts (genes x samples)
#' @param lib_size optional vector of library sizes N_i (length = ncol(counts)); default colSums(counts)
#' @param max_iter max iterations for per-gene lambda/phi iteration (default 50)
#' @param tol convergence tol for lambda/phi (default 1e-6)
#' @param min_phi lower bound for phi to avoid division by zero (default 1e-8)
#' @param subsample_fraction fraction of gene*cell pairs to use to fit the GAM for dropout (0,1], default 1
#' @param verbose logical
#'
#' @return list: Lambda (named vector), Phi (named vector), P_mat (matrix genes x samples), gam_model (fitted mgcv object)
estimate_ztnb_parameters <- function(counts,
                                     lib_size = colSums(counts),
                                     max_iter = 50,
                                     tol = 1e-6,
                                     min_phi = 1e-8,
                                     subsample_fraction = 1,
                                     verbose = TRUE) {
  # Input checks
  if (!is.matrix(counts) && !inherits(counts, "matrix")) {
    stop("counts must be a matrix (genes x samples).")
  }
  G <- nrow(counts)
  S <- ncol(counts)
  if (length(lib_size) != S) stop("lib_size length must equal number of samples (columns).")
  gene_names <- rownames(counts)
  if (is.null(gene_names)) gene_names <- paste0("Gene", seq_len(G))
  sample_names <- colnames(counts)
  if (is.null(sample_names)) sample_names <- paste0("Sample", seq_len(S))

  # Ensure lib_size positive
  lib_size <- as.numeric(lib_size)
  lib_size[lib_size <= 0] <- min(lib_size[lib_size > 0], na.rm = TRUE)

  # Precompute indices of nonzero entries for each gene
  nonzero_idx_list <- lapply(seq_len(G), function(j) which(counts[j, ] > 0))

  # Initialize Lambda and Phi
  # Start: lambda_j = sum(nonzero Y_ij) / sum(N_i over those i)  (per paper idea)
  Lambda <- numeric(G)
  Phi <- numeric(G)
  names(Lambda) <- gene_names
  names(Phi) <- gene_names

  for (j in seq_len(G)) {
    idx <- nonzero_idx_list[[j]]
    if (length(idx) == 0) {
      Lambda[j] <- 0
      Phi[j] <- min_phi
    } else {
      Lambda[j] <- sum(counts[j, idx]) / sum(lib_size[idx])
      # Method of moments init for phi: use sample variance across nonzero normalized by Ni
      mu_i <- Lambda[j] * lib_size[idx]
      obs <- counts[j, idx]
      # If few points, set small phi
      if (length(obs) > 1) {
        m1 <- mean(obs)
        m2 <- mean(obs^2)
        # var = m2 - m1^2 ; for NB var = mu + phi mu^2 => phi ~ (var - mu) / mu^2
        var_obs <- var(obs)
        phi_init <- (var_obs - mean(mu_i)) / (mean(mu_i)^2)
        if (is.na(phi_init) || phi_init <= 0) phi_init <- min_phi
        Phi[j] <- max(phi_init, min_phi)
      } else {
        Phi[j] <- min_phi
      }
    }
  }

  # Iterative updates per gene using formulas (6) and (7)
  if (verbose) cat("Starting iterative ZTNB parameter estimation...\n")
  for (j in seq_len(G)) {
    idx <- nonzero_idx_list[[j]]
    if (length(idx) == 0) next
    lambda_old <- Lambda[j]
    phi_old <- Phi[j]

    # iterate
    for (iter in seq_len(max_iter)) {
      mu_i <- lambda_old * lib_size[idx]  # vector of means for the nonzero cells
      # ensure mu positive
      mu_i[mu_i <= 0] <- 1e-8

      # compute P0 under NB with mean mu_i and dispersion phi_old
      size_param <- 1 / max(phi_old, min_phi)   # size = 1/phi
      # avoid infinite size: cap it
      size_param <- pmin(size_param, 1e12)
      P0 <- dnbinom(0, size = size_param, mu = mu_i)
      # prevent extreme rounding issues
      P0 <- pmin(pmax(P0, 0), 1)

      # update lambda: numerator = sum_i Y_ij * (1 - P0) ; denom = sum_i N_i
      num_lambda <- sum(counts[j, idx] * (1 - P0))
      denom_lambda <- sum(lib_size[idx])
      lambda_new <- ifelse(denom_lambda > 0, num_lambda / denom_lambda, 0)
      lambda_new <- max(lambda_new, 0)

      # update phi:
      # numerator_phi = sum_i Y_ij^2 * (1 - P0) - sum_i (lambda_old * N_i)^2
      # denom_phi = sum_i (lambda_old * N_i)
      num_phi <- sum((counts[j, idx]^2) * (1 - P0)) - sum((lambda_old * lib_size[idx])^2)
      denom_phi <- sum(lambda_old * lib_size[idx])
      if (denom_phi <= 0) {
        phi_new <- min_phi
      } else {
        phi_new <- num_phi / denom_phi
        if (is.na(phi_new) || !is.finite(phi_new) || phi_new <= 0) phi_new <- min_phi
      }

      # check convergence
      if (abs(lambda_new - lambda_old) < tol && abs(phi_new - phi_old) < tol) {
        lambda_old <- lambda_new
        phi_old <- phi_new
        break
      }
      lambda_old <- lambda_new
      phi_old <- phi_new
    } # end iter
    Lambda[j] <- lambda_old
    Phi[j] <- max(phi_old, min_phi)
    if (verbose && (j %% 500 == 0)) cat(sprintf("   ... processed %d / %d genes\n", j, G))
  } # end per gene loop
  if (verbose) cat("ZTNB parameter iteration finished.\n")

  # --- Estimate dropout probabilities p_ij via gam ---
  # z_ij indicator of zero (1 if zero)
  if (verbose) cat("Preparing data for dropout GAM fitting...\n")
  # gene-level A_j: average log CPM (use aveLogCPM approx: log10((counts / lib_size*1e6)+1) mean)
  # compute per-gene average logCPM (A_j)
  # (we use log-scale base e for mgcv; paper uses log-scale, exact base not critical)
  pseudo_cpm <- t(t(counts) / lib_size * 1e6)
  A_j <- rowMeans(log1p(pseudo_cpm))  # average log(1 + CPM)

  # build data frame of (z, A_j, logN) for all gene-sample pairs
  # This can be huge; optionally subsample
  Gs <- rep(seq_len(G), each = S)
  Ss <- rep(seq_len(S), times = G)
  z_vec <- as.vector((counts == 0) * 1)
  A_vec <- A_j[Gs]
  logN_vec <- log(lib_size[Ss])

  df_gam <- data.frame(z = z_vec, A = A_vec, logN = logN_vec)

  # subsample if requested
  if (subsample_fraction < 1 && subsample_fraction > 0) {
    set.seed(123)
    keep_idx <- sample.int(nrow(df_gam), size = ceiling(nrow(df_gam) * subsample_fraction))
    df_fit <- df_gam[keep_idx, , drop = FALSE]
  } else {
    df_fit <- df_gam
  }

  if (verbose) cat("Fitting GAM for dropout p_ij (may be slow for very large datasets)...\n")
  # formula: z ~ s(A) + logN + s(A, by = logN)
  # use binomial family (logit)
  gam_fit <- mgcv::gam(z ~ s(A) + logN + s(A, by = logN), family = binomial(link = "logit"), data = df_fit, method = "REML")
  if (verbose) cat("GAM fitted. Predicting p_ij for all gene-sample pairs...\n")

  # predict for all pairs (use full df_gam)
  p_pred <- predict(gam_fit, newdata = df_gam, type = "response")
  # Bound
  p_pred <- pmin(pmax(p_pred, 1e-8), 1 - 1e-8)
  P_mat <- matrix(p_pred, nrow = G, ncol = S, byrow = FALSE)
  rownames(P_mat) <- gene_names
  colnames(P_mat) <- sample_names

  if (verbose) cat("Done. Returning Lambda, Phi, P_mat, gam_model.\n")
  return(list(
    Lambda = Lambda,
    Phi = Phi,
    P_mat = P_mat,
    gam_model = gam_fit
  ))
}

# ==============================================================================
# === 1. 数据模拟与生成 (Data Simulation & Generation) ===
# ==============================================================================

#' @title 生成模拟基因表达数据 (简化版，按样本数分组)
#'
#' @description
#' 此函数生成模拟的基因表达数据。它首先为每个基因创建基础表达水平，
#' 然后随机选择一个基因集，并将样本分为两组。
#' 最后，对第一组样本中属于基因集的基因进行表达上调。
#'
#' @param num_genes 整数，生成的总基因数。
#' @param num_samples 整数，生成的总样本数。
#' @param CR 数字，作为基因集的基因所占比例（0到1之间）。
#' @param num_diff_samples 整数，指定Group2（非上调组）的样本数量。
#' @param base_mean 数字，基础表达均值的均值。
#' @param base_sd 数字，基础表达均值的标准差。
#' @param expression_sd 数字，每个基因内部表达量的标准差。
#' @param upreg_mean 数字，上调表达量的均值。
#' @param upreg_sd 数字，上调表达量的标准差。
#'
#' @return 一个列表，包含表达矩阵、基因集和样本分组信息。
generate_data_norm_simple <- function(num_genes = 1000,
                                      num_samples = 20,
                                      CR = 0.1,
                                      num_diff_samples = 10,
                                      base_mean = 8,
                                      base_sd = 2,
                                      expression_sd = 1,
                                      upreg_mean = 2,
                                      upreg_sd = 0.5) {


  gene_base_means <- rnorm(num_genes, mean = base_mean, sd = base_sd)
  expression_matrix <- matrix(NA, nrow = num_genes, ncol = num_samples)
  for (i in 1:num_genes) {
    expression_matrix[i, ] <- rnorm(num_samples, mean = gene_base_means[i], sd = expression_sd)
  }
  rownames(expression_matrix) <- paste0("Gene", 1:num_genes)
  colnames(expression_matrix) <- paste0("Sample", 1:num_samples)

  n_gene_set <- floor(num_genes * CR)
  gene_set <- sample(rownames(expression_matrix), size = n_gene_set, replace = FALSE)

  n_group1 <- num_samples - num_diff_samples
  if (n_group1 < 0) {
    stop("Error: num_diff_samples (Group2 size) cannot be larger than num_samples (total size).")
  }
  group1_samples <- sample(colnames(expression_matrix), size = n_group1, replace = FALSE)
  sample_groups <- ifelse(colnames(expression_matrix) %in% group1_samples, "Group1", "Group2")
  names(sample_groups) <- colnames(expression_matrix)

  # --- 4. 对特定基因和样本进行表达上调 ---
  # 找到需要被上调的基因（行）和样本（列）的索引
  gene_set_indices <- which(rownames(expression_matrix) %in% gene_set)
  group1_indices <- which(colnames(expression_matrix) %in% group1_samples)

  if (length(gene_set_indices) > 0 && length(group1_indices) > 0) {
    num_upregulations <- length(gene_set_indices) * length(group1_indices)
    upregulation_values <- rnorm(num_upregulations, mean = upreg_mean, sd = upreg_sd)
    expression_matrix[gene_set_indices, group1_indices] <-
      expression_matrix[gene_set_indices, group1_indices] + upregulation_values
  }

  # --- 5. 返回结果 ---
  return(list(
    expression_data = expression_matrix,
    gene_set_indices = gene_set_indices,
    sample_groups = sample_groups
  ))
}

#' @title Generate Normally Distributed Gene Expression Data
#' @description Creates a matrix of gene expression data following a Normal distribution,
#'              simulating control and treatment groups with differentially expressed (DE) genes.
#'
#' @param num_genes Number of genes (rows) in the matrix.
#' @param num_samples Number of total samples (columns).
#' @param tau0 Baseline log-odds for a gene to be DE.
#' @param tau1 Additional log-odds for a gene in the core gene set to be DE.
#' @param num_diff_samples Number of treatment samples.
#' @param CR "Common Response" rate. The proportion of genes belonging to a core set with a higher chance of being DE.
#' @param sd A vector of standard deviations for each gene across all samples. Length must be `num_genes`.
#' @param Lambda A vector of baseline mean expression levels for each gene.
#'
#' @return A list containing:
#'         1. `gene_expression_matrix`: The generated data matrix.
#'         2. `gene_set_ind`: The indices of genes in the core set.
#'         3. `DE_gene_ind`: The indices of genes that are truly differentially expressed.
generate_data_norm <- function(num_genes, # 矩阵的行
                               num_samples, # 矩阵的列
                               tau0,
                               tau1,
                               num_diff_samples, # 处理组样本数
                               CR , # 决定基因集大小
                               sd, # 向量，每个基因的标准差
                               Lambda) { # 向量，每个基因的平均表达水平

  # 检查sd参数的维度
  if (length(sd) != num_genes) {
    stop("Length of 'sd' vector must be equal to 'num_genes'.")
  }

  # 步骤 1: 选择DE基因
  # ---------------------------------------------------
  gene_set_ind <- sample(1:num_genes, num_genes * CR)
  a <- numeric(num_genes)#指示变量。表示当前基因是否在基因集中
  a[gene_set_ind] <- 1
  p <- exp(tau0 + a * tau1) / (1 + exp(tau0 + a * tau1)) # 成为DE基因的概率
  DE_gene_ind <- which(runif(num_genes) < p)

  # 步骤 2: 计算差异表达的效应量大小
  # ---------------------------------------------------
  beta <- numeric(num_genes)
  num_DE_gene <- length(DE_gene_ind)
  beta0 <- rnorm(num_DE_gene, mean = 0, sd = 3.5) # 效应量的大小可以自行调整
  beta[DE_gene_ind] <- beta0
  fc <- exp(beta) # Fold Change

  # 步骤 3: 生成最终的符合正态分布的表达矩阵
  # ---------------------------------------------------
  num_control_samples <- num_samples - num_diff_samples

  # 构建均值(mean)和标准差(sd)的向量，用于一次性生成所有数据
  # 均值向量
  mu_control_part <- rep(Lambda, times = num_control_samples)
  mu_treatment_part <- rep(Lambda * fc, times = num_diff_samples)
  mu_vector <- c(mu_control_part, mu_treatment_part)

  # 标准差向量
  sd_vector <- rep(sd, times = num_samples)

  # 使用 rnorm 生成所有数据
  all_values <- rnorm(n = num_genes * num_samples,
                      mean = mu_vector,
                      sd = sd_vector)

  gene_expression_matrix <- matrix(all_values, nrow = num_genes, ncol = num_samples)

  # 步骤 4: 添加行名和列名
  # ---------------------------------------------------
  rownames(gene_expression_matrix) <- paste0("Gene", 1:num_genes)
  colnames(gene_expression_matrix) <- paste0("Sample", 1:num_samples)

  # 返回结果
  return(list(
    gene_expression_matrix = gene_expression_matrix,
    gene_set_ind = gene_set_ind,
    DE_gene_ind = DE_gene_ind
  ))
}
