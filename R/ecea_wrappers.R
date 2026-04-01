#' @title 执行GSEA富集分析 (使用clusterProfiler和信噪比S2N排序)
#' @description
#' 该函数封装了完整的GSEA分析流程。它首先使用信噪比(Signal-to-Noise Ratio, S2N)
#' 对基因进行排序，然后利用clusterProfiler包中的GSEA函数执行富集分析。
#' 这是复现GSEA原始论文分析策略的经典方法。
#'
#' @param expression_matrix 一个数值矩阵，行为基因(需要有行名, rownames)，列为样本。
#'        这是标准的基因表达矩阵。
#' @param phenotype_vector 一个字符向量或因子，其长度与表达矩阵的列数相同。
#'        它定义了每个样本所属的组别 (例如, c("Tumor", "Tumor", "Normal"))。
#' @param gene_sets_df 一个数据框 (data.frame或tibble)，至少包含两列：
#'        一列是基因集名称 (gs_name)，另一列是对应的基因符号 (gene_symbol)。
#'        这是进行富集分析的基因集定义文件。
#' @param group1_name 一个字符串，指定表型向量中的第一个组别名称 (例如 "p53_WT")。
#'        这个组将作为S2N计算公式中的分子部分的"正项"，其高表达基因会排在列表顶端。
#' @param group2_name 一个字符串，指定表型向量中的第二个组别名称 (例如 "p53_MUT")。
#'        这个组将作为S2N计算公式中的分子部分的"负项"。
#' @param n_permutations 一个整数，指定GSEA进行置换检验的次数。
#'        默认为1000。增加此值会提高结果的稳定性，但会增加计算时间。
#'
#' @return 一个tibble (一种增强型数据框)。如果GSEA分析成功，返回包含富集结果的表格，
#'         列包括ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue等。
#'         如果分析没有返回任何结果或出错，则返回一个空的tibble。
#'
#' @importFrom dplyr select
#' @importFrom tibble as_tibble
#' @importFrom clusterProfiler GSEA
#'
#' @examples
#' # # 这是一个示例，实际运行时需要真实的表达矩阵、表型和基因集数据
#' # expression_matrix <- matrix(rnorm(1000 * 10), nrow = 1000,
#' #                             dimnames = list(paste0("Gene", 1:1000), paste0("Sample", 1:10)))
#' # phenotype_vector <- c(rep("GroupA", 5), rep("GroupB", 5))
#' # gene_sets <- data.frame(
#' #   gs_name = rep(c("Pathway1", "Pathway2"), each = 50),
#' #   gene_symbol = paste0("Gene", sample(1:1000, 100))
#' # )
#' #
#' # gsea_results <- run_gsea_with_clusterProfiler_s2n(
#' #   expression_matrix = expression_matrix,
#' #   phenotype_vector = phenotype_vector,
#' #   gene_sets_df = gene_sets,
#' #   group1_name = "GroupA",
#' #   group2_name = "GroupB"
#' # )
#' # print(gsea_results)

run_gsea_with_clusterProfiler_s2n <- function(expression_matrix, phenotype_vector, gene_sets_df,
                                              group1_name, group2_name, n_permutations = 1000, minGSSize = 15, maxGSSize = 500) {

  # --- 步骤 1: 计算基因排序列表 (Gene Rank List) ---
  # GSEA的第一步是根据基因在两个表型组之间的差异表达程度，对所有基因进行排序。
  # 这里我们使用经典的"信噪比"(Signal-to-Noise Ratio, S2N)作为排序指标。
  # S2N = (group1的平均表达 - group2的平均表达) / (group1的标准差 + group2的标准差)
  cat("    -> Calculating gene ranks using Signal-to-Noise Ratio (S2N)...\n")

  # 找到每个组别对应的样本索引 (即在矩阵中的列号)
  group1_indices <- which(phenotype_vector == group1_name)
  group2_indices <- which(phenotype_vector == group2_name)
  cat(sprintf("    DEBUG: Found %d samples for group '%s' and %d samples for group '%s'.\n", length(group1_indices), group1_name, length(group2_indices), group2_name))

  # 健壮性检查：确保两个组都有样本，防止后续计算出错
  if (length(group1_indices) == 0 || length(group2_indices) == 0) {
    stop("Error: One or both phenotype groups have zero samples.")
  }

  # 分别计算每个基因 (每一行) 在两组中的平均表达量
  mean_g1 <- rowMeans(expression_matrix[, group1_indices, drop=F], na.rm = TRUE)
  mean_g2 <- rowMeans(expression_matrix[, group2_indices, drop=F], na.rm = TRUE)

  # 分别计算每个基因 (每一行) 在两组中的标准差
  # 特殊处理：如果组内只有一个样本，标准差为0
  sd_g1 <- if(length(group1_indices) > 1) apply(expression_matrix[, group1_indices, drop=F], 1, sd, na.rm = TRUE) else 0
  sd_g2 <- if(length(group2_indices) > 1) apply(expression_matrix[, group2_indices, drop=F], 1, sd, na.rm = TRUE) else 0

  # 为了防止分母为0 (当两个组的标准差都为0时)，加入一个极小值epsilon
  epsilon <- 1e-6
  gene_ranks <- (mean_g1 - mean_g2) / (sd_g1 + sd_g2 + epsilon)

  # 将基因名赋给排序分数值，这是GSEA必须的格式
  names(gene_ranks) <- rownames(expression_matrix)

  # 处理计算过程中可能出现的NA或NaN值，将其替换为0
  gene_ranks[is.na(gene_ranks) | is.nan(gene_ranks)] <- 0
  cat("    DEBUG: Summary of calculated gene ranks (S2N):\n")
  print(summary(gene_ranks))

  # 按S2N分值从高到低对基因进行排序，这是GSEA算法的直接输入
  gene_ranks <- sort(gene_ranks, decreasing = TRUE)

  # --- 步骤 2: 准备GSEA的输入格式 ---
  # clusterProfiler::GSEA函数需要一个特定格式的TERM2GENE数据框，
  # 其中一列是基因集ID(term)，另一列是基因ID(gene)。
  cat("    -> Preparing TERM2GENE and renaming columns to 'term' and 'gene'...\n")
  term2gene <- gene_sets_df %>% dplyr::select(gs_name, gene_symbol)
  colnames(term2gene) <- c("term", "gene")

  # --- 步骤 3: 运行GSEA分析 ---
  cat("    -> Running GSEA using clusterProfiler::GSEA...\n")
  # 设置随机数种子以保证结果的可重复性
  set.seed(42)

  # 调用clusterProfiler的核心GSEA函数
  gsea_res_obj <- GSEA(
    geneList = gene_ranks,      # 输入排序好的基因列表
    TERM2GENE = term2gene,        # 输入基因集定义文件
    minGSSize = minGSSize,             # 过滤掉成员基因数小于minGSSize的基因集 (这是在取交集后的数量)
    maxGSSize = maxGSSize,            # 过滤掉成员基因数大于maxGSSize的基因集 (这也是取交集后的数量)
    pvalueCutoff = 1,           # 设置为1，意味着我们获取所有基因集的计算结果，之后再根据需要进行筛选
    nPermSimple = n_permutations, # 置换检验的次数
    verbose = FALSE,            # 不在控制台打印详细的运行过程
    pAdjustMethod = "BH"          # 使用Benjamini-Hochberg方法对p值进行多重检验校正
  )
  # ================================================================= #

  # --- 步骤 4: 结果处理和返回 ---
  # 检查GSEA是否返回了有效结果
  if (is.null(gsea_res_obj) || nrow(gsea_res_obj@result) == 0) {
    cat("    WARNING: GSEA analysis returned no results.\n")
    return(tibble()) # 如果没有结果，返回一个空的tibble，方便后续处理
  }

  # GSEA函数返回的是一个S4对象，我们将其中的结果表格(@result)提取出来，
  # 并转换为更易于操作的tibble格式。
  return(as_tibble(gsea_res_obj@result))
}









run_ECEA_analysis <- function(gene_expression_matrix,
                              gene_sets_list,       # 接收你的基因集列表 (如 c2_sets)
                              num_diff_samples,
                              VI,                   # 开关：是否启用方差校正因子 (TRUE/FALSE)
                              min_gs_size = 15,     # 最小基因集大小阈值
                              ranking_metric = "t"  # 排序指标: "t", "S2N" 或 "DESeq2"
) {

  # ==============================================================================
  #  计算排序指标 (Metric Calculation)
  # ==============================================================================
  cat(sprintf("\n[ECEA 初始化] 正在进行预计算 (Matrix Size: %d x %d)...\n",
              nrow(gene_expression_matrix), ncol(gene_expression_matrix)))
  cat(sprintf("[ECEA 设置] 排序指标: %s\n", ranking_metric))

  mat1 <- gene_expression_matrix[, 1:num_diff_samples, drop=FALSE]
  mat2 <- gene_expression_matrix[, (num_diff_samples + 1):ncol(gene_expression_matrix), drop=FALSE]

  metric_values <- numeric(nrow(gene_expression_matrix))

  if (ranking_metric == "t") {
    t_test_results_res <- matrixTests::row_t_welch(mat1, mat2)
    metric_values <- t_test_results_res$statistic

  } else if (ranking_metric == "S2N") {
    mean1 <- rowMeans(mat1, na.rm = TRUE)
    mean2 <- rowMeans(mat2, na.rm = TRUE)
    n1 <- ncol(mat1); n2 <- ncol(mat2)
    if (n1 < 2 || n2 < 2) stop("错误：使用 S2N 排序时，每组样本数必须至少为 2。")
    sd1 <- sqrt(rowSums((mat1 - mean1)^2, na.rm = TRUE) / (n1 - 1))
    sd2 <- sqrt(rowSums((mat2 - mean2)^2, na.rm = TRUE) / (n2 - 1))
    metric_values <- (mean1 - mean2) / (sd1 + sd2 + 1e-9)

  } else if (ranking_metric == "DESeq2") {
    message("⚠️ 提示: 使用 DESeq2 排序时，请绝对确保输入的 gene_expression_matrix 为未标准化的 Raw Counts！")
    condition <- factor(rep(c("G1", "G2"), times = c(num_diff_samples, ncol(gene_expression_matrix) - num_diff_samples)))
    coldata <- data.frame(condition = condition, row.names = colnames(gene_expression_matrix))
    dds <- DESeqDataSetFromMatrix(countData = round(gene_expression_matrix),
                                  colData = coldata,
                                  design = ~ condition)
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds)
    metric_values <- res$stat

  } else {
    stop("错误：ranking_metric 参数必须是 't', 'S2N' 或 'DESeq2'")
  }

  names(metric_values) <- rownames(gene_expression_matrix)

  # ==============================================================================
  #  清理无效值 (NaN, Inf) 并全局预排序
  # ==============================================================================
  if (any(is.nan(metric_values))) { metric_values[is.nan(metric_values)] <- 0 }

  if (any(is.infinite(metric_values))) {
    finite_vals <- metric_values[is.finite(metric_values)]
    max_finite_val <- if(length(finite_vals) > 0) max(finite_vals, na.rm = TRUE) else 1
    min_finite_val <- if(length(finite_vals) > 0) min(finite_vals, na.rm = TRUE) else -1
    metric_values[metric_values == Inf] <- max_finite_val + 1
    metric_values[metric_values == -Inf] <- min_finite_val - 1
  }

  gene_data_df <- data.frame(gene_id = names(metric_values), rank_metric = metric_values)

  pre_sorted_data_for_all <- gene_data_df[order(gene_data_df$rank_metric, decreasing = TRUE), ]

  # ==============================================================================
  #  基因集过滤 (Filtering)
  # ==============================================================================
  gs_sizes <- sapply(gene_sets_list, length)
  valid_indices <- gs_sizes >= min_gs_size

  n_original <- length(gene_sets_list)
  gene_sets_list <- gene_sets_list[valid_indices]
  n_final <- length(gene_sets_list)

  cat(sprintf("[ECEA 过滤] 原始通路数: %d | 移除过小通路: %d | 剩余有效: %d\n",
              n_original, n_original - n_final, n_final))

  if (n_final == 0) return(NULL)


  cat("[ECEA 计算] 开始逐个分析通路...\n")

  all_results <- lapply(names(gene_sets_list), function(pathway_name) {
    gene_indices <- gene_sets_list[[pathway_name]]
    # 如果传来的是索引，转成基因名；如果是名字，直接用
    if (is.numeric(gene_indices)) {
      pathway_gene_names <- rownames(gene_expression_matrix)[gene_indices]
    } else {
      pathway_gene_names <- gene_indices
    }
    pathway_size <- length(pathway_gene_names)

    # 根据 VI 开关决定调用哪个核心算法，并提取对应的值
    if(!VI){
      single_result <- ECEA(
        pre_sorted_data = pre_sorted_data_for_all,
        gene_set = pathway_gene_names
      )
      vif_val <- 1.0     # 没开 VIF，就是 1.0
      rho_val <- NA      # 没算相关性，记为 NA
    } else {
      single_result <- ECEA_VIF(
        pre_sorted_data = pre_sorted_data_for_all,
        gene_set = pathway_gene_names,
        gene_expression_matrix = gene_expression_matrix,
        use_correction = TRUE
      )
      vif_val <- single_result$VIF
      rho_val <- single_result$mean_rho
    }

    return(data.frame(
      pathway = pathway_name,
      setSize = pathway_size,
      ES = single_result$ES,
      pvalue = single_result$pvalue,
      VIF = vif_val,
      mean_rho = rho_val,
      direction = single_result$direction,
      stringsAsFactors = FALSE
    ))
  })

  # 整理结果并计算 FDR
  results_df <- do.call(rbind, all_results)
  results_df$FDR <- p.adjust(results_df$pvalue, method = "BH")

  results_df <- results_df[order(results_df$FDR),
                           c("pathway", "setSize", "ES", "pvalue", "FDR", "mean_rho", "VIF", "direction")]
  rownames(results_df) <- NULL

  return(results_df)
}
