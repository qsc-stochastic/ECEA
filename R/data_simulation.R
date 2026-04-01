generate_data_nb <- function(num_genes, num_samples, tau0, tau1,
                             num_diff_samples, CR, dispersions, Lambda,
                             mean_logFC = 0, sd_logFC = 0.8) {

  # 1. 选择 DE 基因
  gene_set_ind <- sample(1:num_genes, num_genes * CR)
  a <- numeric(num_genes)
  a[gene_set_ind] <- 1
  p <- exp(tau0 + a * tau1) / (1 + exp(tau0 + a * tau1))
  DE_gene_ind <- which(runif(num_genes) < p)

  # 2. 效应量大小变化 (使用更真实的参数，正负对称)
  beta <- numeric(num_genes)
  num_DE_gene <- length(DE_gene_ind)
  # 只有 DE 基因才赋予非零的 LogFC
  if(num_DE_gene > 0){
    beta[DE_gene_ind] <- rnorm(num_DE_gene, mean = mean_logFC, sd = sd_logFC)
  }
  fc <- exp(beta)

  # 3. 组装表达均值矩阵 (按列优先填入矩阵)
  num_control_samples <- num_samples - num_diff_samples
  mu_control_part <- rep(Lambda, times = num_control_samples)
  mu_treatment_part <- rep(Lambda * fc, times = num_diff_samples)

  # 直接将均值向量拼成矩阵，方便后续分别计算对照组和处理组的 Dropout
  mu_matrix <- matrix(c(mu_control_part, mu_treatment_part),
                      nrow = num_genes, ncol = num_samples, byrow = FALSE)

  # 4. 动态生成符合当前表达水平的 Dropout 概率矩阵 P
  # (表达量升高，dropout 概率自动下降)
  P_dynamic <- 1 / (1 + exp(2 * (log(mu_matrix + 1e-6) - 1)))

  # 5. 生成负二项分布 counts
  size_matrix <- matrix(rep(1/dispersions, times = num_samples),
                        nrow = num_genes, ncol = num_samples, byrow = FALSE)

  all_counts <- rnbinom(n = num_genes * num_samples,
                        mu = as.vector(mu_matrix),
                        size = as.vector(size_matrix))
  gene_expression_matrix <- matrix(all_counts, nrow = num_genes, ncol = num_samples)

  # 6. 应用 Dropout 掩码
  random_uniform_matrix <- matrix(runif(num_genes * num_samples),
                                  nrow = num_genes, ncol = num_samples)
  dropout_indices <- (random_uniform_matrix < P_dynamic)
  gene_expression_matrix[dropout_indices] <- 0

  rownames(gene_expression_matrix) <- paste0("Gene", 1:num_genes)
  colnames(gene_expression_matrix) <- paste0("Sample", 1:num_samples)

  return(list(gene_expression_matrix, gene_set_ind))
}
