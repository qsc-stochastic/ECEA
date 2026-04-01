read_gct <- function(file) {
  lines <- readLines(file)
  data_df <- read.table(text = lines[-c(1, 2)], header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, quote = "", comment.char = "", na.strings = "na")
  if ("DESCRIPTION" %in% colnames(data_df)) { data_df <- data_df %>% dplyr::select(-DESCRIPTION) }
  data_df[] <- lapply(data_df, function(x) as.numeric(as.character(x)))
  return(as.matrix(data_df))
}




read_cls <- function(file) {
  lines <- readLines(file, n = 3)
  class_levels <- unlist(strsplit(sub("# ?", "", lines[2]), " "))
  phenotypes <- unlist(strsplit(lines[3], " "))
  return(list(levels = class_levels, phenotypes = phenotypes))
}



read_gmt <- function(gmt_file) {
  lines <- readLines(gmt_file)
  gene_sets <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t")[[1]]
    list(name = fields[1],
         description = fields[2],
         genes = fields[-(1:2)])
  })
  names(gene_sets) <- sapply(gene_sets, function(x) x$name)
  return(gene_sets)
}


# --- 计算 GSEA 类的 M, tau 和 W 统计量 ---
calculate_lambda_gc <- function(p_values) {
  p_values_clean <- na.omit(p_values)
  if (length(p_values_clean) == 0) return(NA)
  chi_sq_stats <- qchisq(p_values_clean, df = 1, lower.tail = FALSE)
  lambda <- median(chi_sq_stats) / qchisq(0.5, df = 1)
  return(lambda)
}

#布朗桥最大值的密度函数
f_M <- function(m) {
  result <- ifelse(m > 0, 4 * m * exp(-2 * m^2), NA)

  if (any(m <= 0)) {
    warning("输入值 'm' 必须大于 0。对于非正数将返回 NA。")
  }

  return(result)
}

# 布朗桥最大值的累积分布函数 (CDF)
F_M <- function(m) {
  result <- ifelse(m > 0, 1 - exp(-2 * m^2), 0)

  return(result)
}
