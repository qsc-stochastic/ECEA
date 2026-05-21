#' 为GSEA结果绘制条形图 (终极版：带分类开关)
#'
#' 这个函数接收GSEA分析的结果，并生成一个条形图来可视化富集的通路。
#' 增加了一个 group_by_category 参数，可以灵活控制是否按通路来源进行分类。
#'
#' @param gsea_results 一个数据框，通常是 'run_gsea_analysis' 函数的输出结果。
#' @param gene_sets_list 用于GSEA分析的原始基因集列表。
#' @param top_n 整数，要显示的最显著通路的数量。
#' @param fdr_cutoff 数值，用于筛选显著结果的FDR阈值。
#' @param group_by_direction 逻辑值，是否按UP/DOWN分面（通常保持TRUE）。
#' @param group_by_category 逻辑值。如果为 TRUE (默认)，将按通路来源(KEGG, REACTOME等)进行分面。
#'                          如果为 FALSE，则只按 UP/DOWN 进行分面。
#'
#' @return 返回一个 ggplot 对象。
#'
plot_gsea_barchart <- function(gsea_results,
                               gene_sets_list,
                               top_n = 10,
                               fdr_cutoff = 0.25,
                               group_by_direction = TRUE, # 保留了这个参数
                               group_by_category = TRUE,# <<< 新增的逻辑开关参数，默认为TRUE
                               plot_title = "GSEA Enriched Pathways") {

  # --- 1. 数据准备 ---

  # --- 1. 数据准备 (修改版) ---

  plot_data <- gsea_results %>%
    filter(FDR < fdr_cutoff) %>%
    dplyr::rename(count = setSize) %>%  # <<< 关键修改：直接用结果里的真实数量 setSize
    filter(!is.na(count))

  if (nrow(plot_data) == 0) {
    message("No significant pathways to plot with the given FDR cutoff.")
    return(invisible(NULL))
  }

  # <<< 修改点 1: 只有当开关为TRUE时，才创建category列 >>>
  if (group_by_category) {
    plot_data <- plot_data %>%
      mutate(
        category = str_extract(pathway, "^[A-Z0-9]+"),
        category = ifelse(category == "HALLMARK", "Hallmark", category)
      )
  }

  # --- 2. 筛选Top N个通路 ---

  # <<< 修改点 2: 根据开关决定分组变量 >>>
  grouping_vars <- character()
  if (group_by_category) {
    grouping_vars <- c(grouping_vars, "category")
  }
  if (group_by_direction) {
    grouping_vars <- c(grouping_vars, "direction")
  }

  # 如果有分组变量，才进行分组筛选
  if (length(grouping_vars) > 0) {
    plot_data <- plot_data %>%
      group_by(!!!rlang::syms(grouping_vars)) %>%
      slice_min(order_by = FDR, n = top_n, with_ties = FALSE) %>%
      ungroup()
  } else {
    # 如果没有任何分组，则在全局筛选top_n
    plot_data <- plot_data %>%
      slice_min(order_by = FDR, n = top_n, with_ties = FALSE)
  }


  # --- 3. 排序与截断 ---

  plot_data <- plot_data %>%
    mutate(pathway = str_trunc(pathway, width = 70, side = "right"))%>%
    # <<< 新增: 生成格式化的 FDR 标签 >>>
    # 如果 FDR 极小 (<0.001)，使用科学计数法，否则保留3位小数
    mutate(fdr_label = ifelse(FDR < 0.001,
                              sprintf("%.1e", FDR),
                              sprintf("%.3f", FDR)))

  # <<< 修改点 3: 根据开关决定排序方式 >>>
  if (group_by_category) {
    plot_data <- plot_data %>% arrange(category, direction, count)
  } else {
    plot_data <- plot_data %>% arrange(direction, count)
  }

  plot_data <- plot_data %>%
    mutate(pathway = factor(pathway, levels = unique(.$pathway)))

  # --- 4. 绘图 (基本不变) ---

  p <- ggplot(plot_data, aes(x = count, y = pathway, fill = FDR)) +
    geom_bar(stat = "identity") +

    # <<< 新增: 添加 FDR 数值标签 >>>
    geom_text(aes(label = fdr_label),
              hjust = -0.2,       # -0.2 表示文字在条形图右侧稍微留点空隙
              size = 3,           # 字体大小，可根据需要调整
              color = "black") +  # 字体颜色

    scale_fill_gradient(low = "#e63946", high = "#457b9d") +

    # <<< 新增: 扩展 X 轴范围 >>>
    # mult = c(0, 0.15) 表示左边不扩展，右边扩展 15% 的空间，防止文字被切断
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      x = "Gene Count in Pathway",
      y = NULL,
      title = plot_title,
      subtitle = paste("Top pathways with FDR <", fdr_cutoff),
      fill = "FDR"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y = element_text(hjust = 1),
      strip.text = element_text(face = "bold", size = 10, color = "white"),
      strip.background = element_rect(fill = "#2a3d45", color = "black"),
      panel.border = element_rect(color = "black", fill = NA),
      legend.position = "right"
    )

  # <<< 修改点 4: 根据开关决定分面规则 >>>
  if (group_by_category) {
    p <- p + facet_grid(
      category + direction ~ .,
      scales = "free_y",
      space = "free_y"
    )
  } else if (group_by_direction) {
    p <- p + facet_grid(
      direction ~ .,
      scales = "free_y",
      space = "free_y"
    )
  }

  return(p)
}




# 确保 ggplot2 和 dplyr 包已加载
# library(ggplot2)
# library(dplyr)

#' 为GSEA结果绘制气泡图
#'
#' 这个函数接收GSEA分析的结果，并生成一个气泡图来可视化富集的通路。
#'
#' @param gsea_results 一个数据框，通常是 'run_gsea_analysis' 函数的输出结果。
#'   必须包含 'pathway', 'ES', 'FDR', 'direction' 这几列。
#' @param gene_sets_list 用于GSEA分析的原始基因集列表，用于获取每个通路的基因数量。
#' @param top_n 整数，在每个方向（UP/DOWN）上，要显示的最显著通路的数量。默认为10。
#' @param fdr_cutoff 数值，用于筛选显著结果的FDR阈值。默认为0.25。
#'
#' @return 返回一个 ggplot 对象，可以直接打印显示。如果无显著结果，则返回NULL。
#'
plot_gsea_bubblechart <- function(gsea_results,
                                  gene_sets_list,
                                  top_n = 10,
                                  fdr_cutoff = 0.25) {

  # --- 1. 数据准备 ---

  # 计算每个通路的基因数量
  gene_counts_df <- tibble(
    pathway = names(gene_sets_list),
    count = sapply(gene_sets_list, length)
  )

  # 筛选FDR显著的结果，并合并基因数量信息
  plot_data <- gsea_results %>%
    filter(FDR < fdr_cutoff) %>%
    left_join(gene_counts_df, by = "pathway") %>%
    filter(!is.na(count)) # 确保有计数的通路才被绘制

  # 如果没有显著结果，则打印消息并退出
  if (nrow(plot_data) == 0) {
    message("No significant pathways to plot with the given FDR cutoff.")
    return(invisible(NULL))
  }

  # --- 2. 筛选Top N个通路 ---

  # 分别对 UP 和 DOWN 富集的通路进行筛选
  plot_data <- plot_data %>%
    group_by(direction) %>%
    slice_min(order_by = FDR, n = top_n, with_ties = FALSE) %>%
    ungroup()

  # --- 3. 排序以优化绘图顺序 ---
  # 为了让图中的Y轴标签按富集分数(ES)排序，我们需要对pathway这个因子进行重排
  # 这会让UP和DOWN两边的图呈现出优雅的倾斜效果
  plot_data <- plot_data %>%
    arrange(direction, ES) %>%
    mutate(pathway = factor(pathway, levels = unique(.$pathway)))

  # --- 4. 绘图 ---

  p <- ggplot(plot_data, aes(x = ES, y = pathway)) +
    # 核心：使用 geom_point 来创建气泡，并映射大小和颜色
    geom_point(aes(size = count, color = FDR)) +

    # 关键：创建两个独立的分面，分别用于UP和DOWN富集的通路
    # scales="free"让每个分面的X轴和Y轴都独立，这对于展示ES值很重要
    facet_wrap(~ direction, scales = "free") +

    # 定义颜色渐变：FDR值越小（越显著），颜色越趋向于红色
    scale_color_gradient(low = "#e63946", high = "#457b9d", name = "FDR") +

    # 定义气泡大小的图例标签
    scale_size_continuous(name = "Gene Count") +

    # 添加坐标轴标签和标题
    labs(
      x = "Enrichment Score (ES)",
      y = NULL,
      title = "GSEA Enriched Pathways (Bubble Plot)"
    ) +

    # 使用一个带网格线的、适合发表的主题
    theme_bw(base_size = 11) +
    theme(
      # 让Y轴的文字（通路名）左对齐，方便阅读
      axis.text.y = element_text(hjust = 1),
      # 自定义分面标签的外观
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "lightgray", color = "black"),
      # 调整图例位置
      legend.position = "right"
    )

  return(p)
}


plot_ECEA <- function(result, gene_set_name = "My Gene Set") {

  # --- 依赖检查 (可选) ---
  # library(ggplot2)
  # library(patchwork)

  # --- 步骤 1: 修正变量名，统一使用 'result' ---
  rs_values    <- result$rs_values
  sorted_data  <- result$sorted_data
  p_value      <- result$pvalue
  es_score     <- result$ES
  direction    <- result$direction
  tau_pos      <- result$tau
  W            <- result$W

  # --- 步骤 2: 创建用于绘图的核心数据框 ---
  plot_data <- data.frame(
    rank = 1:nrow(sorted_data),
    running_sum = rs_values,
    rank_metric = sorted_data$rank_metric
  )

  # 找出基因集“命中”的位置
  hit_indices <- which(sorted_data$is_in_set)

  # --- 面板 1: 富集分数图 ---
  p1 <- ggplot(plot_data, aes(x = rank, y = running_sum)) +
    geom_line(color = "#009E73", linewidth = 1.2) +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") +
    geom_vline(xintercept = tau_pos, color = "#D55E00", linetype = "dotted", linewidth = 1) +
    labs(
      title = paste("Enrichment Plot:", gene_set_name),
      subtitle = sprintf("W = %.2f (%s), p-value = %.3g", W, direction, p_value),
      y = "Enrichment Score"
    ) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(face = "bold", size = 14)
    )

  # --- 面板 2: 基因集命中图 ---
  p2 <- ggplot() +
    geom_segment(aes(x = hit_indices, y = 0, xend = hit_indices, yend = 1), color = "black") +
    scale_x_continuous(limits = c(1, nrow(sorted_data)), expand = c(0, 0)) +
    theme_void()

  # --- 面板 3: 排序度量图 ---
  p3 <- ggplot(plot_data, aes(x = rank, y = rank_metric)) +
    geom_area(aes(fill = rank_metric > 0), show.legend = FALSE) +
    scale_fill_manual(values = c("TRUE" = "#E69F00", "FALSE" = "#56B4E9")) +
    labs(
      x = "Rank in Ordered Dataset",
      y = "Ranking Metric"
    ) +
    theme_classic()

  # --- 最终组合 ---
  combined_plot <- p1 / p2 / p3 +
    plot_layout(heights = c(4, 0.8, 2.5))

  return(combined_plot)
}




# ---- 函数1：绘制功效曲线 ----
create_power_curves <- function(data, methods_prefix, name_map, title_info, log_param, ranking_param) {

  # --- 新增：根据新参数构建副标题 ---
  log_text <- ifelse(log_param, "Log Transform: TRUE", "Log Transform: FALSE")
  ranking_text <- paste("Ranking Method:", ranking_param)
  final_subtitle <- paste(
    title_info,
    paste(log_text, ranking_text, sep = " | "),
    sep = "\n"
  )

  long_data <- data %>%
    pivot_longer(
      cols = all_of(methods_prefix),
      names_to = "method",
      values_to = "power"
    ) %>%
    mutate(method = recode(method, !!!name_map))

  # 为 7 条曲线定义高对比度颜色
  custom_colors <- c(
    "ECEA" = "#E69F00", "ECEA_VIF" = "#D55E00",
    "fgsea (param=1)" = "#56B4E9", "fgsea (param=0)" = "#0072B2",
    "camera" = "#009E73", "mroast" = "#CC79A7", "GSVA" = "#F0E442"
  )

  ggplot(long_data, aes(x = tau1, y = power, color = method, group = method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2, alpha = 0.7) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    facet_grid(tau0 ~ CR, labeller = label_both) +
    scale_y_continuous(limits = c(0, 1), name = "Power (Fraction of p < 0.05)") +
    scale_color_manual(values = custom_colors) + # 使用自定义颜色
    labs(
      title = "Power Comparison of Methods",
      subtitle = final_subtitle,
      x = expression(tau[1]~"(Signal Strength)"),
      color = "Analysis Method"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 14, lineheight = 1.2),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      strip.background = element_rect(fill="grey90"),
      strip.text = element_text(face="bold")
    )
}

# ---- 函数2：绘制功效优势热力图 ----
create_advantage_heatmap <- function(data, main_prefix, main_name, competitor_prefix, competitor_name) {

  advantage_col_name <- paste0("advantage_vs_", gsub("[^A-Za-z0-9]", "_", competitor_name))

  plot_data <- data %>%
    mutate(
      !!sym(advantage_col_name) := !!sym(main_prefix) - !!sym(competitor_prefix)
    )

  ggplot(plot_data, aes(x = factor(tau1), y = factor(tau0), fill = !!sym(advantage_col_name))) +
    geom_tile(color = "black", linewidth = 0.4) +
    facet_wrap(~ CR, ncol = 2, labeller = label_bquote(bold(CR) == .(CR))) +
    scale_fill_gradient2(
      name = paste0("Power Advantage\n(", main_name, " - ", competitor_name, ")"),
      low = "#E69F00", mid = "white", high = "#0072B2",
      midpoint = 0, limits = c(-1, 1)
    ) +
    labs(
      title = paste("Power Advantage of", main_name, "vs.", competitor_name),
      x = "τ₁ (Signal Strength)",
      y = "τ₀ (Baseline Effect)"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text = element_text(face = "bold", color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}


# ---- 函数3：绘制性能径向图  ----
create_timing_plot <- function(timing_map, title_info) {

  total_times_data <- tibble(
    method = names(timing_map),
    time = unlist(timing_map)
  ) %>%
    arrange(desc(method)) %>%
    mutate(label_text = paste0(method, "\n", round(time, 1), "s"))

  # 计算 y 轴的最大值，额外增加 20% 的空间确保图形不完全闭合
  max_val <- max(total_times_data$time, na.rm = TRUE) * 1.2
  hole_size <- 2

  ggplot(
    total_times_data,
    aes(x = factor(method, levels = c(paste0("hole_", 1:hole_size), method)), y = time, fill = method)
  ) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(
      aes(y = 0, label = label_text),
      hjust = 1.1,
      color = "gray30",
      size = 3.5,
      fontface = "bold",
      lineheight = 0.8
    ) +
    coord_polar(theta = "y", clip = "off") +
    scale_y_continuous(limits = c(0, max_val)) +
    scale_x_discrete(breaks = total_times_data$method) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 80)
    ) +
    scale_fill_manual(values = c(
      "ECEA" = "#56B4E9", "ECEA_VIF" = "#D55E00",
      "fgsea (param=1)" = "#E69F00", "fgsea (param=0)" = "#56B4E9",
      "camera" = "#009E73", "mroast" = "#CC79A7", "GSVA" = "#F0E442"
    )) +
    labs(title = title_info, x = NULL, y = NULL)
}
