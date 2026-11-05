# =============================================================================
#  生成补充材料图 S1–S6（R 版）
#  S1 留一法三联 / S2 漏斗图三联 / S3 PRESSO（IBS→FM）
#  S4 单 SNP 森林（FM→IBS）/ S5 森林图（IBS→FM）/ S6 工具变量筛选流程
#  依赖：同目录 _common.R
#  运行：Rscript gen_supp_figures.R
# =============================================================================

# ---- 0. 先切 locale（否则 source 中文路径会失败），再载入共用配置 ----------
for (.loc in c("English_United States.65001", "Chinese_China.65001",
               "en_US.UTF-8", "C.UTF-8")) {
  if (suppressWarnings(tryCatch({Sys.setlocale("LC_CTYPE", .loc); TRUE},
                                error = function(e) FALSE, warning = function(w) FALSE))) {
    if (identical(Sys.getlocale("LC_CTYPE"), .loc)) break
  }
}
rm(.loc)

.this_dir <- tryCatch({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
}, error = function(e) getwd())
source(file.path(.this_dir, "_common.R"))

suppressPackageStartupMessages(library(patchwork))   # 三联图拼版用

# 三个分析的统一标签（三个补充图都用同一套，避免一处改了另一处没改）
LABELS <- c(
  d1   = "A  IBS → FM (5 instruments)",
  d2   = "B  FM → IBS (15 instruments)",
  finn = "C  FM → IBS, FinnGen (13 instruments)"
)

# 把三个分析打包成一个列表，后面循环画
ANALYSES <- list(
  d1   = list(j = j_d1,   h = h_ibs_fm, key = "d1"),
  d2   = list(j = j_d2,   h = h_fm_ibs, key = "d2"),
  finn = list(j = j_finn, h = h_finn,   key = "finn")
)

# =============================================================================
#  S1｜留一法（leave-one-out）：剔掉任意一个位点，结论会不会翻
# -----------------------------------------------------------------------------
#  这张图回答什么：把 15 个工具变量轮流剔掉一个重算，如果每个结果都落在
#                  同一条方向、同一个量级上，说明结论不是被某一个位点"绑架"的。
#  横轴 = 剔除该位点后的合并估计 β；纵轴 = 被剔掉的那个位点。
#  竖虚线 = 全量分析的 β（所有点都不该偏离它太远）。
#  三个分析竖排成 A/B/C 三个面板，合成一张 Supplementary Figure S1。
# =============================================================================

plot_loo <- function(entry, title) {
  j <- entry$j
  loo <- do.call(rbind, lapply(j$leave_one_out, function(z)
    data.frame(index = as.integer(z$index),
               beta  = as.numeric(z$beta),
               se    = as.numeric(z$se))))
  # JSON 的 index 从 0 开始，R 从 1 开始，所以取名时要 +1
  loo$snp <- unlist(j$snps)[loo$index + 1]   # j$snps 是 list，必须 unlist 才能当字符向量用
  loo$lo  <- loo$beta - 1.96 * loo$se
  loo$hi  <- loo$beta + 1.96 * loo$se

  full_beta <- as.numeric(primary_ivw(j)$beta)   # 全量分析的 β（参考线）

  ggplot(loo, aes(x = beta, y = reorder(snp, beta))) +
    geom_vline(xintercept = full_beta, linetype = "dashed",
               colour = COL_LINE, linewidth = 0.5) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25,
                   colour = COL_POINT, linewidth = 0.5, alpha = 0.8) +
    geom_point(size = 2.2, colour = COL_POINT) +
    labs(title = title, x = "β after excluding the variant", y = NULL) +
    theme_paper() +
    theme(
      plot.title = element_text(face = "bold", size = BASE_SIZE - 1),
      axis.text.y = element_text(size = BASE_SIZE - 3)   # 位点名多，字号调小
    )
}

s1 <- (plot_loo(ANALYSES$d1,   LABELS[["d1"]]) /
       plot_loo(ANALYSES$d2,   LABELS[["d2"]]) /
       plot_loo(ANALYSES$finn, LABELS[["finn"]])) +
  plot_layout(heights = c(1, 1.7, 1.5))    # d2 有 15 个位点，给它更高的行
save_fig(s1, "FigS1_leave_one_out", w = 7.2, h = 10.5, dir = SUP_OUT)

# =============================================================================
#  S2｜漏斗图：横轴是估计值，纵轴是精度
# -----------------------------------------------------------------------------
#  这张图回答什么：精度低的位点（图下方）散得开、精度高的（上方）该聚拢。
#                  如果明显不对称，可能提示异质性或多效性 —— **但不能直接
#                  当成发表偏倚**（MR 不是发表研究，图注里必须写这句）。
#  横轴 = 该位点的 Wald 比值（结局效应 / 暴露效应）；纵轴 = 1 / SE（精度）。
#  ★ 图注固定句：Asymmetry may reflect heterogeneity or pleiotropy rather than
#    publication bias.（正文与补充说明都用这一句，改的时候两边一起改）
# =============================================================================

plot_funnel <- function(entry, title) {
  j <- entry$j
  h <- entry$h
  d <- data.frame(
    ratio = h$beta_outcome / h$beta_exposure,
    prec  = 1 / h$se_outcome
  )
  est <- as.numeric(primary_ivw(j)$beta)

  ggplot(d, aes(x = ratio, y = prec)) +
    geom_vline(xintercept = est, linetype = "dashed",
               colour = COL_LINE, linewidth = 0.5) +
    geom_point(size = 2.3, colour = COL_POINT, alpha = 0.85) +
    labs(title = title,
         x = "Wald ratio (outcome effect / exposure effect)",
         y = "Precision (1 / SE of outcome effect)") +
    theme_paper() +
    theme(plot.title = element_text(face = "bold", size = BASE_SIZE - 1))
}

s2 <- (plot_funnel(ANALYSES$d1,   LABELS[["d1"]]) /
       plot_funnel(ANALYSES$d2,   LABELS[["d2"]]) /
       plot_funnel(ANALYSES$finn, LABELS[["finn"]])) +
  plot_layout(heights = c(1, 1, 1))
save_fig(s2, "FigS2_funnel", w = 7.2, h = 9.5, dir = SUP_OUT)

# =============================================================================
#  S3｜PRESSO 逐 SNP Wald 比值（IBS → FM）
# -----------------------------------------------------------------------------
#  与正文 Figure 4 同一套画法，只是换成方向一（5 个工具）与它的离群点。
#  离群点判定用官方 R 包（pres$d1_IBS_FM$outliers），不用自建实现。
#  只有 5 个点，横轴标签不必斜排。
# =============================================================================

df_s3 <- h_ibs_fm
df_s3$wald <- df_s3$beta_outcome / df_s3$beta_exposure
df_s3 <- df_s3[order(df_s3$wald), ]
df_s3$SNP <- factor(df_s3$SNP, levels = df_s3$SNP)

out_s3 <- as.character(unlist(pres$d1_IBS_FM$outliers))
df_s3$is_outlier <- df_s3$SNP %in% out_s3
ivw_s3 <- as.numeric(primary_ivw(j_d1)$beta)

s3 <- ggplot(df_s3, aes(x = SNP, y = wald)) +
  geom_hline(yintercept = ivw_s3, linetype = "dashed",
             colour = COL_LINE, linewidth = 0.5) +
  geom_point(data = df_s3[!df_s3$is_outlier, ], size = 3, colour = COL_POINT) +
  geom_point(data = df_s3[df_s3$is_outlier, ], size = 4, colour = COL_OUTLIER) +
  geom_text(data = df_s3[df_s3$is_outlier, ],
            aes(label = SNP), vjust = -1.1, size = 3.4, colour = COL_OUTLIER) +
  # 只有 5 个点，右侧留白不用太多
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.12))) +
  labs(
    title = "Per-variant Wald ratios and MR-PRESSO outliers",
    subtitle = sprintf("IBS → FM, %d instruments; %d outlier(s); dashed line = IVW estimate",
                       nrow(df_s3), length(out_s3)),
    x = NULL,
    y = "Wald ratio (outcome effect / exposure effect)"
  ) +
  theme_paper()
save_fig(s3, "FigS3_presso_IBS_FM", w = 6.4, h = 4.6, dir = SUP_OUT)

# =============================================================================
#  S4｜单 SNP 森林图（FM → IBS）：每个位点自己说了什么
# -----------------------------------------------------------------------------
#  和 Figure 4 的区别：Figure 4 画的是"Wald 比值"（无区间）；
#  这张图画的是**每个位点的 β 与 95% CI**，可以看出各位点的估计有多分散
#  （异质性的直观来源）。
#  SE 用一阶 delta 法近似：se_ratio ≈ se_outcome / |beta_exposure|。
#  ★ 图注必须声明这个近似（暴露侧 F>30 时 se_exposure 可忽略）。
# =============================================================================

df_s4 <- h_fm_ibs
df_s4$beta_snp <- df_s4$beta_outcome / df_s4$beta_exposure
df_s4$se_snp   <- df_s4$se_outcome / abs(df_s4$beta_exposure)   # delta 法近似
df_s4$lo <- df_s4$beta_snp - 1.96 * df_s4$se_snp
df_s4$hi <- df_s4$beta_snp + 1.96 * df_s4$se_snp
df_s4 <- df_s4[order(df_s4$beta_snp), ]
df_s4$SNP <- factor(df_s4$SNP, levels = df_s4$SNP)

peto_s4 <- as.numeric(primary_ivw(j_d2)$beta)

s4 <- ggplot(df_s4, aes(x = beta_snp, y = SNP)) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = COL_GREY, linewidth = 0.4) +
  geom_vline(xintercept = peto_s4, linetype = "dashed",
             colour = COL_LINE, linewidth = 0.5) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25,
                 colour = COL_POINT, linewidth = 0.5, alpha = 0.85) +
  geom_point(size = 2.2, colour = COL_POINT) +
  labs(
    title = "Single-variant estimates for FM → IBS",
    subtitle = "Each row is one instrument; dashed line = pooled IVW estimate",
    x = "Causal estimate (β, log-odds of IBS per log-odds of FM liability)",
    y = NULL
  ) +
  theme_paper() +
  theme(axis.text.y = element_text(size = BASE_SIZE - 3))
save_fig(s4, "FigS4_singlesnp_FM_IBS", w = 7.2, h = 5.6, dir = SUP_OUT)

# =============================================================================
#  S5｜森林图（IBS → FM）
# -----------------------------------------------------------------------------
#  与正文 Figure 2 同一套画法，换成方向一。
#  注意：MR-Egger 在本方向的斜率区间是 −3.69～4.82（I²_GX = 0），
#  画进来会把整张图的 x 轴压扁，所以**这里也不画 Egger**，只在表里报。
#  如果要画（例如审稿人要求），把下面 filter 那一行删掉即可。
# =============================================================================

est_s5 <- est_block(j_d1)
est_s5 <- est_s5[est_s5$estimator != "MR-Egger", ]     # ← 想画 Egger 就删这行
est_s5$estimator <- factor(est_s5$estimator, levels = rev(est_s5$estimator))

s5 <- ggplot(est_s5, aes(x = beta, y = estimator)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = COL_GREY, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18,
                 colour = COL_LINE, linewidth = 0.6) +
  geom_point(size = 2.8, colour = COL_POINT) +
  geom_text(aes(label = paste0("P = ", fmt_p(p))),
            x = max(est_s5$hi) * 1.03, hjust = 0, size = 3.1, colour = "grey30") +
  coord_cartesian(xlim = c(0, max(est_s5$hi) * 1.45), clip = "off") +
  labs(
    title = "MR estimates across estimators",
    subtitle = sprintf("IBS → FM, %d instruments (MR-Egger shown in Table 1 only)",
                       nrow(h_ibs_fm)),
    x = "Causal estimate (β, log-odds of FM per log-odds of IBS liability)",
    y = NULL
  ) +
  theme_paper() +
  theme(plot.margin = margin(5.5, 12, 5.5, 5.5))
save_fig(s5, "FigS5_forest_IBS_FM", w = 7.2, h = 4.0, dir = SUP_OUT)

# =============================================================================
#  S6｜工具变量筛选流程（818 → 749 → 15）
# -----------------------------------------------------------------------------
#  这张图回答什么：从"全基因组显著位点"到"最终进模型的独立工具"，
#                  中间每一步掉了多少、为什么掉。
#  数字不手打：818 与 749 从 04-数据/fm_download/clump_log.md 读，
#  15 从结果 JSON 的 n_snp 读，并做一致性断言（对不上就报错停下）。
# =============================================================================

flow_file <- file.path(BASE, D_04, "fm_download", "clump_log.md")
flow_txt  <- readLines(flow_file, warn = FALSE, encoding = "UTF-8")

n_sig <- as.integer(sub(".*\\*\\*(\\d+)\\*\\*.*", "\\1",
                        grep("p<5e-8", flow_txt, value = TRUE)[1]))
n_ld  <- as.integer(sub(".*\\*\\*(\\d+)/\\d+\\*\\*.*", "\\1",
                        grep("LD 参考面板覆盖", flow_txt, value = TRUE)[1]))
n_iv  <- as.integer(j_d2$n_snp)

# 三处断言：任何一个对不上就立刻停下，避免"静默用错数"
stopifnot(!is.na(n_sig), n_sig == 818L)
stopifnot(!is.na(n_ld),  n_ld  == 749L)
stopifnot(!is.na(n_iv),  n_iv  == 15L)

# 用矩形 + 文字 + 箭头手绘流程图（不依赖 DiagrammeR 之类的额外包）
boxes <- data.frame(
  xmin = 1.6, xmax = 6.4,
  ymin = c(7.4, 4.4, 1.4),
  ymax = c(8.9, 5.9, 2.9),
  label = c(
    sprintf("Genome-wide significant variants\n(P < 5 × 10⁻⁸):  n = %d", n_sig),
    sprintf("With LD reference-panel coverage\nafter two rounds of clumping:  n = %d", n_ld),
    sprintf("Independent instruments\nentering the MR models:  n = %d", n_iv)
  )
)
# 右侧标注每一步"掉了多少、为什么掉"
notes <- data.frame(
  x = 6.7,
  y = c(6.65, 3.65),
  label = c(
    sprintf("−%d  no LD panel coverage", n_sig - n_ld),
    sprintf("−%d  merged as the same signal\n(r² > 0.001 within 10,000 kb)", n_ld - n_iv)
  )
)

s6 <- ggplot() +
  geom_rect(data = boxes,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "white", colour = COL_LINE, linewidth = 0.5) +
  geom_text(data = boxes, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
            size = 3.5, lineheight = 1.15, colour = "grey15") +
  # 向下箭头（三段框之间）
  geom_segment(aes(x = 4, xend = 4, y = 7.25, yend = 6.05),
               arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
               colour = COL_LINE, linewidth = 0.6) +
  geom_segment(aes(x = 4, xend = 4, y = 4.25, yend = 3.05),
               arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
               colour = COL_LINE, linewidth = 0.6) +
  geom_text(data = notes, aes(x = x, y = y, label = label),
            hjust = 0, size = 3.0, colour = COL_OUTLIER, lineheight = 1.15) +
  coord_cartesian(xlim = c(1.2, 11), ylim = c(1, 9.6), clip = "off") +
  labs(title = "Selection of instrumental variables for fibromyalgia",
       subtitle = "From genome-wide significant variants to independent instruments") +
  theme_void(base_size = BASE_SIZE) +
  theme(
    plot.title = element_text(face = "bold", size = BASE_SIZE + 1),
    plot.subtitle = element_text(colour = "grey35", size = BASE_SIZE - 1),
    plot.margin = margin(10, 10, 10, 10)
  )
save_fig(s6, "FigS6_iv_flow", w = 7.6, h = 5.4, dir = SUP_OUT)

# =============================================================================
message("\n补充图全部完成 → ", SUP_OUT)
# =============================================================================
