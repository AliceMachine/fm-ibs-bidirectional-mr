# =============================================================================
#  生成《Genetic liability to IBS and fibromyalgia》正文图表（R 版）
#  正文：Figure 1–5 + Table 1–2
#  依赖：同目录 _common.R（样式、路径、工具函数、三个分析的原始数据）
#  运行：Rscript gen_figures_tables.R
#  说明：只读数据、不改数据；所有数值从原始文件直读，无一处手打。
# =============================================================================

# 【先切 locale，再 source】否则 source 自己就会因为中文路径失败
# （报 "cannot open the connection ... Invalid argument"）。
# _common.R 里还会再设一次，重复设置无害，只是为了让两个脚本都能独立运行。
for (.loc in c("English_United States.65001", "Chinese_China.65001",
               "en_US.UTF-8", "C.UTF-8")) {
  if (suppressWarnings(tryCatch({Sys.setlocale("LC_CTYPE", .loc); TRUE},
                                error = function(e) FALSE, warning = function(w) FALSE))) {
    if (identical(Sys.getlocale("LC_CTYPE"), .loc)) break
  }
}
rm(.loc)

# 找到本脚本所在目录，source 同目录的 _common.R
# （Rscript 与 RStudio 两种运行方式的取法不同，这里都兼容）
.this_dir <- tryCatch({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
}, error = function(e) getwd())
source(file.path(.this_dir, "_common.R"))

# 本脚本额外需要的两块数据（_common.R 里没有）
coloc <- sens$colocalisation$per_iv          # 共定位 20 个区域的四后验
j_mv <- list(                                # MVMR 四个模型
  "d1-mv4" = fromJSON(file.path(BASE, D_05, "d1-mv4-mvmr.json"), simplifyVector = FALSE),
  "d1-mv3" = fromJSON(file.path(BASE, D_05, "d1-mv3-mvmr.json"), simplifyVector = FALSE),
  "d2-mv4" = fromJSON(file.path(BASE, D_05, "d2-mv4-mvmr.json"), simplifyVector = FALSE),
  "d2-mv3" = fromJSON(file.path(BASE, D_05, "d2-mv3-mvmr.json"), simplifyVector = FALSE)
)

# =============================================================================
#  Figure 1｜散点图 + 各估计器拟合线（FM → IBS，15 个工具）
# -----------------------------------------------------------------------------
#  这张图回答什么：每个工具变量自己算出的"暴露效应 vs 结局效应"落在哪儿，
#                  以及几种估计器给出的斜率是否一致（线越靠拢，结论越稳）。
#  横轴 = 工具变量对 FM 的效应；纵轴 = 同一工具变量对 IBS 的效应。
#  为什么过原点：MR 理论要求"暴露效应为 0 时结局效应也为 0"，所以不设截距。
# =============================================================================

# 一个分析的四个估计器：IVW（按预注册规则自动选固定/随机）+ Egger + 加权中位数 + 加权众数
est_block <- function(j) {
  p <- primary_ivw(j)
  r <- rbind(
    data.frame(estimator = "IVW", beta = as.numeric(p$beta), se = as.numeric(p$se),
               lo = as.numeric(p$ci_lower), hi = as.numeric(p$ci_upper),
               p = as.numeric(p$p), stringsAsFactors = FALSE),
    egger_row(j),
    est_row(j, "weighted_median", "Weighted median"),
    est_row(j, "weighted_mode",   "Weighted mode")
  )
  r[!is.na(r$beta), ]
}
est_fm <- est_block(j_d2)
est_fm$estimator <- factor(est_fm$estimator, levels = names(EST_LTY))

p1 <- ggplot(h_fm_ibs, aes(x = beta_exposure, y = beta_outcome)) +
  # 点：每个工具变量一个。点大小固定；想让点更大改 size
  geom_point(size = 2.4, colour = COL_POINT, alpha = 0.85) +
  # 误差棒：横 = 暴露侧 SE，竖 = 结局侧 SE。想去掉就删这两行 geom_errorbarh/errorbar
  geom_errorbar(aes(ymin = beta_outcome - 1.96 * se_outcome,
                    ymax = beta_outcome + 1.96 * se_outcome),
                width = 0, colour = COL_POINT, alpha = 0.35) +
  geom_errorbarh(aes(xmin = beta_exposure - 1.96 * se_exposure,
                     xmax = beta_exposure + 1.96 * se_exposure),
                 height = 0, colour = COL_POINT, alpha = 0.35) +
  # 四条拟合线：斜率 = 该估计器的 β，截距强制 0
  geom_abline(data = est_fm,
              aes(slope = beta, intercept = 0, linetype = estimator, colour = estimator),
              linewidth = 0.9) +
  scale_linetype_manual(values = EST_LTY, name = NULL) +
  scale_colour_manual(values = EST_COL, name = NULL) +
  # 参考线：x=0 与 y=0，方便看象限
  geom_hline(yintercept = 0, colour = COL_GREY, linewidth = 0.3, linetype = "dashed") +
  geom_vline(xintercept = 0, colour = COL_GREY, linewidth = 0.3, linetype = "dashed") +
  labs(
    title = "FM liability and IBS risk: per-variant effects",
    subtitle = sprintf("FM → IBS, %d instruments; line slope = MR estimate", nrow(h_fm_ibs)),
    x = "SNP effect on fibromyalgia (log-odds)",
    y = "SNP effect on irritable bowel syndrome (log-odds)"
  ) +
  theme_paper()
save_fig(p1, "Figure1_scatter_FM_IBS")

# =============================================================================
#  Figure 2｜森林图：各估计器对比（FM → IBS）
# -----------------------------------------------------------------------------
#  这张图回答什么：同一个问题用不同方法问，答案（点和区间）是否落在同一处。
#  点的位置 = β，横线 = 95% CI，竖虚线 = 0（无效线）。
#  区间跨过 0 就说明该估计器不能排除"无效应"。
# =============================================================================

# 森林图里要多两行：PRESSO 剔除离群点后的校正估计
# 注意：PRESSO 校正行用**重拟合**的 IVW 结果（refit），不是 corrected 的 b ± 1.96×se。
# 前者是"剔除离群点后重新跑一遍 IVW"，后者是用 PRESSO 自己的 se 反推区间，
# 两者数值不同（正文用的是前者）。
refit_d2 <- pres$ivw_refit_on_official_outlier_sets$d2_FM_IBS
corr_d2 <- data.frame(
  estimator = "IVW (PRESSO-corrected)",
  beta = as.numeric(refit_d2$beta),
  se   = as.numeric(refit_d2$se),
  lo   = as.numeric(refit_d2$ci[[1]]),
  hi   = as.numeric(refit_d2$ci[[2]]),
  p    = as.numeric(refit_d2$p),
  stringsAsFactors = FALSE
)

forest_d2 <- rbind(
  est_fm[, c("estimator", "beta", "se", "lo", "hi", "p")],
  corr_d2
)
# 排列顺序：从下往上画，所以反转一次让 IVW 在最上面
forest_d2$estimator <- factor(forest_d2$estimator, levels = rev(forest_d2$estimator))

p2 <- ggplot(forest_d2, aes(x = beta, y = estimator)) +
  # 无效线（β = 0）
  geom_vline(xintercept = 0, linetype = "dashed", colour = COL_GREY, linewidth = 0.4) +
  # 区间线 + 点
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, colour = COL_LINE, linewidth = 0.6) +
  geom_point(size = 2.8, colour = COL_POINT) +
  # P 值：放在最右侧的**固定** x 位置（不是跟着每个点挪），
  # 否则区间线长的那些行，文字会压在线上看不清。
  # 想让 P 值列更靠右/靠左，改下面的 1.03；想让留白更宽，改 coord_cartesian 里的 1.45。
  geom_text(aes(label = paste0("P = ", fmt_p(p))),
            x = max(forest_d2$hi) * 1.03, hjust = 0,
            size = 3.1, colour = "grey30") +
  scale_y_discrete(labels = function(x) gsub(" \\(", "\n(", x)) +  # 长标签换行
  labs(
    title = "MR estimates across estimators",
    subtitle = sprintf("FM → IBS, %d instruments", nrow(h_fm_ibs)),
    x = "Causal estimate (β, log-odds of IBS per log-odds of FM liability)",
    y = NULL
  ) +
  # x 轴向右扩到 1.45 倍，给 P 值一列专用空间（clip="off" 让超出的文字照常显示）
  coord_cartesian(xlim = c(0, max(forest_d2$hi) * 1.45), clip = "off") +
  theme_paper() +
  theme(plot.margin = margin(5.5, 12, 5.5, 5.5))   # 右 = 12pt
save_fig(p2, "Figure2_forest_FM_IBS", h = 4.2)

# =============================================================================
#  Figure 3｜散点图（IBS → FM，5 个工具）
# -----------------------------------------------------------------------------
#  与 Figure 1 同构，只是方向反过来。代码几乎一样，唯一区别是数据源。
#  只有 5 个点，所以点画大一点、透明度低一点，避免重叠加重视觉错觉。
# =============================================================================

est_ibs <- est_block(j_d1)
# Egger 的斜率在本方向极不稳定（I²_GX = 0），区间会拉到 ±4，
# 如果把它画进来整张图会被压扁，所以**图里只画三个稳健估计器**，Egger 只进表。
est_ibs_plot <- est_ibs[est_ibs$estimator != "MR-Egger", ]
est_ibs_plot$estimator <- factor(est_ibs_plot$estimator, levels = names(EST_LTY))

p3 <- ggplot(h_ibs_fm, aes(x = beta_exposure, y = beta_outcome)) +
  geom_point(size = 3.2, colour = COL_POINT, alpha = 0.8) +
  geom_errorbar(aes(ymin = beta_outcome - 1.96 * se_outcome,
                    ymax = beta_outcome + 1.96 * se_outcome),
                width = 0, colour = COL_POINT, alpha = 0.3) +
  geom_errorbarh(aes(xmin = beta_exposure - 1.96 * se_exposure,
                     xmax = beta_exposure + 1.96 * se_exposure),
                 height = 0, colour = COL_POINT, alpha = 0.3) +
  geom_abline(data = est_ibs_plot,
              aes(slope = beta, intercept = 0, linetype = estimator, colour = estimator),
              linewidth = 0.9) +
  scale_linetype_manual(values = EST_LTY, name = NULL) +
  scale_colour_manual(values = EST_COL, name = NULL) +
  geom_hline(yintercept = 0, colour = COL_GREY, linewidth = 0.3, linetype = "dashed") +
  geom_vline(xintercept = 0, colour = COL_GREY, linewidth = 0.3, linetype = "dashed") +
  # 只有 5 个点，如果不收紧范围，过原点的拟合线会把面板拉到很远，
  # 数据挤在左边一小块。这里按数据范围向外扩 8%。
  # 注意：不能用 range(x) * 1.15 —— 当数值是负数时，乘 1.15 会把右边界"往小改"，
  # 结果是把点裁掉。正确做法是按区间宽度双向加边。
  coord_cartesian(
    xlim = pad_range(h_ibs_fm$beta_exposure),
    ylim = pad_range(c(h_ibs_fm$beta_outcome - 1.96 * h_ibs_fm$se_outcome,
                       h_ibs_fm$beta_outcome + 1.96 * h_ibs_fm$se_outcome))
  ) +
  labs(
    title = "IBS liability and FM risk: per-variant effects",
    subtitle = sprintf("IBS → FM, %d instruments (MR-Egger omitted: I²_GX = 0)",
                       nrow(h_ibs_fm)),
    x = "SNP effect on irritable bowel syndrome (log-odds)",
    y = "SNP effect on fibromyalgia (log-odds)"
  ) +
  theme_paper()
save_fig(p3, "Figure3_scatter_IBS_FM")

# =============================================================================
#  Figure 4｜逐变异的 Wald 比值 + MR-PRESSO 离群点（FM → IBS）
# -----------------------------------------------------------------------------
#  这张图回答什么：哪一个工具变量"不听话"。
#  Wald 比值 = 结局效应 / 暴露效应，每个点是一个工具变量自己的因果估计。
#  离群点（红色 + 标注 rsID）就是被 PRESSO 判定偏离总体模式的位点。
#  水平虚线 = IVW 主估计，可以看离群点把整体拉高还是压低。
#  注意：离群点判定以**官方 R 包**结果为准（pres$d2_FM_IBS$outliers）。
# =============================================================================

df4 <- h_fm_ibs
df4$wald <- df4$beta_outcome / df4$beta_exposure        # 每个位点自己的因果估计
df4 <- df4[order(df4$wald), ]                            # 排序，便于观察
df4$SNP <- factor(df4$SNP, levels = df4$SNP)

outliers_d2 <- as.character(unlist(pres$d2_FM_IBS$outliers))
df4$is_outlier <- df4$SNP %in% outliers_d2              # 布尔标记，决定颜色

ivw_beta_d2 <- as.numeric(j_d2$ivw_random$beta)         # 主估计（参考线）

p4 <- ggplot(df4, aes(x = SNP, y = wald)) +
  # 参考线：IVW 主估计
  geom_hline(yintercept = ivw_beta_d2, linetype = "dashed",
             colour = COL_LINE, linewidth = 0.5) +
  # 所有位点：先画灰点作背景
  geom_point(data = df4[!df4$is_outlier, ], size = 2.4, colour = COL_POINT) +
  # 离群点：红色 + 更大
  geom_point(data = df4[df4$is_outlier, ], size = 3.4, colour = COL_OUTLIER) +
  # 给离群点贴 rsID（想让标签更靠上，改 vjust）
  geom_text(data = df4[df4$is_outlier, ],
            aes(label = SNP), vjust = -0.9, size = 3.0, colour = COL_OUTLIER) +
  labs(
    title = "Per-variant Wald ratios and MR-PRESSO outliers",
    subtitle = sprintf("FM → IBS, %d instruments; %d outlier(s); dashed line = IVW estimate",
                       nrow(df4), length(outliers_d2)),
    x = NULL,
    y = "Wald ratio (outcome effect / exposure effect)"
  ) +
  # 右侧多留 8% 空间，让最右那个离群点的 rsID 标签不被裁掉
  scale_x_discrete(expand = expansion(mult = c(0.02, 0.08))) +
  theme_paper() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = BASE_SIZE - 4))  # 横轴 rsID 斜排（角度/字号都可改）
save_fig(p4, "Figure4_presso_FM_IBS", h = 4.8)

# =============================================================================
#  Figure 5｜共定位后验概率 PP.H4（20 个工具变量区域）
# -----------------------------------------------------------------------------
#  这张图回答什么：MR 说"有关联"，共定位问"是不是同一个因果变异"。
#  横轴 = PP.H4（两性状共享同一因果变异的后验概率），0–1。
#  阈值 0.8：预注册标准，超过才算支持共定位。
#  **重要**：PP.H4 低不等于"没有共定位"——如果该区域只有一侧有信号（H1/H2 主导），
#            共定位根本不可评估。本图只展示 PP.H4，判断"可不可评估"要看附录表 S3。
# =============================================================================

# 把 20 个区域拉成一张表
df5 <- do.call(rbind, lapply(names(coloc), function(rs) {
  z <- coloc[[rs]]
  data.frame(
    iv   = rs,
    tag  = ifelse(is.null(z$tag), NA, z$tag),          # d1 = IBS→FM 的工具，d2 = FM→IBS
    H1   = as.numeric(z$PP$H1),
    H2   = as.numeric(z$PP$H2),
    H3   = as.numeric(z$PP$H3),
    H4   = as.numeric(z$PP$H4),
    stringsAsFactors = FALSE
  )
}))
# 按 PP.H4 从大到小排，最大的排最上面
df5 <- df5[order(df5$H4), ]
df5$iv <- factor(df5$iv, levels = df5$iv)
df5$pass <- ifelse(df5$H4 > 0.8, "PP.H4 > 0.8", "PP.H4 ≤ 0.8")   # 分色依据

p5 <- ggplot(df5, aes(x = H4, y = iv, fill = pass)) +
  # 条形：横向，长度 = PP.H4
  geom_col(width = 0.68) +
  # 0.8 阈值线
  geom_vline(xintercept = 0.8, linetype = "dashed",
             colour = COL_LINE, linewidth = 0.5) +
  # 只给超过阈值的那个标数值（其余数值看附录表 S3，避免图上拥挤）
  geom_text(data = df5[df5$H4 > 0.8, ],
            aes(label = sprintf("%.3f", H4)),
            hjust = -0.15, size = 3.2, colour = COL_LINE) +
  scale_fill_manual(values = c("PP.H4 > 0.8" = COL_HI, "PP.H4 ≤ 0.8" = COL_POINT),
                    name = NULL) +
  scale_x_continuous(limits = c(0, 1.05),                # 坐标范围：想放大局部改这里
                     breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Colocalisation across 20 instrument regions",
    subtitle = "Posterior probability of a shared causal variant (PP.H4)",
    x = "PP.H4 (1.0 = same causal variant)",
    y = NULL
  ) +
  theme_paper() +
  theme(axis.text.y = element_text(size = BASE_SIZE - 2),
        legend.position = "inside",
        legend.position.inside = c(0.40, 0.42))   # 右侧中部空白处：既不挡底部条形，也不压 0.8 阈值线
save_fig(p5, "Figure5_coloc", h = 5.6)

# =============================================================================
#  Table 1｜单变量 MR 估计（全部估计器，三个分析）
# -----------------------------------------------------------------------------
#  内容与正文 Table 1 完全一致：三个分析 × 各估计器 + PRESSO 校正行。
#  为什么用 gt：可以直接出 HTML（贴进 Word 或转 PDF），且能分组、能加脚注。
# =============================================================================

# 从估算器 JSON 里组装一个分析的所有行
tab1_block <- function(j, label, n_iv) {
  # 主分析是固定还是随机，由 JSON 的 primary_analysis 决定（预注册规则自动选）
  is_fixed <- identical(j$primary_analysis, "ivw_fixed")
  r <- rbind(
    data.frame(estimator = if (is_fixed) "IVW (fixed)" else "IVW (random)",
               beta = as.numeric(primary_ivw(j)$beta), se = as.numeric(primary_ivw(j)$se),
               lo = as.numeric(primary_ivw(j)$ci_lower), hi = as.numeric(primary_ivw(j)$ci_upper),
               p = as.numeric(primary_ivw(j)$p), stringsAsFactors = FALSE),
    egger_row(j, "MR-Egger (slope)"),
    est_row(j, "weighted_median", "Weighted median"),
    est_row(j, "weighted_mode",   "Weighted mode")
  )
  r$analysis <- label
  r$n_iv <- n_iv
  r[, c("analysis", "estimator", "n_iv", "beta", "lo", "hi", "p")]
}

t1 <- rbind(
  tab1_block(j_d1,   "IBS → FM", 5),
  tab1_block(j_d2,   "FM → IBS", 15),
  tab1_block(j_finn, "FM → IBS, FinnGen", 13)
)

# 追加 PRESSO 校正行（这两个数来自官方包复核，不在 univariate JSON 里）
# PRESSO 校正行：取"剔除离群点后重拟合 IVW"的 β / CI / P 与剩余位点数
rf1 <- pres$ivw_refit_on_official_outlier_sets$d1_IBS_FM
rf2 <- pres$ivw_refit_on_official_outlier_sets$d2_FM_IBS
presso_rows <- data.frame(
  analysis  = c("IBS → FM", "FM → IBS"),
  estimator = c("IVW (fixed, PRESSO-corrected)", "IVW (fixed, PRESSO-corrected)"),
  n_iv      = c(as.integer(rf1$n_iv), as.integer(rf2$n_iv)),   # 剔除离群点后的位点数
  beta      = c(as.numeric(rf1$beta), as.numeric(rf2$beta)),
  lo        = c(as.numeric(rf1$ci[[1]]), as.numeric(rf2$ci[[1]])),
  hi        = c(as.numeric(rf1$ci[[2]]), as.numeric(rf2$ci[[2]])),
  p         = c(as.numeric(rf1$p), as.numeric(rf2$p)),
  stringsAsFactors = FALSE
)

t1_all <- rbind(t1, presso_rows)
# 按分析分组 + 组内固定顺序（IVW 在前，校正行在后）
t1_all$analysis <- factor(t1_all$analysis,
                          levels = c("IBS → FM", "FM → IBS", "FM → IBS, FinnGen"))
# 组内固定顺序：IVW（不论固定/随机）→ Egger → 加权中位数 → 加权众数 → PRESSO 校正
# 用显式排名表，不用 match()——match 遇到表中没有的字符串会返回 NA，静默排到最后。
rank_map <- c("IVW (random)" = 1, "IVW (fixed)" = 1, "MR-Egger (slope)" = 2,
              "Weighted median" = 3, "Weighted mode" = 4,
              "IVW (fixed, PRESSO-corrected)" = 5)
t1_all <- t1_all[order(t1_all$analysis, rank_map[t1_all$estimator]), ]

# 加 OR 与格式化列
t1_all$or     <- exp(t1_all$beta)
t1_all$ci_txt <- sprintf("%.3f–%.3f", t1_all$lo, t1_all$hi)
t1_all$or_txt <- sprintf("%.3f", t1_all$or)
t1_all$p_txt  <- fmt_p(t1_all$p)
t1_all$label  <- sprintf("%s (%d IVs)", t1_all$analysis, t1_all$n_iv)

gt_t1 <- t1_all |>
  select(analysis, label, estimator, beta, ci_txt, or_txt, p_txt) |>
  gt(groupname_col = "analysis") |>
  cols_hide(columns = label) |>
  cols_label(estimator = "Estimator", beta = "β", ci_txt = "95% CI",
             or_txt = "OR", p_txt = "P") |>
  fmt_number(columns = beta, decimals = 3) |>
  tab_header(
    title = "Table 1. Univariable MR estimates (all estimators)",
    subtitle = "β = log-odds of the outcome per log-odds of genetic liability"
  ) |>
  tab_footnote(
    footnote = "PRESSO-corrected rows refit IVW (fixed effects) on variants remaining after outlier removal.",
    locations = cells_body(columns = estimator, rows = grepl("PRESSO", estimator))
  ) |>
  tab_options(table.font.size = px(12), heading.title.font.size = px(14))
gtsave(gt_t1, file.path(TAB_OUT, "Table1_univariable.html"))

# 同时导一份 CSV 便于核对（不含格式）
write.csv(t1_all[, c("analysis", "estimator", "n_iv", "beta", "lo", "hi", "or", "p")],
          file.path(TAB_OUT, "Table1_univariable.csv"), row.names = FALSE)

# =============================================================================
#  Table 2｜多变量 MR 净效应（探索性）
# -----------------------------------------------------------------------------
#  四行 = 两个方向 × 两种协变量组合（4 协变量 / 2 协变量）。
#  P 值不在 JSON 里，按 β/SE 用正态近似现算（与正文方法一致）。
#  **这张表是探索性的**：模型几乎全靠协变量工具识别，标题里必须带 exploratory。
# =============================================================================

mv_row <- function(j, model_label, direction_label) {
  core <- j$exposures[[1]]                     # 第一个暴露 = 主暴露（注意是 [[ ]]，JSON 里是 list）
  e <- j$estimates[[core]]
  beta <- as.numeric(e$beta); se <- as.numeric(e$se)
  data.frame(
    model     = model_label,
    direction = direction_label,
    exposure  = core,
    variants  = as.integer(j$n_snps),
    beta      = beta,
    ci_txt    = sprintf("%.3f–%.3f", as.numeric(e$ci_lower), as.numeric(e$ci_upper)),
    or_txt    = sprintf("%.3f", exp(beta)),
    p_txt     = fmt_p(2 * pnorm(-abs(beta / se))),   # 双侧正态近似
    stringsAsFactors = FALSE
  )
}

t2 <- rbind(
  mv_row(j_mv[["d1-mv4"]], "4-covariate (BMI + depression + insomnia)", "IBS → FM"),
  mv_row(j_mv[["d2-mv4"]], "4-covariate (BMI + depression + insomnia)", "FM → IBS"),
  mv_row(j_mv[["d1-mv3"]], "2-covariate (BMI + depression)",            "IBS → FM"),
  mv_row(j_mv[["d2-mv3"]], "2-covariate (BMI + depression)",            "FM → IBS")
)

gt_t2 <- t2 |>
  select(model, direction, exposure, variants, beta, ci_txt, or_txt, p_txt) |>
  gt(groupname_col = "model") |>
  cols_label(direction = "Exposure (net effect)", exposure = "Term",
             variants = "Variants", beta = "β", ci_txt = "95% CI",
             or_txt = "OR", p_txt = "P") |>
  fmt_number(columns = beta, decimals = 3) |>
  tab_header(
    title = "Table 2. Multivariable MR estimates (exploratory)",
    subtitle = "Net effect after simultaneous adjustment; see section 3.4 for identification limits"
  ) |>
  tab_footnote(
    footnote = "P values are two-sided normal approximations of β/SE.",
    locations = cells_column_labels(columns = p_txt)
  ) |>
  tab_options(table.font.size = px(12), heading.title.font.size = px(14))
gtsave(gt_t2, file.path(TAB_OUT, "Table2_mvmr.html"))
write.csv(t2, file.path(TAB_OUT, "Table2_mvmr.csv"), row.names = FALSE)

# =============================================================================
#  补充图（S1–S6）怎么做：复用上面同一套函数
# -----------------------------------------------------------------------------
#  正文 5 图之外，补充材料还有 6 张。它们不是新类型的图，只是同一套图的变体：
#
#    S1 留一法三联  → 每个 SNP 剔除一次后的 β，用 geom_pointrange 画，三分析并排
#    S2 漏斗图三联  → x = β/se，y = 1/se，加垂直线 = IVW 估计
#    S3 PRESSO（IBS→FM）→ 复制 Figure 4 的代码，把 df4 换成 h_ibs_fm、离群换 d1 的
#    S4 单 SNP 森林  → 逐 SNP 的 β（不是 Wald 比值），横线用 delta 法近似 SE
#    S5 森林图（IBS→FM）→ 复制 Figure 2，把 est_fm 换成 est_ibs
#    S6 工具变量流程 → 818 → 749 → 15 的漏斗，用 geom_rect + geom_text 手绘
#
#  下面给一个"留一法"的可复用函数，说明扩展方式（其余各图同理）：
#
#  plot_loo <- function(j, title) {
#    loo <- do.call(rbind, lapply(j$leave_one_out, function(z)
#      data.frame(index = z$index, beta = z$beta, se = z$se)))
#    loo$lo <- loo$beta - 1.96 * loo$se
#    loo$hi <- loo$beta + 1.96 * loo$se
#    loo$snp <- j$snps[loo$index + 1]            # index 从 0 开始，R 从 1 开始
#    ggplot(loo, aes(x = beta, y = reorder(snp, beta))) +
#      geom_vline(xintercept = as.numeric(j$ivw_random$beta), linetype = "dashed") +
#      geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
#      geom_point() + labs(title = title, x = "β (excluding one SNP)", y = NULL) +
#      theme_paper()
#  }
#  save_fig(plot_loo(j_d1, "Leave-one-out: IBS → FM"), "SuppS1a_loo_IBS_FM", h = 3.6)
#
#  三联图（S1/S2）的拼版用 patchwork：
#    library(patchwork)
#    (plot_loo(j_d1, "A  IBS → FM") | ... ) / ( ... ) + plot_layout(heights = c(1,1,1))
# =============================================================================

message("\n全部完成。图 → ", FIG_OUT, "\n           表 → ", TAB_OUT)
