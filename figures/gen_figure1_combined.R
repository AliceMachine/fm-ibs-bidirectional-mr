# =============================================================================
#  为投 Clinical Gastroenterology and Hepatology（CGH）合并两张散点图
#
#  背景：CGH 规定"图 + 表合计不超过 6 个"。本稿原有 5 图 + 2 表 = 7 个，
#        超出 1 个。CGH 同时允许"每张图最多 6 个 panel（A–F）"，
#        因此把 Figure 1（FM → IBS 散点）与 Figure 3（IBS → FM 散点）
#        拼成一张双 panel 图，图从 5 张减到 4 张，信息一点不删。
#
#  做法：source 正文作图脚本，复用它算好的 p1 与 p3（保证与正文图完全一致），
#        只改标题文字并把两张并排，不改任何数据、不改任何拟合线。
#
#  运行：Rscript gen_figure1_combined.R
#  输出：figures/Figure1_scatter_both_directions.pdf / .png（300 dpi）
# =============================================================================

# 【先切 locale，再 source】否则 source 会因为中文路径失败
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

# 复用正文作图脚本：它会算出 p1（FM→IBS 散点）与 p3（IBS→FM 散点），
# 并顺带把 Figure1–5 全部重存一遍（内容一致，无副作用）。
source(file.path(.this_dir, "gen_figures_tables.R"))

suppressPackageStartupMessages(library(patchwork))

# 两张图原本各带一条长标题（"…: per-variant effects" + "line slope = MR estimate"），
# 并排后会换行挤在一起，所以这里只改标题文字，不动图形本身。
p_A <- p1 + labs(
  title = "Fibromyalgia liability and IBS risk",
  subtitle = sprintf("FM \u2192 IBS, %d instruments; line slope = MR estimate",
                     nrow(h_fm_ibs))
)
p_B <- p3 + labs(
  title = "IBS liability and fibromyalgia risk",
  # 副标题不能太长，否则会被面板右边界截断（MR-Egger 不可估计一事放图注里说明）
  subtitle = sprintf("IBS \u2192 FM, %d instruments", nrow(h_ibs_fm))
)

# 关键一步：两张图的估计器个数不同（右图没有 MR-Egger），
# 图例的 scale 不一致时 patchwork 无法合并图例，会各画一套。
# 所以这里把两张图的图例刻度统一成同一组 4 个估计器（drop = FALSE 保留空档），
# 合并后才只有一个图例。右图缺的那条线由副标题注明"not estimable"。
LEV4 <- names(EST_LTY)
for (.p in c("p_A", "p_B")) {
  assign(.p, suppressWarnings(
    get(.p) +
      scale_linetype_manual(values = EST_LTY, name = NULL, limits = LEV4, drop = FALSE) +
      scale_colour_manual(values = EST_COL, name = NULL, limits = LEV4, drop = FALSE)
  ))
}
rm(.p)

# 图例只保留一套（放在整图左下角 panel 下方）。
# 说明：patchwork 的 guides = "collect" 在本机 1.3.1 上对"两图估计器个数不同"
# 的情形不生效，会各画一套且右侧被裁，所以改为显式只留一套。
p_A <- p_A + theme(
  legend.position   = "bottom",
  legend.direction  = "horizontal",
  legend.text       = element_text(size = 9),
  legend.key.width  = unit(0.55, "cm"),
  legend.spacing.x  = unit(0.12, "cm")
)
p_B <- p_B + theme(legend.position = "none")

p_comb <- (p_A | p_B) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 15))
  )

# 9.0 × 4.4 英寸：CGH 双栏图印刷宽度约 7 英寸，按 78% 缩放后
# 正文 11 pt 约剩 8.6 pt，仍在可读区间；再窄轴标签会被压。
save_fig(p_comb, "Figure1_scatter_both_directions", w = 9.0, h = 4.4)

message("完成：Figure1_scatter_both_directions（A = FM \u2192 IBS，B = IBS \u2192 FM）")
