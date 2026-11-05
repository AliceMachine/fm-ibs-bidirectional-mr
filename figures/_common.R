# =============================================================================
#  _common.R —— 正文图与补充图共用的样式、路径与工具函数
#  被 gen_figures_tables.R 与 gen_supp_figures.R 同时 source。
#  改配色 / 字号 / 输出尺寸，只改这一个文件，两张脚本立刻同步。
# =============================================================================

# ---- 0. 环境 ---------------------------------------------------------------
# 【别删】Windows 上 R 启动会尝试 C.UTF-8（系统没有），退回 C locale 后
# 含中文的路径既读不到也写不出（"unable to translate ... to native encoding"）。
# 切到 UTF-8 代码页（65001）后，中文路径与中文注释都正常。按顺序试，哪个成功用哪个。
for (.loc in c("English_United States.65001",
               "Chinese_China.65001",
               "en_US.UTF-8",
               "C.UTF-8")) {
  if (suppressWarnings(tryCatch({Sys.setlocale("LC_CTYPE", .loc); TRUE},
                                error = function(e) FALSE,
                                warning = function(w) FALSE))) {
    if (identical(Sys.getlocale("LC_CTYPE"), .loc)) break
  }
}
rm(.loc)

# ---- 1. 依赖 ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
  library(dplyr)
  library(scales)
  library(gt)
})

# ---- 2. 路径（只改 D_ROOT 一处就能换电脑） ----------------------------------
D_ROOT  <- "D:/2027/workbuddy working/021"   # ← 改这里
D_PROJ  <- "MR演练"
D_04    <- "04-数据"
D_05    <- "05-分析"
D_06    <- "06-敏感性"
D_07    <- "07-图表"

BASE    <- file.path(D_ROOT, D_PROJ)
FIG_OUT <- file.path(BASE, D_07, "R", "figures")   # 图输出
TAB_OUT <- file.path(BASE, D_07, "R", "tables")    # 表输出
SUP_OUT <- file.path(BASE, D_07, "R", "figures_supp")  # 补充图输出
dir.create(FIG_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(SUP_OUT, recursive = TRUE, showWarnings = FALSE)

# ---- 3. 全局样式（改风格只改这一段） ----------------------------------------
COL_POINT   <- "#4C72B0"   # 散点 / 主色（蓝）
COL_OUTLIER <- "#C44E52"   # 离群点 / 强调（红）
COL_LINE    <- "#333333"   # 拟合线 / 参考线（深灰）
COL_HI      <- "#DD8452"   # 高亮（橙，PP.H4 > 0.8）
COL_GREY    <- "#999999"   # 次要文字 / 网格
COL_ALT     <- "#55A868"   # 第二个分类色（补充图的分面用）

EST_LTY <- c(
  "IVW"             = "solid",
  "MR-Egger"        = "dashed",
  "Weighted median" = "dotdash",
  "Weighted mode"   = "dotted"
)
EST_COL <- c(
  "IVW"             = COL_LINE,
  "MR-Egger"        = "#8172B3",
  "Weighted median" = "#55A868",
  "Weighted mode"   = "#937860"
)

BASE_SIZE  <- 11    # 全局基准字号
FIG_WIDTH  <- 7.2   # 单栏约 3.5、双栏约 7.2 英寸
FIG_HEIGHT <- 5.4
DPI        <- 300

theme_paper <- function(base_size = BASE_SIZE) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      axis.text  = element_text(colour = "grey20"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(colour = "grey35", size = base_size - 1),
      legend.position = "inside",
      legend.position.inside = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
      legend.key.width = unit(1.6, "lines")
    )
}

# 存图：PDF（矢量投稿用）+ PNG（预览用），输出目录可选
save_fig <- function(p, name, w = FIG_WIDTH, h = FIG_HEIGHT, dir = FIG_OUT) {
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = DPI)
  message("  saved: ", name, ".pdf / .png")
}

# ---- 4. 工具函数 -----------------------------------------------------------

# 从 JSON 里安全取一个数
jnum <- function(x, path) {
  v <- tryCatch(Reduce(`[[`, strsplit(path, "/")[[1]], x), error = function(e) NULL)
  if (is.null(v)) NA_real_ else as.numeric(v)
}

# 读一个估计器块 → 一行 data.frame
est_row <- function(j, key, label) {
  e <- j[[key]]
  if (is.null(e) || is.null(e$beta)) return(NULL)
  data.frame(estimator = label, beta = as.numeric(e$beta), se = as.numeric(e$se),
             lo = as.numeric(e$ci_lower), hi = as.numeric(e$ci_upper),
             p = as.numeric(e$p), stringsAsFactors = FALSE)
}

# MR-Egger 的 β 在 mr_egger$slope 下面，不在顶层（直接取 $beta 会得到 NULL，
# rbind 时静默丢行 —— 图和表就都会少一条）
egger_row <- function(j, label = "MR-Egger") {
  s <- j$mr_egger$slope
  if (is.null(s) || is.null(s$beta)) return(NULL)
  data.frame(estimator = label, beta = as.numeric(s$beta), se = as.numeric(s$se),
             lo = as.numeric(s$ci_lower), hi = as.numeric(s$ci_upper),
             p = as.numeric(s$p), stringsAsFactors = FALSE)
}

# 主分析用哪个 IVW：JSON 的 primary_analysis 写着 ivw_random 或 ivw_fixed，
# 这是预注册规则自动选的，不能统一写死 random（FinnGen 复制是 fixed）
primary_ivw <- function(j) j[[j$primary_analysis]]

# 一个分析的四个估计器
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

# P 值排版
fmt_p <- function(p) {
  ifelse(is.na(p), "—",
         ifelse(p < 0.001,
                formatC(p, format = "e", digits = 1),
                formatC(p, format = "f", digits = 3)))
}

# 读 harmonised.csv
read_harm <- function(sub) {
  f <- file.path(BASE, D_04, "harmonise", sub, "harmonised.csv")
  read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
}

# 按区间宽度向外扩（用于收紧坐标）
#
# 注意不能用 range(x) * 1.15 —— 数值为负时乘 1.15 会把右边界"往小改"，
# 结果把点裁掉（本项目踩过：Figure 3 一度只剩 2 个点）。
pad_range <- function(x, pct = 0.08) {
  r <- range(x, na.rm = TRUE)
  r + c(-1, 1) * diff(r) * pct
}

# ---- 5. 数据（两个脚本都要用） ----------------------------------------------
h_fm_ibs <- read_harm("FM-IBS")            # 15 工具：FM → IBS
h_ibs_fm <- read_harm("IBS-FM")            #  5 工具：IBS → FM
h_finn   <- read_harm("FM-IBS-finnrep")    # 13 工具：FinnGen 复制

j_d2   <- fromJSON(file.path(BASE, D_05, "d2-FM-IBS-univariate.json"), simplifyVector = FALSE)
j_d1   <- fromJSON(file.path(BASE, D_05, "d1-IBS-FM-univariate.json"), simplifyVector = FALSE)
j_finn <- fromJSON(file.path(BASE, D_05, "d2-finnrep-FM-IBS.json"), simplifyVector = FALSE)

pres <- fromJSON(file.path(BASE, D_06, "presso-official-10k.json"), simplifyVector = FALSE)
sens <- fromJSON(file.path(BASE, D_06, "sensitivity.json"), simplifyVector = FALSE)
