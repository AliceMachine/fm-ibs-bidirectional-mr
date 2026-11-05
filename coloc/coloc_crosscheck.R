# -*- coding: utf-8 -*-
# R coloc 官方包交叉验证：对 3 个代表区域独立复算 coloc.abf
# 运行方式：工作目录 = coloc-r-crosscheck/，全部相对路径（避开中文绝对路径坑）
# 口径：type="cc"，p1=p2=1e-4，p12=1e-5 —— 与 Python 引擎完全一致
#       （R 源码 approx.bf.estimates：cc 时 sd.prior=0.2 → W=0.04，V=varbeta）

library(coloc)

REGIONS <- c("rs2587363", "rs10156602", "rs426444")

ref <- read.csv("python_reference.csv", stringsAsFactors = FALSE)
ref <- ref[match(REGIONS, ref$rsid), ]

rows <- data.frame(rsid = character(), tag = character(), n_snp = integer(),
                   PP.H0 = numeric(), PP.H1 = numeric(), PP.H2 = numeric(),
                   PP.H3 = numeric(), PP.H4 = numeric(),
                   py.H4 = numeric(), abs.diff = numeric(), stringsAsFactors = FALSE)

for (i in seq_along(REGIONS)) {
  rs <- REGIONS[i]
  fm <- read.csv(paste0(rs, "_FM.csv"), stringsAsFactors = FALSE)
  ibs <- read.csv(paste0(rs, "_IBS.csv"), stringsAsFactors = FALSE)

  d1 <- list(beta = fm$beta, varbeta = fm$se^2, snp = fm$rsid, type = "cc")
  d2 <- list(beta = ibs$beta, varbeta = ibs$se^2, snp = ibs$rsid, type = "cc")

  res <- coloc.abf(dataset1 = d1, dataset2 = d2,
                   p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  s <- res$summary
  pp <- as.numeric(s[grep("^PP.H", names(s))])   # 按 H0..H4 顺序
  nsnps <- as.integer(s[["nsnps"]])              # summary 是原子向量，须用 [[ ]]

  pyH4 <- ref$py_H4[i]
  rH4 <- pp[5]
  rows[nrow(rows) + 1, ] <- list(rs, ref$tag[i], nsnps,
                                 pp[1], pp[2], pp[3], pp[4], pp[5],
                                 pyH4, abs(rH4 - pyH4))

  cat("=====", rs, paste0("(", ref$tag[i], ")"), "=====\n")
  cat("R  coloc nsnps =", nsnps, " | PP.H4 =", sprintf("%.6f", rH4), "\n")
  cat("Python     n =", ref$n_snp[i], " | PP.H4 =", sprintf("%.6f", pyH4),
      " | |diff| =", sprintf("%.2e", abs(rH4 - pyH4)), "\n")
}

write.csv(rows, "r_vs_python.csv", row.names = FALSE)
cat("\nsaved: r_vs_python.csv\n")
maxdiff <- max(rows$abs.diff)
cat("MAX |PP.H4(R) - PP.H4(Python)| =", sprintf("%.3e", maxdiff), "\n")
if (maxdiff < 1e-3) cat("CROSSCHECK_PASS\n") else cat("CROSSCHECK_REVIEW_NEEDED\n")
