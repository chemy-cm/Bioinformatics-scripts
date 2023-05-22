#!/usr/bin/Rscript

library(ggplot2)
library(ggthemes)

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
output_file <- args[2]

# 读取数据
dt <- read.table(file = input_file, sep = "\t", header = TRUE)

# 设置颜色向量
col <- rep(c("#E41A1C", "#377EB8"), length.out = length(unique(dt$CHROM)))

# 绘制曼哈顿图
popPlot <- ggplot(dt) +
  geom_line(aes(BIN_START, SNP_COUNT, color = CHROM), size = 0.8, alpha = 1) +
  theme_classic() +  # 使用 theme_classic() 函数来获得经典的主题
  labs(x = "Chromosome", y = "No. of SNPs/1kbp") +
  scale_color_manual(values = col) +
  scale_x_continuous(
    limits = c(-100000, 9200000),
    breaks = c(437623, 1371267, 2418803.5, 3524450.5, 4627710, 5829966.5, 7138210, 8460556.5),
    labels = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, ceiling(60* 1.1)),  # 将纵坐标最大值扩大10%的空间
    expand = c(0, 0),
    breaks = seq(0, ceiling(max(dt$SNP_COUNT) * 1.1), by = 5),  # 设置刻度为5的倍数
    labels = function(x) as.character(x),  # 将刻度标签转为字符类型
    minor_breaks = seq(0, ceiling(max(dt$SNP_COUNT) * 1.1), by = 5)  # 添加纵坐标次要网格线
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major.x = element_blank(),  # 移除 x 轴主要网格线
    panel.grid.minor.x = element_blank(),  # 移除 x 轴次要网格线
    panel.grid.major.y = element_line(color = "lightgray", linetype = "dashed"),  # 添加纵坐标主要网格线
    panel.grid.minor.y = element_blank(),  # 移除纵坐标次要网格线
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14, color = "black"),  # 调整字体大小
    axis.text = element_text(size = 12, color = "black", angle = 0),  # 调整字体大小
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),  # 调整图片边距
    legend.position = "none",
    plot.title = element_text(size = 16)  # 调整字体大小
  )

# 保存图片
ggsave(file = output_file, plot = popPlot, width = 12, height = 6, dpi = 600)  # 调整图片大小
