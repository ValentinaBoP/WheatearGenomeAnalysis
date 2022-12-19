

piechart_data = read.table("pieChart_data", sep = "\t", stringsAsFactors=FALSE, header = FALSE)
piechart_data$V1 = factor(x = piechart_data$V1, levels = c("DNA", "LINE", "LTR", "SINE", "Unknown", "Host genome"))

colors = c("#e41a1c","#377eb8","#4daf4a","#ff7f00","#984ea3", "grey80")
pie_plot = ggplot(data = piechart_data, aes(x = "", y = V2, fill = V1)) + geom_bar(stat = "identity", width = 1, color="white") + coord_polar("y", start = 0, direction = -1)  + theme_void() +  scale_fill_manual(values = colors)

ggsave(plot = pie_plot, filename = "pieChart_1.0.pdf", units = "cm", device = "pdf", width = 15, height = 15, dpi = 300)