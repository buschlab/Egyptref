library(ggplot2)
library(IRanges)
library(data.table)

setwd("C:\\Users\\Matthias\\Documents\\Postdoc\\Projects\\2019-02-08_-_Egyptian_genome\\sv")

df = read.table("transloc.txt", header=TRUE, sep="\t", as.is=TRUE, na.strings = c("NA"))

dt = as.data.table(df)

res = dt[, list(tsv_id=list(tsv_id)),
         by=list(sample, chr, pos, chr2, end)]

res$sv = paste0(res$chr, ":", res$pos, "-", res$chr2, ":", res$end);

c(round(mean(table(res$sample))), round(sd(table(res$sample))), round(min(table(res$sample))), round(max(table(res$sample))))


fwrite(res, "transloc.collapsed.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)


sv_count = as.data.frame(table(res$sample))
df_plot = data.frame(sample=factor(res$sample, levels = sv_count[order(sv_count$Freq, decreasing=TRUE), 1]), 
                     chr=factor(res$chr, levels = paste0("chr", c(1:22, "X", "Y"))), sv=res$sv)


ggplot(df_plot, aes(x=sample, fill=chr)) +
  geom_bar() +
  scale_y_continuous(expand = c(0,0.05)) +
  xlab("") + 
  ylab(" Count") +
  theme_bw() +
  guides(fill=guide_legend(title="Chromosome")) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(angle = 45, vjust = 1, hjust=1, size=4))

ggsave("transloc.collapsed.counts_per_ind_chr.png", plot=last_plot(), dpi=300, width=20, height=10, scale=1.2, units="cm")


ggplot(df_plot[!duplicated(df_plot$sv),], aes(x=chr)) +
  geom_bar() +
  scale_y_continuous(expand = c(0,0.05)) +
  xlab("") + 
  ylab("Count") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))

ggsave("transloc.collapsed.counts_per_chr.png", plot=last_plot(), dpi=300, width=20, height=10, scale=1.2, units="cm")

