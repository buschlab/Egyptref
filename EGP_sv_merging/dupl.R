library(ggplot2)
library(IRanges)
library(data.table)

setwd("C:\\Users\\Matthias\\Documents\\Postdoc\\Projects\\2019-02-08_-_Egyptian_genome\\sv")

df = read.table("dupl.txt", header=TRUE, sep="\t", as.is=TRUE, na.strings = c("NA"))

dt = as.data.table(df)

dt[,group := { 
  ir =  IRanges(pos, end);
  subjectHits(findOverlaps(ir, reduce(ir)))
},by=list(sample,chr)]


res = dt[, list(pos=min(pos), end=max(end), core_length=min(end)-max(pos), tsv_id=list(tsv_id), 
                irange=paste0("IRanges(c(", do.call(paste, c(as.list(pos), sep = ",")), "), c(", do.call(paste, c(as.list(end), sep = ",")), "))")),
         by=list(group,sample,chr)]

res$length = res$end-res$pos
res$sv = paste0(res$chr, ":", res$pos, "-", res$end)

res$core_length[which(res$core_length < 0)] = 0


# Filter
# res = res[!which(res$length < 300 | res$chr == "chrX" | res$chr == "chrY"),]


c(round(mean(table(res$sample))), round(sd(table(res$sample))), round(min(table(res$sample))), round(max(table(res$sample))))


fwrite(res, "dupl.collapsed.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)


sv_count = as.data.frame(table(res$sample))
df_plot = data.frame(sample=factor(res$sample, levels = sv_count[order(sv_count$Freq, decreasing=TRUE), 1]), 
                     chr=factor(res$chr, levels = paste0("chr", c(1:22, "X", "Y"))), #
                     length=res$length,
                     core_length=res$core_length,
                     sv=res$sv)


ggplot(df_plot, aes(x=chr, y=core_length/length)) +
  geom_boxplot(outlier.alpha = 0.2) +
  xlab("") + 
  ylab("core_length/length") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))

ggsave("dupl.collapsed.core_length_per_chr.png", plot=last_plot(), dpi=300, width=20, height=10, scale=1.2, units="cm")


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

ggsave("dupl.collapsed.counts_per_ind_chr.png", plot=last_plot(), dpi=300, width=20, height=10, scale=1.2, units="cm")


ggplot(df_plot, aes(x=sample, y=log10(length))) +
  geom_boxplot(outlier.alpha = 0.2) +
  xlab("") + 
  ylab("log10(length)") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(angle = 45, vjust = 1, hjust=1, size=4))

ggsave("dupl.collapsed.lengths_per_ind.png", plot=last_plot(), dpi=300, width=20, height=10, scale=1.2, units="cm")



uniq_df_plot = df_plot[!duplicated(df_plot$sv),]
log_count = data.frame(chr=c(), count=c(), log10_thr=c())

tmp = as.data.frame(table(uniq_df_plot[which(log10(uniq_df_plot$length) < 2), "chr"]))
tmp$log10_thr = "log10(length) < 2"
colnames(tmp) = c("chr", "count", "log10_thr")
log_count = rbind(log_count, tmp)

for (i in 2:6){
  tmp = as.data.frame(table(uniq_df_plot[which(log10(uniq_df_plot$length) >= i & log10(uniq_df_plot$length) < i+1), "chr"]))
  tmp$log10_thr = paste(i, "<= log10(length) <", i+1, sep=" ")
  colnames(tmp) = c("chr", "count", "log10_thr")
  log_count = rbind(log_count, tmp)
}

tmp = as.data.frame(table(uniq_df_plot[which(log10(uniq_df_plot$length) >= 7), "chr"]))
tmp$log10_thr = "7 <= log10(length)"
colnames(tmp) = c("chr", "count", "log10_thr")
log_count = rbind(log_count, tmp)


ggplot(log_count, aes(x=chr, y=count)) +
  geom_bar(stat = "identity", aes(fill=factor(log10_thr, levels=unique(log10_thr)))) +
  scale_y_continuous(expand = c(0,0.05)) +
  xlab("") + 
  ylab("Count") +
  theme_bw() +
  guides(fill=guide_legend(title="Groups")) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))

ggsave("dupl.collapsed.counts_per_chr.png", plot=last_plot(), dpi=300, width=20, height=10, scale=1.2, units="cm")



sample_count = data.frame(sample=c(), count=c(), log10_thr=c())

tmp = as.data.frame(table(df_plot[which(log10(df_plot$length) < 2), "sample"]))
tmp$log10_thr = "log10(length) < 2"
colnames(tmp) = c("sample", "count", "log10_thr")
sample_count = rbind(sample_count, tmp)

for (i in 2:6){
  tmp = as.data.frame(table(df_plot[which(log10(df_plot$length) >= i & log10(df_plot$length) < i+1), "sample"]))
  tmp$log10_thr = paste(i, "<= log10(length) <", i+1, sep=" ")
  colnames(tmp) = c("sample", "count", "log10_thr")
  sample_count = rbind(sample_count, tmp)
}

tmp = as.data.frame(table(df_plot[which(log10(df_plot$length) >= 7), "sample"]))
tmp$log10_thr = "7 <= log10(length)"
colnames(tmp) = c("sample", "count", "log10_thr")
sample_count = rbind(sample_count, tmp)

ggplot(sample_count, aes(x=sample, y=count)) +
  geom_bar(stat = "identity", aes(fill=factor(log10_thr, levels=unique(log10_thr)))) +
  scale_y_continuous(expand = c(0,0.05)) +
  xlab("") + 
  ylab("Count") +
  theme_bw() +
  guides(fill=guide_legend(title="Groups")) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(angle = 45, vjust = 1, hjust=1, size=4))

ggsave("dupl.collapsed.counts_per_ind_log10.png", plot=last_plot(), dpi=300, width=20, height=10, scale=1.2, units="cm")
