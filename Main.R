library("tidyverse")

# 1. Read expression data into a data frame
alz.expr.data <- read.csv("alz_expr_data.csv");
genes <- alz.expr.data$X;

temp <- alz.expr.data[,-1];
rownames(temp) <- alz.expr.data[,1];
alz.expr.data <- as.data.frame(temp);
samples <- colnames(alz.expr.data);

# 2.
# 2.1
norm.expr.data <- apply(X = alz.expr.data[,1:100], MARGIN=2, FUN=scale);
boxplot(norm.expr.data, col="cyan", cex.axis = 0.35, las = 2);

# 2.2
log2.expr.data <- apply(X = alz.expr.data[,1:100], MARGIN=2, FUN=log2);
boxplot(log2.expr.data, col="pink", cex.axis = 0.35, las = 2);

# 2.3
minmax.expr.data <- apply(alz.expr.data[,1:100], MARGIN=2, function(x) {
  (x - min(x)) / (max(x) - min(x))
});
boxplot(minmax.expr.data, col="yellow", cex.axis = 0.35, las = 2)

# 3.
heatmap.data <- head(log2.expr.data, 50);
heatmap(as.matrix(heatmap.data));

# 4.
alz.condition.data <- read.table("alz_condition.csv", header = TRUE);
condition <- as.vector(alz.condition.data[,1]);
log2.expr.cond.data <- add_column(as.data.frame(log2.expr.data), condition=condition, .before=1)

# 5.
log2.expr.control <- log2.expr.cond.data[log2.expr.cond.data$condition == "control",-1];
log2.expr.disease <- log2.expr.cond.data[log2.expr.cond.data$condition == "alzheimer",-1];

control.mean.scores <- as.data.frame(apply(log2.expr.control[,1:300], MARGIN=2, mean));
colnames(control.mean.scores) <- "mean control";

disease.mean.scores <- as.data.frame(apply(log2.expr.disease[,1:300], MARGIN=2, mean));
colnames(disease.mean.scores) <- "mean alzheimer";

mean.scores <- cbind(control.mean.scores, disease.mean.scores);
df.log2FC <- t(as.data.frame(mean.scores$"mean alzheimer" - mean.scores$"mean control" ))
rownames(df.log2FC) <- "Log2FC";
colnames(df.log2FC) <- genes;
log2.expr.data <- rbind(log2.expr.data, df.log2FC);

# 6
pvalues <- c();
for (c in 1:300) {
  ttest <- t.test(log2.expr.data[1:100, c]);
  pvalues <- c(pvalues, ttest$p.value);
}
df.pvalues <- t(as.data.frame(pvalues));
colnames(df.pvalues) <- genes;
log2.expr.data <- rbind(log2.expr.data, df.pvalues);

df.fdr <- t(as.data.frame(p.adjust(log2.expr.data[102,], method="fdr")));
rownames(df.fdr) <- "p-adjusted";
colnames(df.fdr) <- genes;

log2.expr.data <- rbind(log2.expr.data, df.fdr);

t.log2.expr.data <- t(log2.expr.data);

# 7
library("dplyr");
df.deg <- subset(t.log2.expr.data, abs(t.log2.expr.data[,101]) >= 0.3 & t.log2.expr.data[,103] <= 0.8)

#8
df.deg <- df.deg[order(df.deg[,101]),]
heatmap(df.deg[,1:100])
