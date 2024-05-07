# 6.
# library("stats")
# pvalues <- c()
# fdrs <- c()
# for (r in 1:300) {
#   if (var(as.numeric(mean.scores[r, 1:2])) == 0) {
#     pvalues <- c(pvalues, NA);
#     fdrs <- c(fdrs, NA)
#     next;
#   }
#   ttest <- t.test(mean.scores[r, 1:2]);
#   pvalues <- c(pvalues, ttest$p.value);
#   false.discovery.rate <- p.adjust(ttest$p.value, method="fdr");
#   fdrs <- c(fdrs, false.discovery.rate);
# }
# df.pvalues <- as.data.frame(pvalues);
# df.fdrs <- as.data.frame(fdrs);
# mean.scores <- cbind(mean.scores, df.pvalues, df.fdrs);


# iii. 6
#pvalues <- c();
#for (c in 1:100) {
#  if (var(alz.expr.data[,c]) == 0) {
#    pvalues <- c(pvalues, NA)
#    next;
#  }
#  pvalues <- c(pvalues, t.test(alz.expr.data[,c])$p.value);
#}
#df.pvalues <- t(as.data.frame(pvalues))
#colnames(df.pvalues) <- samples
#alz.expr.data <- rbind(alz.expr.data, df.pvalues);

### - iii

# ii.
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
#library("DESeq2")
#condition
#DESeqDataSetFromMatrix()
#

#control.mean.scores <- as.data.frame(apply(log2.expr.control[,1:300], MARGIN=2, mean));
#colnames(control.mean.scores) <- "mean control";

#disease.mean.scores <- as.data.frame(apply(log2.expr.disease[,1:300], MARGIN=2, mean));
#colnames(disease.mean.scores) <- "mean alzheimer";

#
#control.means <- t(as.data.frame(colMeans(log2.expr.control)));
#rownames(control.means) <- "mean control"
#log2.expr.control <- rbind(log2.expr.control, control.means);
#
#disease.means <- t(as.data.frame(colMeans(log2.expr.disease)));
#rownames(disease.means) <- "mean disease"
#log2.expr.disease <- rbind(log2.expr.disease, disease.means);
#
#pvalues.control <- c();
#for (c in 1:300) {
#  if (var(log2.expr.control[-51,c]) == 0) {
#    pvalues.control <- c(pvalues.control, NA)
#    next;
#  }
#  pvalues.control <- c(pvalues.control, t.test(log2.expr.control[-51,c])$p.value);
#}
#
#df.pvalues.control <- t(as.data.frame(pvalues.control));
#rownames(df.pvalues.control) <- "pvalue control";
#colnames(df.pvalues.control) <- genes;
#log2.expr.control <- rbind(log2.expr.control, df.pvalues.control);
#
#pvalues.disease <- c()
#for (c in 1:300) {
#  if (var(log2.expr.disease[-51,c]) == 0) {
#    pvalues.disease <- c(pvalues.disease, NA);
#    next;
#  }
#  pvalues.disease <- c(pvalues.disease, t.test(log2.expr.disease[-51,c])$p.value);
#}
#
#df.pvalues.disease <- t(as.data.frame(pvalues.disease));
#rownames(df.pvalues.disease) <- "pvalue disease";
#colnames(df.pvalues.disease) <- genes;
#log2.expr.disease <- rbind(log2.expr.disease, df.pvalues.disease);

#
# 6.1 -
pvalues <- c();
fdrs <- c();
for (r in 1:300) {
  ttest <- t.test(as.numeric(log2.expr.data[1:100, r]));
  pvalues <- c(pvalues, ttest$p.value);
}
pvalues
p.adjust(pvalues, method="fdr");

t.test(as.numeric(alz.expr.data[1, 1:100]))

pvalues.control <- c();
for (c in 1:300) {
  if (var(log2.expr.control[-51,c]) == 0) {
    pvalues.control <- c(pvalues.control, NA)
    next;
  }
  pvalues.control <- c(pvalues.control, t.test(log2.expr.control[-51,c])$p.value);
}

df.pvalues.control <- as.data.frame(pvalues.control);
rownames(df.pvalues.control) <- genes;

pvalues.disease <- c()
for (c in 1:300) {
  if (var(log2.expr.disease[-51,c]) == 0) {
    pvalues.disease <- c(pvalues.disease, NA);
    next;
  }
  pvalues.disease <- c(pvalues.disease, t.test(log2.expr.disease[-51,c])$p.value);
}

df.pvalues.disease <- as.data.frame(pvalues.disease);
rownames(df.pvalues.disease) <- genes;

mean.scores <- cbind(mean.scores, df.pvalues.control, df.pvalues.disease);

fdrs <- c();
for (c in 1:300) {
  print(p.adjust(as.numeric(mean.scores[c, 4:5]), method="fdr"))
}
