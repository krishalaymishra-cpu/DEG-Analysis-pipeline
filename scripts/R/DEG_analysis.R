#THE DEG ANALYSIS WAS DONE USING THREE LIBRARIES-DESeq2, edgeR and Limma-voom

#================================DESEQ2=======================================================

library(DESeq2)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("apeglm")
}
BiocManager::install("apeglm")
library(apeglm)
all(colnames(count_data)==rownames(col_data))
all(rownames(col_data)%in%colnames(count_data))

dds_analysis<-DESeqDataSetFromMatrix(countData = count_data,
                                     colData = col_data,
                                     design = ~ gender + condition)
dds<-estimateSizeFactors(dds_analysis)
dds<-DESeq(dds)
levels(dds$condition)
result<-results(dds,contrast = c("condition","Benign adjacent","Primary tumor"))
resultsNames(dds)
result <- lfcShrink(dds,
                        coef = "condition_Primary.tumor_vs_Benign.adjacent",
                        type = "apeglm")
sig_res <- subset(result, padj < 0.05 & log2FoldChange > 0.5)


sig_res<-subset(result, padj < 0.05 & log2FoldChange < -0.5) 



# Convert results to a data frame
res_df <- as.data.frame(result)
res_df$Gene <- rownames(res_df)

# Add regulation status
res_df$Regulation <- ifelse(res_df$log2FoldChange > 0 & res_df$padj < 0.05, "Upregulated",
                            ifelse(res_df$log2FoldChange < 0 & res_df$padj < 0.05, "Downregulated",
                                   "Not Significant"))


#=================================limma(voom)=======================================================
# 1. Load libraries
library(edgeR)
library(limma)

# 2. Create DGEList and normalize
dge <- DGEList(counts = count_data)
dge <- calcNormFactors(dge)

# 3. Build design matrix (two factors: gender + condition)
design <- model.matrix(~ gender + condition, data = col_data)
colnames(design) <- make.names(colnames(design))
colnames(design)

# 4. Apply voom transformation
v <- voom(dge, design, plot = TRUE)

# 5. Fit linear model
fit <- lmFit(v, design)
fit <- eBayes(fit)

# 6. Extract results for condition effect (Primary tumor vs Benign adjacent)
res_condition <- topTable(fit, coef="conditionPrimary.tumor", number=Inf)

# 7. Extract results for gender effect (male vs female)
res_gender <- topTable(fit, coef="gendermale", number=Inf)


#Condition
res_down<-subset(res_condition, adj.P.Val < 0.05 & logFC < -0.5)


res_up<-subset(res_condition, adj.P.Val < 0.05 & logFC > 0.5)


#Gender
res_down_gender<-subset(res_gender, adj.P.Val < 0.05 & logFC < -0.5)


res_up<-subset(res_gender, adj.P.Val < 0.05 & logFC > 0.5)


#======================================EdgeR=========================================================
library(edgeR)

# 1. Create DGEList
dge <- DGEList(counts = count_data, samples = col_data)

# 2. Make sure metadata columns are factors
col_data$condition <- factor(col_data$condition,
                             levels = c("Benign adjacent","Primary tumor"))
col_data$gender <- factor(col_data$gender,
                          levels = c("female","male"))

# 3. Build design matrix with both condition + gender (this is your "batch")
design <- model.matrix(~ gender + condition, data = col_data)
colnames(design) <- make.names(colnames(design))
colnames(design)
# "X.Intercept." "gendermale" "conditionPrimary.tumor"

# 4. Normalize and estimate dispersion
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, design)

# 5. Fit GLM
fit <- glmFit(dge, design)

# 6. Test contrasts
lrt_condition <- glmLRT(fit, coef="conditionPrimary.tumor")
lrt_gender    <- glmLRT(fit, coef="gendermale")

# 7. Extract DE genes
res_condition_edge <- topTags(lrt_condition, n=Inf)$table
res_gender_edge    <- topTags(lrt_gender, n=Inf)$table

#Condition
edge_res_down<-subset(res_condition_edge, FDR < 0.05 & logFC < -0.5)


edge_res_up<-subset(res_condition_edge, FDR < 0.05 & logFC > 0.5)
