rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(1234)

###############################################################################
## 0. Working directory and output folders
###############################################################################

setwd("Your path")

out_dirs <- c(
    "01_metadata",
    "02_expression",
    "03_PCA",
    "04_DEG",
    "05_volcano",
    "06_heatmap",
    "07_selected_genes",
    "07_selected_genes/individual",
    "08_GO_KEGG",
    "09_GSEA"
)

invisible(lapply(out_dirs, dir.create, showWarnings = FALSE, recursive = TRUE))

###############################################################################
## 1. Install/load packages
###############################################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

cran_pkgs <- c(
    "ggplot2",
    "ggrepel",
    "pheatmap",
    "dplyr",
    "stringr",
    "msigdbr",
    "ggridges"
)

bioc_pkgs <- c(
    "GEOquery",
    "limma",
    "clusterProfiler",
    "org.Hs.eg.db",
    "enrichplot"
)

for (pkg in cran_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}

for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
}

suppressPackageStartupMessages({
    library(GEOquery)
    library(limma)
    library(ggplot2)
    library(ggrepel)
    library(pheatmap)
    library(dplyr)
    library(stringr)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(msigdbr)
    library(ggridges)
})

###############################################################################
## 2. Download/read GSE159104 supplementary files
###############################################################################

dge_file <- file.path("GSE159104", "GSE159104_DGE_raw_154samples.txt.gz")
sample_file <- file.path("GSE159104", "GSE159104_samples.txt.gz")

if (!file.exists(dge_file) || !file.exists(sample_file)) {
    getGEOSuppFiles(
        "GSE159104",
        makeDirectory = TRUE,
        baseDir = "."
    )
}

stopifnot(file.exists(dge_file))
stopifnot(file.exists(sample_file))

cat("Reading GSE159104 expression matrix...\n")

dge <- read.delim(
    gzfile(dge_file),
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
)

sample_info <- read.delim(
    gzfile(sample_file),
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
)

cat("DGE dimension:", paste(dim(dge), collapse = " x "), "\n")
cat("Sample-info dimension:", paste(dim(sample_info), collapse = " x "), "\n")

required_sample_cols <- c("Sample Name", "ADHD_Status", "Channel")
if (!all(required_sample_cols %in% colnames(sample_info))) {
    stop(
        "The sample annotation file does not contain the expected columns: ",
        paste(required_sample_cols, collapse = ", ")
    )
}

###############################################################################
## 3. Match expression columns to official sample mapping
###############################################################################

matched <- sample_info$Channel %in% colnames(dge)
cat("Matched channels:", sum(matched), "/", nrow(sample_info), "\n")

sample_info <- sample_info[matched, , drop = FALSE]

if (nrow(sample_info) == 0) {
    stop("No sample channels matched the DGE matrix.")
}

expr_raw <- as.matrix(
    dge[, sample_info$Channel, drop = FALSE]
)
storage.mode(expr_raw) <- "numeric"

rownames(expr_raw) <- make.unique(as.character(dge[[1]]))
colnames(expr_raw) <- sample_info$`Sample Name`

cat("Raw expression dimension:", paste(dim(expr_raw), collapse = " x "), "\n")
cat("Raw NA count:", sum(is.na(expr_raw)), "\n")

###############################################################################
## 4. Library-size normalization to CPM
###############################################################################

lib_size <- colSums(expr_raw, na.rm = TRUE)
if (any(!is.finite(lib_size)) || any(lib_size <= 0)) {
    stop("At least one sample has an invalid library size.")
}

expr_cpm <- sweep(expr_raw, 2, lib_size, "/") * 1e6

###############################################################################
## 5. Transcript -> gene-level expression
##    Mean CPM across rows sharing a gene symbol.
###############################################################################

gene_symbol <- as.character(dge[[2]])

valid_gene <- (
    !is.na(gene_symbol) &
        gene_symbol != "" &
        gene_symbol != "NA"
)

expr_cpm_gene_rows <- expr_cpm[valid_gene, , drop = FALSE]
gene_symbol_valid <- gene_symbol[valid_gene]

expr_gene_sum <- rowsum(
    expr_cpm_gene_rows,
    group = gene_symbol_valid,
    reorder = FALSE
)

gene_n <- table(gene_symbol_valid)
expr_gene <- expr_gene_sum /
    as.numeric(gene_n[rownames(expr_gene_sum)])
expr_gene <- as.matrix(expr_gene)

cat("Gene-level matrix:", paste(dim(expr_gene), collapse = " x "), "\n")
cat("Gene-level NA count:", sum(is.na(expr_gene)), "\n")

###############################################################################
## 6. Build subject metadata and merge repeated measurements
###############################################################################

meta <- data.frame(
    sample_name = sample_info$`Sample Name`,
    group = toupper(as.character(sample_info$ADHD_Status)),
    stringsAsFactors = FALSE
)

meta$sample_base <- sub(
    "_rep[0-9]+$",
    "",
    meta$sample_name,
    ignore.case = TRUE
)

meta$subject_id <- sub(
    "_(ADHD|CTRL)$",
    "",
    meta$sample_base,
    ignore.case = TRUE
)

subject_order <- unique(meta$subject_id)

expr_list <- lapply(
    subject_order,
    function(id) {
        idx <- which(meta$subject_id == id)
        if (length(idx) == 1) {
            expr_gene[, idx]
        } else {
            rowMeans(expr_gene[, idx, drop = FALSE], na.rm = TRUE)
        }
    }
)

expr_subject <- do.call(cbind, expr_list)
expr_subject <- as.matrix(expr_subject)
rownames(expr_subject) <- rownames(expr_gene)
colnames(expr_subject) <- subject_order

subject_group <- vapply(
    subject_order,
    function(id) {
        g <- unique(meta$group[meta$subject_id == id])
        g <- g[!is.na(g) & g != ""]
        if (length(g) == 0) NA_character_ else g[1]
    },
    FUN.VALUE = character(1)
)

sample_meta <- data.frame(
    subject_id = subject_order,
    group = subject_group,
    stringsAsFactors = FALSE
)
rownames(sample_meta) <- sample_meta$subject_id

if (any(is.na(sample_meta$group))) {
    stop("Some subjects have missing ADHD/CTRL group labels.")
}

stopifnot(identical(colnames(expr_subject), rownames(sample_meta)))

cat("\nSubject-level groups:\n")
print(table(sample_meta$group))
cat("Repeated-measure distribution:\n")
print(table(table(meta$subject_id)))

write.csv(
    sample_meta,
    "01_metadata/final_sample_metadata.csv",
    row.names = FALSE
)

write.csv(
    data.frame(
        Gene = rownames(expr_subject),
        expr_subject,
        check.names = FALSE
    ),
    "02_expression/GSE159104_gene_CPM_subject_level.csv",
    row.names = FALSE
)

###############################################################################
## 7. Log2(CPM + 1), remove invariant genes, filter low expression
###############################################################################

expr_log <- log2(expr_subject + 1)

keep_gene <- apply(
    expr_log,
    1,
    function(x) {
        sum(is.finite(x)) >= 2 && sd(x, na.rm = TRUE) > 0
    }
)
expr_log <- expr_log[keep_gene, , drop = FALSE]

min_samples <- max(5, ceiling(ncol(expr_log) * 0.20))
keep_expression <- rowSums(expr_log > 1, na.rm = TRUE) >= min_samples
expr_filt <- expr_log[keep_expression, , drop = FALSE]

expr_filt <- avereps(expr_filt, ID = rownames(expr_filt))

cat("Genes retained for analysis:", nrow(expr_filt), "\n")

write.csv(
    data.frame(
        Feature = rownames(expr_filt),
        expr_filt,
        check.names = FALSE
    ),
    "02_expression/GSE159104_final_log2_expression.csv",
    row.names = FALSE
)

###############################################################################
## 8. PCA 
###############################################################################

gene_var <- apply(expr_filt, 1, var, na.rm = TRUE)
gene_var <- sort(gene_var, decreasing = TRUE)
pca_genes <- names(gene_var)[seq_len(min(1000, length(gene_var)))]

pca <- prcomp(
    t(expr_filt[pca_genes, , drop = FALSE]),
    center = TRUE,
    scale. = FALSE
)

pca_var <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

pca_df <- data.frame(
    Sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Group = sample_meta[rownames(pca$x), "group"],
    stringsAsFactors = FALSE
)

pca_df$Group <- factor(pca_df$Group, levels = c("CTRL", "ADHD"))

PCA_COLORS <- c(
    "CTRL" = "#3C8DBC",
    "ADHD" = "#E76F51"
)

p_pca <- ggplot(
    pca_df,
    aes(x = PC1, y = PC2)
) +
    stat_ellipse(
        aes(fill = Group, color = Group),
        geom = "polygon",
        type = "norm",
        level = 0.95,
        alpha = 0.07,
        linewidth = 0.6,
        show.legend = FALSE
    ) +
    stat_ellipse(
        aes(color = Group),
        type = "norm",
        level = 0.95,
        linewidth = 1.0,
        show.legend = FALSE
    ) +
    geom_point(
        aes(fill = Group),
        shape = 21,
        size = 3.8,
        alpha = 0.90,
        color = "grey15",
        stroke = 0.50
    ) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey90") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey90") +
    scale_fill_manual(values = PCA_COLORS) +
    scale_color_manual(values = PCA_COLORS) +
    labs(
        title = "GSE159104",
        x = paste0("PC1 (", sprintf("%.1f", pca_var[1]), "%)"),
        y = paste0("PC2 (", sprintf("%.1f", pca_var[2]), "%)"),
        fill = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        axis.line = element_line(linewidth = 0.8, color = "black"),
        legend.position = "right",
        legend.text = element_text(size = 12)
    )

print(p_pca)

ggsave(
    "03_PCA/PCA_GSE159104_publication.pdf",
    p_pca,
    device = cairo_pdf,
    width = 7.2,
    height = 6.0
)

ggsave(
    "03_PCA/PCA_GSE159104_publication.png",
    p_pca,
    width = 7.2,
    height = 6.0,
    dpi = 600
)

###############################################################################
## 9. limma differential-expression analysis
###############################################################################

stopifnot(identical(colnames(expr_filt), rownames(sample_meta)))

group <- factor(
    sample_meta$group,
    levels = c("CTRL", "ADHD")
)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
rownames(design) <- rownames(sample_meta)

fit <- lmFit(expr_filt, design)
contrast_matrix <- makeContrasts(
    ADHD_vs_CTRL = ADHD - CTRL,
    levels = design
)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

deg <- topTable(
    fit2,
    coef = "ADHD_vs_CTRL",
    number = Inf,
    sort.by = "P",
    adjust.method = "BH"
)

deg$Feature <- rownames(deg)
deg <- deg[, c(
    "Feature",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B"
)]

###############################################################################
## 10. DEG selection: BH-adjusted FDR < 0.05 and |logFC| > log2(1.25)
###############################################################################

FDR_CUTOFF <- 0.05
FC_CUTOFF <- 1.25
LOGFC_CUTOFF <- log2(FC_CUTOFF)

sig_idx <- (
    !is.na(deg$adj.P.Val) &
        deg$adj.P.Val < FDR_CUTOFF &
        abs(deg$logFC) > LOGFC_CUTOFF
)

cat("\n========================================\n")
cat("DEG rule: BH-adjusted FDR < 0.05 & FC > 1.25 / < 0.80\n")
cat("DEGs:", sum(sig_idx, na.rm = TRUE), "\n")
cat("========================================\n")

deg$Change <- "NS"
deg$Change[sig_idx & deg$logFC > LOGFC_CUTOFF] <- "Up"
deg$Change[sig_idx & deg$logFC < -LOGFC_CUTOFF] <- "Down"
deg$Change <- factor(deg$Change, levels = c("Down", "NS", "Up"))

deg_sig <- deg[sig_idx, , drop = FALSE]
deg_sig <- deg_sig[order(deg_sig$adj.P.Val), , drop = FALSE]

deg_up <- deg_sig[deg_sig$logFC > LOGFC_CUTOFF, , drop = FALSE]
deg_down <- deg_sig[deg_sig$logFC < -LOGFC_CUTOFF, , drop = FALSE]

write.csv(deg, "04_DEG/DEG_all_limma.csv", row.names = FALSE)
write.csv(deg_sig, "04_DEG/DEG_FDR0.05_FC1.25.csv", row.names = FALSE)
write.csv(deg_up, "04_DEG/DEG_FDR0.05_FC1.25_UP.csv", row.names = FALSE)
write.csv(deg_down, "04_DEG/DEG_FDR0.05_FC1.25_DOWN.csv", row.names = FALSE)

write.csv(
    data.frame(
        DEG_rule = "BH-adjusted FDR < 0.05 & |logFC| > log2(1.25)",
        FC_cutoff = FC_CUTOFF,
        log2FC_cutoff = LOGFC_CUTOFF,
        N_total = nrow(deg_sig),
        N_up = nrow(deg_up),
        N_down = nrow(deg_down)
    ),
    "04_DEG/DEG_selection_summary.csv",
    row.names = FALSE
)

###############################################################################
## 11. Volcano plot (FDR-based)
###############################################################################

deg$plot_sig_value <- deg$adj.P.Val
deg$minus_log10_sig <- -log10(pmax(deg$plot_sig_value, 1e-300))

n_up <- sum(deg$Change == "Up", na.rm = TRUE)
n_down <- sum(deg$Change == "Down", na.rm = TRUE)
n_ns <- sum(deg$Change == "NS", na.rm = TRUE)

down_label <- paste0("Down (", n_down, ")")
ns_label <- paste0("NS (", n_ns, ")")
up_label <- paste0("Up (", n_up, ")")

deg$Change_label <- factor(
    as.character(deg$Change),
    levels = c("Down", "NS", "Up"),
    labels = c(down_label, ns_label, up_label)
)

VOLCANO_COLORS <- c("#3C78D8", "#C8C8C8", "#D94841")
names(VOLCANO_COLORS) <- c(down_label, ns_label, up_label)

TOP_LABEL_N <- 8
label_up <- deg[deg$Change == "Up", , drop = FALSE]
label_down <- deg[deg$Change == "Down", , drop = FALSE]
label_up <- head(label_up[order(label_up$plot_sig_value), , drop = FALSE], TOP_LABEL_N)
label_down <- head(label_down[order(label_down$plot_sig_value), , drop = FALSE], TOP_LABEL_N)
label_deg <- rbind(label_up, label_down)

p_volcano <- ggplot(
    deg,
    aes(x = logFC, y = minus_log10_sig)
) +
    geom_point(
        aes(fill = Change_label),
        shape = 21,
        size = 2.1,
        alpha = 0.80,
        color = "white",
        stroke = 0.12
    ) +
    geom_vline(
        xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF),
        linetype = "dashed",
        linewidth = 0.55,
        color = "grey35"
    ) +
    geom_hline(
        yintercept = -log10(FDR_CUTOFF),
        linetype = "dashed",
        linewidth = 0.55,
        color = "grey35"
    ) +
    ggrepel::geom_text_repel(
        data = label_deg,
        aes(x = logFC, y = minus_log10_sig, label = Feature),
        inherit.aes = FALSE,
        size = 3.1,
        fontface = "italic",
        box.padding = 0.45,
        point.padding = 0.25,
        min.segment.length = 0,
        segment.color = "grey45",
        segment.size = 0.30,
        max.overlaps = Inf,
        seed = 1234
    ) +
    scale_fill_manual(values = VOLCANO_COLORS, drop = FALSE) +
    labs(
        title = "GSE159104: ADHD vs CTRL",
        subtitle = paste0(
            "BH-adjusted FDR < ",
            FDR_CUTOFF,
            " | FC > ",
            FC_CUTOFF,
            " or FC < ",
            sprintf("%.2f", 1 / FC_CUTOFF)
        ),
        x = expression(log[2]~"Fold Change"),
        y = expression(-log[10]~"FDR"),
        fill = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10.5, hjust = 0.5, color = "grey30"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        axis.line = element_line(linewidth = 0.8),
        legend.position = "top",
        legend.text = element_text(size = 10.5)
    )

print(p_volcano)

ggsave(
    "05_volcano/Volcano_GSE159104_publication.pdf",
    p_volcano,
    device = cairo_pdf,
    width = 7.6,
    height = 6.6
)

ggsave(
    "05_volcano/Volcano_GSE159104_publication.png",
    p_volcano,
    width = 7.6,
    height = 6.6,
    dpi = 600
)

###############################################################################
## 12. Top DEG heatmap
###############################################################################

if (nrow(deg_sig) > 0) {
    heat_genes <- head(deg_sig$Feature, min(50, nrow(deg_sig)))
} else {
    warning("No FDR-significant DEGs; heatmap uses top 50 genes by FDR.")
    heat_genes <- head(deg$Feature[order(deg$adj.P.Val)], min(50, nrow(deg)))
}

heat_genes <- intersect(heat_genes, rownames(expr_filt))
heat_expr <- expr_filt[heat_genes, , drop = FALSE]

annotation_col <- data.frame(
    Group = factor(
        sample_meta[colnames(heat_expr), "group"],
        levels = c("CTRL", "ADHD")
    )
)
rownames(annotation_col) <- colnames(heat_expr)

ann_colors <- list(
    Group = c(
        "CTRL" = "#3C8DBC",
        "ADHD" = "#E76F51"
    )
)

heatmap_title <- if (nrow(deg_sig) > 0) {
    "Top DEGs (BH-adjusted FDR < 0.05)"
} else {
    "Top 50 genes by FDR (no DEG passed FDR < 0.05)"
}

pdf("06_heatmap/Heatmap_Top50_DEG.pdf", width = 10, height = 9)
pheatmap(
    heat_expr,
    scale = "row",
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 7,
    border_color = NA,
    clustering_method = "complete",
    main = heatmap_title
)
dev.off()

png("06_heatmap/Heatmap_Top50_DEG.png", width = 10, height = 9, units = "in", res = 600)
pheatmap(
    heat_expr,
    scale = "row",
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 7,
    border_color = NA,
    clustering_method = "complete",
    main = heatmap_title
)
dev.off()

###############################################################################
## 13. Selected genes: one boxplot PDF per gene  
###############################################################################

target_genes <- c(
    "GAB1",
    "PIK3CA",
    "PIK3R1",
    "AKT1",
    "NFE2L2",
    "TRIB3",
    "GPX3",
    "APPL2",
    "GCLC",
    "AKT2",
    "NQO1"
)

gene_check <- data.frame(
    Gene = target_genes,
    In_expression_matrix = target_genes %in% rownames(expr_filt),
    In_limma_results = target_genes %in% deg$Feature
)
write.csv(gene_check, "07_selected_genes/selected_gene_check.csv", row.names = FALSE)
print(gene_check)

genes_present <- target_genes[
    target_genes %in% rownames(expr_filt) &
        target_genes %in% deg$Feature
]

if (length(genes_present) == 0) {
    stop("None of the selected genes were found in both expr_filt and limma results.")
}

gene_stats <- deg[
    match(genes_present, deg$Feature),
    c("Feature", "logFC", "P.Value", "adj.P.Val"),
    drop = FALSE
]
colnames(gene_stats)[1] <- "Gene"
gene_stats$FC <- 2^gene_stats$logFC
gene_stats$Direction <- ifelse(gene_stats$logFC > 0, "Up", "Down")
gene_stats$PlotStat <- gene_stats$adj.P.Val
gene_stats$Significant_FDR <- gene_stats$adj.P.Val < FDR_CUTOFF

format_stat_value <- function(x) {
    if (is.na(x)) return("NA")
    if (x < 0.001) {
        format(x, scientific = TRUE, digits = 2)
    } else {
        sprintf("%.3f", x)
    }
}

gene_stats$StatLabel <- vapply(
    gene_stats$PlotStat,
    function(x) paste0("FDR = ", format_stat_value(x)),
    FUN.VALUE = character(1)
)

write.csv(
    gene_stats,
    "07_selected_genes/selected_gene_limma_statistics.csv",
    row.names = FALSE
)

plot_data <- do.call(
    rbind,
    lapply(
        genes_present,
        function(gene) {
            data.frame(
                Subject = colnames(expr_filt),
                Group = sample_meta[colnames(expr_filt), "group"],
                Gene = gene,
                Expression = as.numeric(expr_filt[gene, colnames(expr_filt)]),
                stringsAsFactors = FALSE
            )
        }
    )
)

plot_data$Group <- factor(plot_data$Group, levels = c("CTRL", "ADHD"))
plot_data <- plot_data[
    is.finite(plot_data$Expression) & !is.na(plot_data$Group),
    ,
    drop = FALSE
]

write.csv(
    plot_data,
    "07_selected_genes/selected_gene_expression_long.csv",
    row.names = FALSE
)

GENE_COLORS <- c(
    "CTRL" = "#3C8DBC",
    "ADHD" = "#E76F51"
)

for (gene in genes_present) {
    df_g <- plot_data[plot_data$Gene == gene, , drop = FALSE]
    st_g <- gene_stats[gene_stats$Gene == gene, , drop = FALSE]

    y_min <- min(df_g$Expression, na.rm = TRUE)
    y_max <- max(df_g$Expression, na.rm = TRUE)
    y_rng <- y_max - y_min
    if (!is.finite(y_rng) || y_rng <= 0) y_rng <- 0.5

    y_bracket <- y_max + 0.08 * y_rng
    y_text <- y_max + 0.18 * y_rng

    p_gene <- ggplot(
        df_g,
        aes(x = Group, y = Expression, fill = Group)
    ) +
        geom_boxplot(
            width = 0.56,
            outlier.shape = NA,
            alpha = 0.72,
            linewidth = 0.65
        ) +
        geom_jitter(
            aes(color = Group),
            width = 0.12,
            size = 2.0,
            alpha = 0.62,
            show.legend = FALSE
        ) +
        geom_segment(
            aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket),
            inherit.aes = FALSE,
            linewidth = 0.6
        ) +
        geom_segment(
            aes(x = 1, xend = 1, y = y_bracket, yend = y_bracket - 0.035 * y_rng),
            inherit.aes = FALSE,
            linewidth = 0.6
        ) +
        geom_segment(
            aes(x = 2, xend = 2, y = y_bracket, yend = y_bracket - 0.035 * y_rng),
            inherit.aes = FALSE,
            linewidth = 0.6
        ) +
        annotate(
            "text",
            x = 1.5,
            y = y_text,
            label = st_g$StatLabel[1],
            size = 4.1
        ) +
        scale_fill_manual(values = GENE_COLORS) +
        scale_color_manual(values = GENE_COLORS) +
        scale_y_continuous(
            expand = expansion(mult = c(0.05, 0.18))
        ) +
        labs(
            title = gene,
            x = NULL,
            y = expression(log[2]~"(CPM + 1)")
        ) +
        theme_classic(base_size = 13) +
        theme(
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            axis.title.y = element_text(size = 13, face = "bold"),
            axis.text.x = element_text(size = 11.5, face = "bold"),
            axis.text.y = element_text(size = 10.5, color = "black"),
            axis.line = element_line(linewidth = 0.7),
            legend.position = "none",
            plot.margin = margin(12, 12, 12, 12)
        )

    print(p_gene)

    safe_gene <- gsub("[^A-Za-z0-9_.-]", "_", gene)

    ggsave(
        file.path("07_selected_genes/individual", paste0(safe_gene, "_boxplot.pdf")),
        p_gene,
        device = cairo_pdf,
        width = 4.5,
        height = 5.0
    )

    ggsave(
        file.path("07_selected_genes/individual", paste0(safe_gene, "_boxplot.png")),
        p_gene,
        width = 4.5,
        height = 5.0,
        dpi = 600
    )
}

expression_summary <- plot_data %>%
    group_by(Gene, Group) %>%
    summarise(
        N = sum(!is.na(Expression)),
        Mean = mean(Expression, na.rm = TRUE),
        SD = sd(Expression, na.rm = TRUE),
        Median = median(Expression, na.rm = TRUE),
        Q1 = quantile(Expression, 0.25, na.rm = TRUE),
        Q3 = quantile(Expression, 0.75, na.rm = TRUE),
        .groups = "drop"
    )

write.csv(
    expression_summary,
    "07_selected_genes/selected_gene_expression_summary.csv",
    row.names = FALSE
)

###############################################################################
## 14. GO / KEGG ORA with BH-adjusted FDR < 0.05
###############################################################################

deg_symbols <- unique(as.character(deg_sig$Feature))
deg_symbols <- deg_symbols[!is.na(deg_symbols) & deg_symbols != ""]
background_symbols <- unique(rownames(expr_filt))
background_symbols <- background_symbols[
    !is.na(background_symbols) & background_symbols != ""
]

if (length(deg_symbols) == 0) {
    warning("No FDR-significant DEGs; GO/KEGG ORA is skipped.")
} else {
    deg_map <- bitr(
        deg_symbols,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
    )

    background_map <- bitr(
        background_symbols,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
    )

    deg_entrez <- unique(deg_map$ENTREZID)
    background_entrez <- unique(background_map$ENTREZID)

    write.csv(deg_map, "08_GO_KEGG/DEG_SYMBOL_ENTREZ_mapping.csv", row.names = FALSE)
    write.csv(background_map, "08_GO_KEGG/Background_SYMBOL_ENTREZ_mapping.csv", row.names = FALSE)

    ego_BP <- enrichGO(
        gene = deg_entrez,
        universe = background_entrez,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        readable = TRUE
    )

    ego_CC <- enrichGO(
        gene = deg_entrez,
        universe = background_entrez,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "CC",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        readable = TRUE
    )

    ego_MF <- enrichGO(
        gene = deg_entrez,
        universe = background_entrez,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "MF",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        readable = TRUE
    )

    ekegg <- enrichKEGG(
        gene = deg_entrez,
        universe = background_entrez,
        organism = "hsa",
        keyType = "ncbi-geneid",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1
    )

    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
        ekegg <- setReadable(
            ekegg,
            OrgDb = org.Hs.eg.db,
            keyType = "ENTREZID"
        )
    }

    go_bp_df <- as.data.frame(ego_BP)
    go_cc_df <- as.data.frame(ego_CC)
    go_mf_df <- as.data.frame(ego_MF)
    kegg_df <- as.data.frame(ekegg)

    write.csv(go_bp_df, "08_GO_KEGG/GO_BP_all.csv", row.names = FALSE)
    write.csv(go_cc_df, "08_GO_KEGG/GO_CC_all.csv", row.names = FALSE)
    write.csv(go_mf_df, "08_GO_KEGG/GO_MF_all.csv", row.names = FALSE)
    write.csv(kegg_df, "08_GO_KEGG/KEGG_all.csv", row.names = FALSE)
}

###############################################################################
## 15. Whole-transcriptome GSEA: Hallmark + Reactome
###############################################################################

gene_list <- deg$t
names(gene_list) <- deg$Feature

keep_gsea <- (
    is.finite(gene_list) &
        !is.na(names(gene_list)) &
        names(gene_list) != ""
)
gene_list <- gene_list[keep_gsea]

# One row per gene: keep the largest |t| when a symbol appears twice.
gsea_rank_df <- data.frame(
    Gene = names(gene_list),
    t = as.numeric(gene_list),
    stringsAsFactors = FALSE
) %>%
    group_by(Gene) %>%
    slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
    ungroup()

gene_list <- gsea_rank_df$t
names(gene_list) <- gsea_rank_df$Gene
gene_list <- sort(gene_list, decreasing = TRUE)

write.csv(
    data.frame(Gene = names(gene_list), t = gene_list),
    "09_GSEA/GSEA_ranked_gene_list_limma_t.csv",
    row.names = FALSE
)

hallmark_msig <- msigdbr(
    db_species = "HS",
    species = "Homo sapiens",
    collection = "H"
)

hallmark_t2g <- hallmark_msig %>%
    dplyr::select(gs_name, gene_symbol) %>%
    distinct()

reactome_msig <- msigdbr(
    db_species = "HS",
    species = "Homo sapiens",
    collection = "C2",
    subcollection = "CP:REACTOME"
)

reactome_t2g <- reactome_msig %>%
    dplyr::select(gs_name, gene_symbol) %>%
    distinct()

set.seed(1234)
gsea_hallmark <- GSEA(
    geneList = gene_list,
    TERM2GENE = hallmark_t2g,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    eps = 0,
    verbose = FALSE,
    seed = TRUE
)

set.seed(1234)
gsea_reactome <- GSEA(
    geneList = gene_list,
    TERM2GENE = reactome_t2g,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    eps = 0,
    verbose = FALSE,
    seed = TRUE
)

hallmark_res <- as.data.frame(gsea_hallmark)
reactome_res <- as.data.frame(gsea_reactome)

hallmark_res$Direction <- ifelse(
    hallmark_res$NES > 0,
    "ADHD enriched",
    "CTRL enriched / ADHD suppressed"
)

reactome_res$Direction <- ifelse(
    reactome_res$NES > 0,
    "ADHD enriched",
    "CTRL enriched / ADHD suppressed"
)

hallmark_res <- hallmark_res[order(hallmark_res$p.adjust, hallmark_res$pvalue), , drop = FALSE]
reactome_res <- reactome_res[order(reactome_res$p.adjust, reactome_res$pvalue), , drop = FALSE]

write.csv(hallmark_res, "09_GSEA/Hallmark_GSEA_all.csv", row.names = FALSE)
write.csv(reactome_res, "09_GSEA/Reactome_GSEA_all.csv", row.names = FALSE)

###############################################################################
## 16. Single-pathway GSEA plots
###############################################################################

format_gsea_p <- function(x) {
    if (is.na(x)) return("NA")
    if (x < 0.001) {
        format(x, scientific = TRUE, digits = 2)
    } else {
        sprintf("%.3f", x)
    }
}

plot_single_gsea <- function(
    gsea_object,
    result_df,
    pathway_id,
    display_title,
    output_prefix
) {
    pathway_res <- result_df[result_df$ID == pathway_id, , drop = FALSE]

    if (nrow(pathway_res) == 0) {
        message("Pathway not found: ", pathway_id)
        return(NULL)
    }

    NES_value <- pathway_res$NES[1]
    P_value <- pathway_res$pvalue[1]
    FDR_value <- pathway_res$p.adjust[1]

    line_color <- ifelse(NES_value > 0, "#D64F4B", "#4472C4")
    direction_text <- ifelse(
        NES_value > 0,
        "Enriched in ADHD",
        "Enriched in CTRL / relatively suppressed in ADHD"
    )

    full_title <- paste0(
        display_title,
        "\nNES = ", sprintf("%.2f", NES_value),
        "    |    P = ", format_gsea_p(P_value),
        "    |    FDR = ", format_gsea_p(FDR_value),
        "\n", direction_text
    )

    p <- gseaplot2(
        gsea_object,
        geneSetID = pathway_id,
        title = full_title,
        base_size = 13,
        color = line_color,
        pvalue_table = FALSE,
        subplots = 1:3,
        rel_heights = c(1.6, 0.45, 1),
        ES_geom = "line"
    )

    print(p)

    grDevices::cairo_pdf(
        filename = file.path("09_GSEA/single_pathway", paste0(output_prefix, ".pdf")),
        width = 8.8,
        height = 6.8
    )
    print(p)
    dev.off()

    png(
        filename = file.path("09_GSEA/single_pathway", paste0(output_prefix, ".png")),
        width = 8.8,
        height = 6.8,
        units = "in",
        res = 600
    )
    print(p)
    dev.off()

    invisible(p)
}

plot_single_gsea(
    gsea_hallmark,
    hallmark_res,
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
    "HALLMARK PI3K/AKT/mTOR Signaling",
    "GSEA_HALLMARK_PI3K_AKT_MTOR_SIGNALING"
)

plot_single_gsea(
    gsea_hallmark,
    hallmark_res,
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK mTORC1 Signaling",
    "GSEA_HALLMARK_MTORC1_SIGNALING"
)

plot_single_gsea(
    gsea_hallmark,
    hallmark_res,
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "HALLMARK Reactive Oxygen Species Pathway",
    "GSEA_HALLMARK_ROS_PATHWAY"
)

plot_single_gsea(
    gsea_hallmark,
    hallmark_res,
    "HALLMARK_PEROXISOME",
    "HALLMARK Peroxisome",
    "GSEA_HALLMARK_PEROXISOME"
)

plot_single_gsea(
    gsea_reactome,
    reactome_res,
    "REACTOME_KEAP1_NFE2L2_PATHWAY",
    "REACTOME KEAP1-NFE2L2 Pathway",
    "GSEA_REACTOME_KEAP1_NFE2L2"
)

plot_single_gsea(
    gsea_reactome,
    reactome_res,
    "REACTOME_NUCLEAR_EVENTS_MEDIATED_BY_NFE2L2",
    "Nuclear Events Mediated by NFE2L2",
    "GSEA_REACTOME_NFE2L2_NUCLEAR_EVENTS"
)

plot_single_gsea(
    gsea_reactome,
    reactome_res,
    "REACTOME_PI3K_AKT_ACTIVATION",
    "REACTOME PI3K/AKT Activation",
    "GSEA_REACTOME_PI3K_AKT_ACTIVATION"
)


)
