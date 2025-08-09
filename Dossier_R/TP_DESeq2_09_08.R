# -----------------------------------------
# TP Analyse d'expression différentielle avec DESeq2
# -----------------------------------------
# Ce TP guide l'utilisateur de l'import des données RNA-seq bulk 
# à la visualisation et l'interprétation des résultats,
# incluant un volcano plot personnalisé avec ggplot2.
#
# Références principales :
# Love et al., Genome Biol. 2014 (DESeq2)
# Himes et al., PLoS ONE 2014 (jeu airway)
#
# Assurez-vous d'avoir installé les packages nécessaires avant de lancer ce TP.
# -----------------------------------------


# 0. Installation et chargement des packages
# -------------------------------------------
# Décommentez si nécessaire
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "airway", "AnnotationDbi", "org.Hs.eg.db", 
#                        "pheatmap", "RColorBrewer", "EnhancedVolcano", "IHW", "clusterProfiler",
#                        "enrichplot", "DOSE", "pathview", "ReactomePA", "ggplot2", "dplyr", "ggrepel"))

library(DESeq2)
library(airway)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(IHW)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(pathview)
library(ReactomePA)
library(ggplot2)
library(dplyr)
library(ggrepel)

# 1. Chargement et inspection des données
# ---------------------------------------
# Chargement du jeu airway, un jeu de données RNA-seq bulk sur cellules bronchiques humaines
data(airway)
se <- airway

# Afficher un résumé de la structure des données (dimensions, types)
print("Structure et métadonnées des données airway :")
print(paste("Dimensions : ", nrow(se), " gènes x ", ncol(se), " échantillons"))
print(paste("Type d'objet : ", class(se)[1]))
metadata_df <- as.data.frame(colData(se))
print("Exemple de métadonnées d'échantillons :")
print(head(metadata_df))

# Accès aux comptages bruts (raw counts)
raw_counts <- assays(se)[[1]]
print("Premiers comptages bruts (extrait 6x3):")
print(raw_counts[1:6, 1:3])

# 2. Préparation de l'objet DESeqDataSet
# --------------------------------------
# Création d'un objet DESeqDataSet avec design prenant en compte les lignées cellulaires et le traitement dexaméthasone
dds <- DESeqDataSet(se, design = ~ cell + dex)

# Relevel du facteur dex pour que "untrt" soit la référence
dds$dex <- relevel(dds$dex, ref = "untrt")

# Filtrer les gènes faiblement exprimés 
# Critère : au moins 10 lectures dans au moins 4 échantillons
keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep,]
print(paste("Dimensions après filtrage des gènes faiblement exprimés:", dim(dds)[1], "gènes x", dim(dds)[2], "échantillons"))

# 3. Analyse d'expression différentielle
# ---------------------------------------
# Lancer la fonction DESeq pour normalisation, estimation de la dispersion et tests statistiques
dds <- DESeq(dds)

# Extraction des facteurs de normalisation pour vérifier
print("Facteurs de normalisation (size factors) par échantillon:")
print(sizeFactors(dds))

# Visualiser la dispersion estimée (qualité du modèle)
plotDispEsts(dds)

# 4. Extraction et annotation des résultats
# ----------------------------------------
# Extraction des résultats contrastés entre traité (trt) et non traité (untrt)
res <- results(dds, contrast = c("dex", "trt", "untrt"))

# Shrinkage des log2 fold changes pour stabilité statistique
res <- lfcShrink(dds, contrast = c("dex", "trt", "untrt"), res = res, type = "ashr")

# Annotation des gènes (symboles, noms, entrez) pour faciliter l'interprétation
res$symbol <- mapIds(org.Hs.eg.db, keys=row.names(res), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
res$gene_name <- mapIds(org.Hs.eg.db, keys=row.names(res), column="GENENAME", keytype="ENSEMBL", multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db, keys=row.names(res), column="ENTREZID", keytype="ENSEMBL", multiVals="first")

# Trier les résultats par p-adj croissant
resOrdered <- res[order(res$padj),]
print("Top 10 gènes différentiellement exprimés :")
print(head(as.data.frame(resOrdered), 10))

# 5. Visualisations avancées
# --------------------------
# Transformation variance stabilisante pour visualisation (log transformation corrigée)
vsd <- vst(dds, blind = FALSE)

# Analyse en composantes principales pour vérifier la séparation des groupes
pcaData <- plotPCA(vsd, intgroup=c("dex", "cell"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=dex, shape=cell)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  theme_bw()

# Heatmap des 50 gènes les plus variables
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat <- assay(vsd)[topVarGenes, ]
mat <- mat - rowMeans(mat)  # centrage par ligne
annotation_col <- as.data.frame(colData(vsd)[, c("dex", "cell")])
pheatmap(mat,
         annotation_col = annotation_col,
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, main = "Heatmap : 50 gènes les plus variables")

# Visualisation des comptages normalisés pour un gène d’intérêt le plus significatif
topGene <- rownames(resOrdered)[1]
topSymbol <- resOrdered$symbol[1]
normCounts <- counts(dds, normalized=TRUE)
barData <- data.frame(
  Sample = colnames(dds),
  Count = normCounts[topGene,],
  Condition = dds$dex,
  Cell = dds$cell
)
ggplot(barData, aes(x=Sample, y=Count, fill=Condition)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste("Comptages normalisés pour ", topSymbol),
       x = "Échantillon",
       y = "Comptages normalisés")

# 6. Fonction personnalisée pour Volcano Plot en ggplot2
# ------------------------------------------------------
my_ggvolcano_new_blue2 <- function(data, p.thresh, FC, annot.pos, title="") {
  options(ggrepel.max.overlaps = Inf)
  
  data <- data %>%
    na.omit() %>%
    mutate(color = case_when(
      log2FoldChange > FC & padj < p.thresh ~ "Increased",
      log2FoldChange < -FC & padj < p.thresh ~ "Decreased",
      abs(log2FoldChange) < FC ~ "nonsignificant",
      padj > p.thresh ~ "nonsignificant"
    ))
  
  vol <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj), color = color)) + 
    geom_point(size = 2.5, alpha = 0.8, na.rm = TRUE) +
    scale_color_manual(
      name = "status",
      values = c(
        Increased = "gold4",
        Decreased = "blue4",
        nonsignificant = "darkgray"
      )
    ) +
    theme_light(base_size = 12) + 
    theme(legend.position = "right") + 
    xlab("Log2 Fold Change") + 
    ylab(expression(-log[10]("adjusted p-value"))) +
    geom_hline(yintercept = -log10(p.thresh), colour = "darkgrey") + 
    annotate(geom="text", label=paste0("p.adj ", p.thresh), x=min(data$log2FoldChange + annot.pos), y=-log10(p.thresh), vjust=-1) +
    geom_vline(xintercept = FC, linetype="dotted") +
    annotate(geom="text", label=paste0("lFC >", FC), x=1 + FC, y=0, vjust=-1) +
    geom_vline(xintercept = -FC, linetype="dotted") +
    annotate(geom="text", label=paste0("lFC < -", FC), x=-1 - FC, y=0, vjust=-1) +
    scale_y_continuous(trans = "log1p") + 
    ggtitle(title) +
    geom_text_repel(data = data[na.omit(which(data$padj < p.thresh & data$log2FoldChange > FC))[1:5], ], 
                    ggplot2::aes(label = data[na.omit(which(data$padj < p.thresh & data$log2FoldChange > FC))[1:5], "gene_name"])) +
    geom_text_repel(data = data[na.omit(which(data$padj < p.thresh & data$log2FoldChange < (-FC)))[1:5], ], 
                    ggplot2::aes(label = data[na.omit(which(data$padj < p.thresh & data$log2FoldChange < (-FC)))[1:5], "gene_name"]))
  return(vol)
}

# Appel de la fonction pour le volcano plot
volcano_plot <- my_ggvolcano_new_blue2(as.data.frame(res), p.thresh = 0.05, FC = 1, annot.pos = -5, title = "Volcano plot DESeq2 airway")
print(volcano_plot)

# 7. Amélioration statistique avec IHW (Independent Hypothesis Weighting)
# ------------------------------------------------------------------------
resIHW <- results(dds, filterFun = ihw, contrast = c("dex", "trt", "untrt"))
print("Résumé des résultats avec IHW :")
summary(resIHW)

print(paste("Nombre de gènes détectés avec p-adj < 0.1 : méthode standard =", sum(res$padj < 0.1, na.rm=TRUE)))
print(paste("Nombre de gènes détectés avec p-adj < 0.1 : méthode IHW =", sum(resIHW$padj < 0.1, na.rm=TRUE)))

# 8. Export des résultats significatifs
# -------------------------------------
resSig <- subset(res, padj < 0.05)
resSig <- resSig[order(resSig$padj),]
resSigDF <- as.data.frame(resSig)

# write.csv(resSigDF, file="DESeq2_results_significant.csv")  # Décommentez pour exporter

# 9. Analyse fonctionnelle (exemple GO enrichment)
# -----------------------------------------------
sig_entrez <- resSig$entrez[!is.na(resSig$entrez)]

all_genes <- res[!is.na(res$entrez),]
gene_list <- all_genes$log2FoldChange
names(gene_list) <- all_genes$entrez
gene_list <- sort(gene_list, decreasing = TRUE)

ego <- enrichGO(
  gene = sig_entrez,
  universe = all_genes$entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE
)

if (nrow(as.data.frame(ego)) > 0) {
  print("Top termes GO enrichis :")
  print(head(as.data.frame(ego), 10))
  barplot(ego, showCategory = 10) + ggtitle("Enrichissement GO : Processus biologiques")
} else {
  print("Aucun terme GO significativement enrichi trouvé.")
}

# Fin du TP
# ---------
# Ce script guide les étapes essentielles pour faire une analyse d'expression différentielle
# complète avec DESeq2, depuis le chargement des données brutes jusqu’à l’interprétation
# biologique par l’enrichissement fonctionnel, en intégrant une visualisation personnalisée
# avec un volcano plot élégant et annoté.

