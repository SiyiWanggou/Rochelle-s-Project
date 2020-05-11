#for Rochelle, Mki67, Sox2 (or Atoh1)
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Cairo)
library(tidyverse)
library(pheatmap)
library(patchwork)
library(grid)

#Alter "Sox2" to "Atoh1", we could generate the results for Mki67-Atoh1 combination.
#--------------P21 Math1-cre;SmoM2 Model----------------
setwd("E:/share18/share1/Xin/analysis/data_pre_monocle")
load("MB.Robj")
setwd("E:/share18/share1/Rochelle/P21_Math1_Cre_SmoM2")

#Violin-Plot "Mki67" and "Sox2"
Cairo(file="Violinplot-Mki67 and Sox2 in P21 math1-cre;smoM2 version 1.png",type="png",units="in",bg="white",width=8,height=4,pointsize=114,dpi=300)
VlnPlot(MB, feature = c("Mki67","Sox2"),
        pt.size = 0.1,sort = FALSE, slot = "counts", ncol = 2)
dev.off()
Cairo(file="Violinplot-Mki67 and Sox2 in P21 math1-cre;smoM2 version 2.png",type="png",units="in",bg="white",width=8,height=4,pointsize=114,dpi=300)
VlnPlot(MB, feature = c("Mki67","Sox2"),
        pt.size = 0,sort = FALSE, slot = "counts", ncol = 2)
dev.off()

#Feature-Plot "Mki67" and "Sox2"
Cairo(file="Featureplot-Mki67 and Sox2 in P21 math1-cre;smoM2 version 1.png",type="png",units="in",bg="white",width=18,height=8,pointsize=44,dpi=300)
FeaturePlot(MB, features = c("Mki67","Sox2"),pt.size = 1.5,
            cols = c("White","Red"),
            order = TRUE,
            min.cutoff = NA, max.cutoff = NA, reduction = NULL,
            split.by = NULL, shape.by = NULL, slot = "counts", blend = FALSE, #change "counts" to "data", will get gene expression plot in normalized log data
            blend.threshold = 0, label = TRUE, label.size = 3,
            repel = FALSE, ncol = 2, combine = TRUE, coord.fixed = FALSE,
            by.col = TRUE, sort.cell = FALSE) 
dev.off()
Cairo(file="Featureplot-Mki67 and Sox2 in P21 math1-cre;smoM2 version 2.png",type="png",units="in",bg="white",width=18,height=8,pointsize=44,dpi=300)
FeaturePlot(MB, features = c("Mki67","Sox2"),pt.size = 1.5,
            cols = c("White","Red"),
            order = TRUE,
            min.cutoff = NA, max.cutoff = NA, reduction = NULL,
            split.by = NULL, shape.by = NULL, slot = "counts", blend = FALSE, #change "counts" to "data", will get gene expression plot in normalized log data
            blend.threshold = 0, label = FALSE, label.size = 3,
            repel = FALSE, ncol = 2, combine = TRUE, coord.fixed = FALSE,
            by.col = TRUE, sort.cell = FALSE) 
dev.off()

#Co-expression of "Mki67" and "Sox2"
#Co-labeling of Mki67 and Sox2
gene_exprs <- t(as.matrix(MB@assays$spliced@counts[which(rownames(MB@assays$spliced@counts) %in% c("Mki67","Sox2")),]))
MB@meta.data <- cbind(MB@meta.data,gene_exprs)
#Regrouping Cells based Sox2, Mki67 colabelling
MB@meta.data$Mki67_Grouping <- ifelse(MB@meta.data$Mki67 >= 1,"Mki67_pos","Mki67_neg")
MB@meta.data$Sox2_Grouping <- ifelse(MB@meta.data$Sox2 >= 1,"Sox2_pos","Sox2_neg")
MB@meta.data$Mki67_Sox2_Grouping <- ifelse(MB@meta.data$Mki67 < 1 & MB@meta.data$Sox2 < 1, "Double negtive",
                                            ifelse(MB@meta.data$Mki67 >= 1 & MB@meta.data$Sox2 >= 1,
                                                               "Mki67_pos_Sox2_pos",
                                            ifelse(MB@meta.data$Mki67 >= 1 & MB@meta.data$Sox2 < 1,
                                                                                          "Mki67_pos_Sox2_neg","Mki67_neg_Sox2_pos")))
#Statistics of cell number
table(MB@meta.data$Mki67_Grouping,MB@meta.data$Sox2_Grouping,MB@meta.data$Mki67_Sox2_Grouping)

#Plot colabeling
p1 <- DimPlot(MB, dims = c(1,2), reduction = "umap",
              pt.size = 1.5, split.by = NULL, group.by = "Mki67_Sox2_Grouping", cols = c("Grey","Green","Red","Yellow"),
              shape.by = NULL, order = NULL,
              label = FALSE, label.size = 3) 
Cairo(file="Mki67_Sox2_coexpression_in P21 math1-cre;smoM2 version 1.png",type="png",units="in",bg="white",width=8.5,height=6,pointsize=114,dpi=300)
plot_grid(p1)
dev.off()

#Violin plot of Mki67 and Sox2 in each cell subtype.
a<-table(MB@meta.data$Ident,MB@meta.data$Mki67_Sox2_Grouping)
write.table(a,"Mki67 and Sox2 cell counts P21 math1-cre;smoM2.txt",row.names = TRUE,col.names = TRUE,sep='\t',quote = FALSE)
a<-read.delim("Mki67 and Sox2 cell counts P21 math1-cre;smoM2.txt",header = TRUE)
Cairo(file="Mki67 and Sox2 cell counts P21 math1-cre;smoM2 heatmap.png",type="png",units="in",bg="white",width=5,height=5,pointsize=114,dpi=300)
pheatmap::pheatmap(a,border_color = "black", cellwidth = 15, cellheight = 15,
                   cutree_rows = 1, cutree_cols = 1,
                   legend = TRUE)
dev.off()
#--------------End. P21 Math1-cre;SmoM2 Model----------------

#--------------P15 Math1-cre;SmoM2 Model----------------
setwd("E:/share18/share1/Math1_Cre_SmoM2_model/")
load("integrated_all.Robj")
setwd("E:/share18/share1/Rochelle/P15_Math1_Cre_SmoM2")

#Violin-Plot "Mki67" and "Sox2"
Cairo(file="Violinplot-Mki67 and Sox2 in P15 math1-cre;smoM2 version 1.png",type="png",units="in",bg="white",width=8,height=4,pointsize=114,dpi=300)
VlnPlot(integrated, feature = c("Mki67","Sox2"),
        pt.size = 0.1,sort = FALSE, slot = "counts", ncol = 2)
dev.off()
Cairo(file="Violinplot-Mki67 and Sox2 in P15 math1-cre;smoM2 version 2.png",type="png",units="in",bg="white",width=8,height=4,pointsize=114,dpi=300)
VlnPlot(integrated, feature = c("Mki67","Sox2"),
        pt.size = 0,sort = FALSE, slot = "counts", ncol = 2)
dev.off()

#Feature-Plot "Mki67" and "Sox2"
Cairo(file="Featureplot-Mki67 and Sox2 in P15 math1-cre;smoM2 version 1.png",type="png",units="in",bg="white",width=18,height=8,pointsize=44,dpi=300)
FeaturePlot(integrated, features = c("Mki67","Sox2"),pt.size = 1.5,
            cols = c("White","Red"),
            order = TRUE,
            min.cutoff = NA, max.cutoff = NA, reduction = NULL,
            split.by = NULL, shape.by = NULL, slot = "counts", blend = FALSE, #change "counts" to "data", will get gene expression plot in normalized log data
            blend.threshold = 0, label = TRUE, label.size = 3,
            repel = FALSE, ncol = 2, combine = TRUE, coord.fixed = FALSE,
            by.col = TRUE, sort.cell = FALSE) 
dev.off()
Cairo(file="Featureplot-Mki67 and Sox2 in P15 math1-cre;smoM2 version 2.png",type="png",units="in",bg="white",width=18,height=8,pointsize=44,dpi=300)
FeaturePlot(integrated, features = c("Mki67","Sox2"),pt.size = 1.5,
            cols = c("White","Red"),
            order = TRUE,
            min.cutoff = NA, max.cutoff = NA, reduction = NULL,
            split.by = NULL, shape.by = NULL, slot = "counts", blend = FALSE, #change "counts" to "data", will get gene expression plot in normalized log data
            blend.threshold = 0, label = FALSE, label.size = 3,
            repel = FALSE, ncol = 2, combine = TRUE, coord.fixed = FALSE,
            by.col = TRUE, sort.cell = FALSE) 
dev.off()

#Co-expression of "Mki67" and "Sox2"
#Co-labeling of Mki67 and Sox2
gene_exprs <- t(as.matrix(integrated@assays$RNA@counts[which(rownames(integrated@assays$RNA@counts) %in% c("Mki67","Sox2")),]))
integrated@meta.data <- cbind(integrated@meta.data,gene_exprs)
#Regrouping Cells based Sox2, Mki67 colabelling
integrated@meta.data$Mki67_Grouping <- ifelse(integrated@meta.data$Mki67 >= 1,"Mki67_pos","Mki67_neg")
integrated@meta.data$Sox2_Grouping <- ifelse(integrated@meta.data$Sox2 >= 1,"Sox2_pos","Sox2_neg")
integrated@meta.data$Mki67_Sox2_Grouping <- ifelse(integrated@meta.data$Mki67 < 1 & integrated@meta.data$Sox2 < 1, "Double negtive",
                                            ifelse(integrated@meta.data$Mki67 >= 1 & integrated@meta.data$Sox2 >= 1,
                                                   "Mki67_pos_Sox2_pos",
                                                   ifelse(integrated@meta.data$Mki67 >= 1 & integrated@meta.data$Sox2 < 1,
                                                          "Mki67_pos_Sox2_neg","Mki67_neg_Sox2_pos")))
#Statistics of cell nuintegrateder
table(integrated@meta.data$Mki67_Grouping,integrated@meta.data$Sox2_Grouping,integrated@meta.data$Mki67_Sox2_Grouping)

#Plot colabeling
p1 <- DimPlot(integrated, dims = c(1,2), reduction = "umap",
              pt.size = 1.0, split.by = NULL, group.by = "Mki67_Sox2_Grouping", cols = c("Grey","Green","Red","Yellow"),
              shape.by = NULL, order = NULL,
              label = FALSE, label.size = 3) 
Cairo(file="Mki67_Sox2_coexpression_in P15 math1-cre;smoM2 version 1.png",type="png",units="in",bg="white",width=8.5,height=6,pointsize=114,dpi=300)
plot_grid(p1)
dev.off()

#Violin plot of Mki67 and Sox2 in each cell subtype.
a<-table(integrated@meta.data$Ident,integrated@meta.data$Mki67_Sox2_Grouping)
write.table(a,"Mki67 and Sox2 cell counts P15 math1-cre;smoM2.txt",row.names = TRUE,col.names = TRUE,sep='\t',quote = FALSE)
a<-read.delim("Mki67 and Sox2 cell counts P15 math1-cre;smoM2.txt",header = TRUE)
Cairo(file="Mki67 and Sox2 cell counts P15 math1-cre;smoM2 heatmap.png",type="png",units="in",bg="white",width=5,height=5,pointsize=114,dpi=300)
pheatmap::pheatmap(a,border_color = "black", cellwidth = 15, cellheight = 15,
                   cutree_rows = 1, cutree_cols = 1,
                   legend = TRUE)
dev.off()
#--------------End. P15 Math1-cre;SmoM2 Model----------------

#--------------GFAP-Ptch Model----------------
setwd("E:/share18/share1/Ptch_model/")
load("MB.Robj")
setwd("E:/share18/share1/Rochelle/Ptch")
#Violin-Plot "Mki67" and "Sox2"
Cairo(file="Violinplot-Mki67 and Sox2 in GFAP-Ptch version 1.png",type="png",units="in",bg="white",width=8,height=4,pointsize=114,dpi=300)
VlnPlot(MB, feature = c("Mki67","Sox2"),
        pt.size = 0.1,sort = FALSE, slot = "counts", ncol = 2)
dev.off()
Cairo(file="Violinplot-Mki67 and Sox2 in GFAP-Ptch version 2.png",type="png",units="in",bg="white",width=8,height=4,pointsize=114,dpi=300)
VlnPlot(MB, feature = c("Mki67","Sox2"),
        pt.size = 0,sort = FALSE, slot = "counts", ncol = 2)
dev.off()

#Feature-Plot "Mki67" and "Sox2"
Cairo(file="Featureplot-Mki67 and Sox2 in GFAP-Ptch version 1.png",type="png",units="in",bg="white",width=18,height=8,pointsize=44,dpi=300)
FeaturePlot(MB, features = c("Mki67","Sox2"),pt.size = 1.5,
            cols = c("White","Red"),
            order = TRUE,
            min.cutoff = NA, max.cutoff = NA, reduction = NULL,
            split.by = NULL, shape.by = NULL, slot = "counts", blend = FALSE, #change "counts" to "data", will get gene expression plot in normalized log data
            blend.threshold = 0, label = TRUE, label.size = 3,
            repel = FALSE, ncol = 2, combine = TRUE, coord.fixed = FALSE,
            by.col = TRUE, sort.cell = FALSE) 
dev.off()
Cairo(file="Featureplot-Mki67 and Sox2 in GFAP-Ptch version 2.png",type="png",units="in",bg="white",width=18,height=8,pointsize=44,dpi=300)
FeaturePlot(MB, features = c("Mki67","Sox2"),pt.size = 1.5,
            cols = c("White","Red"),
            order = TRUE,
            min.cutoff = NA, max.cutoff = NA, reduction = NULL,
            split.by = NULL, shape.by = NULL, slot = "counts", blend = FALSE, #change "counts" to "data", will get gene expression plot in normalized log data
            blend.threshold = 0, label = FALSE, label.size = 3,
            repel = FALSE, ncol = 2, combine = TRUE, coord.fixed = FALSE,
            by.col = TRUE, sort.cell = FALSE) 
dev.off()

#Co-expression of "Mki67" and "Sox2"
#Co-labeling of Mki67 and Sox2
gene_exprs <- t(as.matrix(MB@assays$RNA@counts[which(rownames(MB@assays$RNA@counts) %in% c("Mki67","Sox2")),]))
MB@meta.data <- cbind(MB@meta.data,gene_exprs)
#Regrouping Cells based Sox2, Mki67 colabelling
MB@meta.data$Mki67_Grouping <- ifelse(MB@meta.data$Mki67 >= 1,"Mki67_pos","Mki67_neg")
MB@meta.data$Sox2_Grouping <- ifelse(MB@meta.data$Sox2 >= 1,"Sox2_pos","Sox2_neg")
MB@meta.data$Mki67_Sox2_Grouping <- ifelse(MB@meta.data$Mki67 < 1 & MB@meta.data$Sox2 < 1, "Double negtive",
                                            ifelse(MB@meta.data$Mki67 >= 1 & MB@meta.data$Sox2 >= 1,
                                                   "Mki67_pos_Sox2_pos",
                                                   ifelse(MB@meta.data$Mki67 >= 1 & MB@meta.data$Sox2 < 1,
                                                          "Mki67_pos_Sox2_neg","Mki67_neg_Sox2_pos")))
#Statistics of cell number
table(MB@meta.data$Mki67_Grouping,MB@meta.data$Sox2_Grouping,MB@meta.data$Mki67_Sox2_Grouping)

#Plot colabeling
p1 <- DimPlot(MB, dims = c(1,2), reduction = "umap",
              pt.size = 1.5, split.by = NULL, group.by = "Mki67_Sox2_Grouping", cols = c("Grey","Green","Red","Yellow"),
              shape.by = NULL, order = NULL,
              label = FALSE, label.size = 3) 
Cairo(file="Mki67_Sox2_coexpression_in GFAP-Ptch version 1.png",type="png",units="in",bg="white",width=8.5,height=6,pointsize=114,dpi=300)
plot_grid(p1)
dev.off()

#Violin plot of Mki67 and Sox2 in each cell subtype.
a<-table(MB@meta.data$Ident,MB@meta.data$Mki67_Sox2_Grouping)
write.table(a,"Mki67 and Sox2 cell counts GFAP-Ptch.txt",row.names = TRUE,col.names = TRUE,sep='\t',quote = FALSE)
a<-read.delim("Mki67 and Sox2 cell counts GFAP-Ptch.txt",header = TRUE)
Cairo(file="Mki67 and Sox2 cell counts GFAP-Ptchheatmap.png",type="png",units="in",bg="white",width=5,height=5,pointsize=114,dpi=300)
pheatmap::pheatmap(a,border_color = "black", cellwidth = 15, cellheight = 15,
                   cutree_rows = 1, cutree_cols = 1,
                   legend = TRUE)
dev.off()
#--------------End. GFAP-Ptch Model----------------

#--------------Human SHH MB samples----------------
setwd("E:/share18/share1/Patient_sample_from_Northcortt/data/SHH_patients/")
load("integrated_MB.Robj")
setwd("E:/share18/share1/Rochelle/Human_SHH_MB")
Idents(integrated) <- "Ident"


#Violin-Plot "Mki67" and "Sox2"
Cairo(file="Violinplot-MKI67 and SOX2 in Human_SHH_MB_patients version 1.png",type="png",units="in",bg="white",width=8,height=4,pointsize=114,dpi=300)
VlnPlot(integrated, feature = c("MKI67","SOX2"),
        pt.size = 0.1,sort = FALSE, slot = "counts", ncol = 2)
dev.off()
Cairo(file="Violinplot-MKI67 and SOX2 in Human_SHH_MB_patients version 2.png",type="png",units="in",bg="white",width=8,height=4,pointsize=114,dpi=300)
VlnPlot(integrated, feature = c("MKI67","SOX2"),
        pt.size = 0,sort = FALSE, slot = "counts", ncol = 2)
dev.off()

#Feature-Plot "MKI67" and "SOX2"
Cairo(file="Featureplot-MKI67 and SOX2 in Human_SHH_MB_patients version 1.png",type="png",units="in",bg="white",width=18,height=8,pointsize=44,dpi=300)
FeaturePlot(integrated, features = c("MKI67","SOX2"),pt.size = 2.5,
            cols = c("White","Red"),
            order = TRUE,
            min.cutoff = NA, max.cutoff = NA, reduction = NULL,
            split.by = NULL, shape.by = NULL, slot = "counts", blend = FALSE, #change "counts" to "data", will get gene expression plot in normalized log data
            blend.threshold = 0, label = TRUE, label.size = 3,
            repel = FALSE, ncol = 2, combine = TRUE, coord.fixed = FALSE,
            by.col = TRUE, sort.cell = FALSE) 
dev.off()
Cairo(file="Featureplot-MKI67 and SOX2 in Human_SHH_MB_patients version 2.png",type="png",units="in",bg="white",width=18,height=8,pointsize=44,dpi=300)
FeaturePlot(integrated, features = c("MKI67","SOX2"),pt.size = 2.5,
            cols = c("White","Red"),
            order = TRUE,
            min.cutoff = NA, max.cutoff = NA, reduction = NULL,
            split.by = NULL, shape.by = NULL, slot = "counts", blend = FALSE, #change "counts" to "data", will get gene expression plot in normalized log data
            blend.threshold = 0, label = FALSE, label.size = 3,
            repel = FALSE, ncol = 2, combine = TRUE, coord.fixed = FALSE,
            by.col = TRUE, sort.cell = FALSE) 
dev.off()

#Co-expression of "MKI67" and "SOX2"
#Co-labeling of MKI67 and SOX2
gene_exprs <- t(as.matrix(integrated@assays$RNA@counts[which(rownames(integrated@assays$RNA@counts) %in% c("MKI67","SOX2")),]))
integrated@meta.data <- cbind(integrated@meta.data,gene_exprs)
#Regrouping Cells based Sox2, Mki67 colabelling
integrated@meta.data$MKI67_Grouping <- ifelse(integrated@meta.data$MKI67 >= 1,"MKI67_pos","MKI67_neg")
integrated@meta.data$SOX2_Grouping <- ifelse(integrated@meta.data$SOX2 >= 1,"SOX2_pos","SOX2_neg")
integrated@meta.data$MKI67_SOX2_Grouping <- ifelse(integrated@meta.data$MKI67 < 1 & integrated@meta.data$SOX2 < 1, "Double negtive",
                                                    ifelse(integrated@meta.data$MKI67 >= 1 & integrated@meta.data$SOX2 >= 1,
                                                           "MKI67_pos_SOX2_pos",
                                                           ifelse(integrated@meta.data$MKI67 >= 1 & integrated@meta.data$SOX2 < 1,
                                                                  "MKI67_pos_SOX2_neg","MKI67_neg_SOX2_pos")))
#Statistics of cell nuintegrateder
table(integrated@meta.data$MKI67_Grouping,integrated@meta.data$SOX2_Grouping,integrated@meta.data$MKI67_SOX2_Grouping)

#Plot colabeling
p1 <- DimPlot(integrated, dims = c(1,2), reduction = "umap",
              pt.size = 2.0, split.by = NULL, group.by = "MKI67_SOX2_Grouping", cols = c("Grey","Green","Red","Yellow"),
              shape.by = NULL, order = NULL,
              label = FALSE, label.size = 3) 
Cairo(file="MKI67_SOX2_coexpression_in Human_SHH_MB_patients version 1.png",type="png",units="in",bg="white",width=8.5,height=6,pointsize=114,dpi=300)
plot_grid(p1)
dev.off()

#Violin plot of MKI67 and SOX2 in each cell subtype.
a<-table(integrated@meta.data$Ident,integrated@meta.data$MKI67_SOX2_Grouping)
write.table(a,"MKI67 and SOX2 cell counts Human_SHH_MB_patients.txt",row.names = TRUE,col.names = TRUE,sep='\t',quote = FALSE)
a<-read.delim("MKI67 and SOX2 cell counts Human_SHH_MB_patients.txt",header = TRUE)
Cairo(file="MKI67 and SOX2 cell counts Human_SHH_MB_patients heatmap.png",type="png",units="in",bg="white",width=5,height=5,pointsize=114,dpi=300)
pheatmap::pheatmap(a,border_color = "black", cellwidth = 15, cellheight = 15,
                   cutree_rows = 1, cutree_cols = 1,
                   legend = TRUE)
dev.off()
#------------End. Human SHH MB patients----------


