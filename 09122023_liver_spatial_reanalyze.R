setwd("~/Lab/LiverAtlas/spatial_data/filtered/")
f = list.files(pattern = "-")
library(Seurat)
library(dplyr)
library(ggplot2)

for (i in 1:length(f)) {
    cat("filtering", f[i], "..\n")
    
    st = Load10X_Spatial(f[i])
    
    genes = rownames(st@assays$Spatial@counts)
    mt_rps_genes = grep("^MT-|^RPS", genes, value = T)
    low_expr_genes = apply(st@assays$Spatial@counts, 1, function(x) {
        sum(x > 0)
    })
    low_expr_genes = names(low_expr_genes[low_expr_genes < 3])
    rm_genes = union(mt_rps_genes, low_expr_genes)
    
    filterd_genes = genes[! genes %in% rm_genes]
    filtered_st = subset(st, features = filterd_genes)
    filtered_st@meta.data$Study = "Wu.SciAdv.2021"
    filtered_st@meta.data$Cancer_Subtype = unlist(strsplit(f[i], "-"))[1]
    filtered_st@meta.data$Patient = substr(f[i], 1, nchar(f[i])-1)
    filtered_st@meta.data$Tissue_Section = substr(f[i], nchar(f[i]), nchar(f[i]))
    saveRDS(filtered_st, paste0(f[i], "/", f[i], "_filtered_srt.rds"))
}

patient = unique(substr(f, 1, 5))

for (p in patient) {
    cat("processing", p, "..\n")
    
    samples = grep(p, f, value = T)
    srt_lst = lapply(paste0(samples, "/", samples, "_filtered_srt.rds"), 
                     readRDS)
    names(srt_lst) = samples
    srt = Reduce(function(x, y) merge(x, y), x = srt_lst)
    cnt = srt@assays$Spatial@counts
    meta = srt@meta.data
    srt = CreateSeuratObject(counts = cnt, meta.data = meta)
    
    srt = srt %>%
        NormalizeData() %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA()
    
    if (length(samples) > 1) {
        srt = harmony::RunHarmony(srt, group.by.vars = "Tissue_Section")
        srt = srt %>%
            FindNeighbors(reduction = "harmony", dims = 1:20) %>%
            FindClusters(resolution = seq(0.1, 0.5, 0.05)) %>%
            RunUMAP(reduction = "harmony", dims = 1:20)
    } else {
        srt = srt %>%
            FindNeighbors(reduction = "pca", dims = 1:20) %>%
            FindClusters(resolution = seq(0.1, 0.5, 0.05)) %>%
            RunUMAP(reduction = "pca", dims = 1:20)
    }
    saveRDS(srt, paste0("../patient/", p, "_processed_srt.rds"))
}


author_col = c("#F5C95D", "#517AF6", "#82FB7C", "#EC724F", 
                        "#DA8FC0", "#73D6FA", "#FCED53", "#B0D767")
names(author_col) = 1:8

# ---- cHC-1 ----
p = "cHC-1"
res = 0.3
res_var = paste0("RNA_snn_res.", res)

srt = readRDS(paste0("../patient/", p, "_processed_srt.rds"))
srt@meta.data$Tissue_Section = factor(srt@meta.data$Tissue_Section, 
                                      levels = c("N", "L", "T"))
srt@meta.data$Sample = paste0(srt@meta.data$Patient, srt@meta.data$Tissue_Section)
DimPlot(srt, split.by = "Tissue_Section", group = res_var)

samples = grep(p, f, value = T)
srt_lst = lapply(paste0(samples, "/", samples, "_filtered_srt.rds"), 
                 readRDS)
names(srt_lst) = samples
for (s in samples) {
    srt_lst[[s]]@meta.data$cluster_patient = srt@meta.data[[res_var]][srt@meta.data$Sample == s]
}

for (s in samples) {
    cluster_author = srt_lst[[s]]@meta.data$cluster_patient
    
    srt_lst[[s]]@meta.data$cluster_author = cluster_author
}
SpatialDimPlot(srt_lst$`cHC-1N`, group.by = "cluster_patient") + NoLegend()
SpatialDimPlot(srt_lst$`cHC-1L`, group.by = "cluster_patient") + NoLegend()
SpatialDimPlot(srt_lst$`cHC-1T`, group.by = "cluster_patient") + NoLegend()

cluster_author = srt@meta.data[[res_var]]
mapping = data.frame(
    patient = c(0, 1, 2, 3, 4, 5),
    author  = c(1, 3, 4, 2, 5, 3)
)
cluster_author = as.factor(mapping$author[match(cluster_author, mapping$patient)])
for (s in samples) {
    srt_lst[[s]]@meta.data$cluster_author = cluster_author[srt@meta.data$Sample == s]
}
SpatialDimPlot(srt_lst$`cHC-1N`, group.by = "cluster_author", image.alpha = 0, stroke = NA) + 
    NoLegend() + 
    scale_fill_manual(values = author_col[unique(srt_lst$`cHC-1N`@meta.data$cluster_author)])
SpatialDimPlot(srt_lst$`cHC-1L`, group.by = "cluster_author", image.alpha = 0, stroke = NA) + 
    NoLegend() + 
    scale_fill_manual(values = author_col[unique(srt_lst$`cHC-1L`@meta.data$cluster_author)])
SpatialDimPlot(srt_lst$`cHC-1T`, group.by = "cluster_author", image.alpha = 0, stroke = NA) + 
    NoLegend() + 
    scale_fill_manual(values = author_col[unique(srt_lst$`cHC-1T`@meta.data$cluster_author)])
for (s in samples) {
    saveRDS(srt_lst[[s]], paste0(s, "/", s, "_filtered_srt.rds"))
}


# ---- ICC-1 ----
p = "ICC-1"
res = 0.3
res_var = paste0("RNA_snn_res.", res)

srt = readRDS(paste0("../patient/", p, "_processed_srt.rds"))
srt@meta.data$Tissue_Section = factor(srt@meta.data$Tissue_Section, 
                                      levels = c("N", "L", "T"))
srt@meta.data$Sample = paste0(srt@meta.data$Patient, srt@meta.data$Tissue_Section)
DimPlot(srt, split.by = "Tissue_Section", group = res_var)

samples = grep(p, f, value = T)
srt_lst = lapply(paste0(samples, "/", samples, "_filtered_srt.rds"), 
                 readRDS)
names(srt_lst) = samples
for (s in samples) {
    srt_lst[[s]]@meta.data$cluster_patient = srt@meta.data[[res_var]][srt@meta.data$Sample == s]
}

for (s in samples) {
    cluster_author = srt_lst[[s]]@meta.data$cluster_patient
    
    srt_lst[[s]]@meta.data$cluster_author = cluster_author
}
SpatialDimPlot(srt_lst$`ICC-1L`, group.by = "cluster_patient") + NoLegend()

cluster_author = srt@meta.data[[res_var]]
mapping = data.frame(
    patient = c(0, 1, 2, 3),
    author  = c(1, 2, 4, 3)
)
cluster_author = as.factor(mapping$author[match(cluster_author, mapping$patient)])
for (s in samples) {
    srt_lst[[s]]@meta.data$cluster_author = cluster_author[srt@meta.data$Sample == s]
}
SpatialDimPlot(srt_lst$`ICC-1L`, group.by = "cluster_author", image.alpha = 0, stroke = NA) + 
    NoLegend() + 
    scale_fill_manual(values = author_col[unique(srt_lst$`ICC-1L`@meta.data$cluster_author)])
for (s in samples) {
    saveRDS(srt_lst[[s]], paste0(s, "/", s, "_filtered_srt.rds"))
}


# ---- HCC-1 ----
p = "HCC-1"
res = 0.3
res_var = paste0("RNA_snn_res.", res)

srt = readRDS(paste0("../patient/", p, "_processed_srt.rds"))
srt@meta.data$Tissue_Section = factor(srt@meta.data$Tissue_Section, 
                                      levels = c("N", "L", "T"))
srt@meta.data$Sample = paste0(srt@meta.data$Patient, srt@meta.data$Tissue_Section)
DimPlot(srt, split.by = "Tissue_Section", group = res_var)

samples = grep(p, f, value = T)
srt_lst = lapply(paste0(samples, "/", samples, "_filtered_srt.rds"), 
                 readRDS)
names(srt_lst) = samples
for (s in samples) {
    srt_lst[[s]]@meta.data$cluster_patient = srt@meta.data[[res_var]][srt@meta.data$Sample == s]
}

for (s in samples) {
    cluster_author = srt_lst[[s]]@meta.data$cluster_patient
    
    srt_lst[[s]]@meta.data$cluster_author = cluster_author
}
SpatialDimPlot(srt_lst$`HCC-1N`, group.by = "cluster_patient") + NoLegend()
SpatialDimPlot(srt_lst$`HCC-1L`, group.by = "cluster_patient") + NoLegend()
SpatialDimPlot(srt_lst$`HCC-1T`, group.by = "cluster_patient") + NoLegend()

cluster_author = srt@meta.data[[res_var]]
mapping = data.frame(
    patient = c(0, 1, 2, 3, 4, 5, 6, 7),
    author  = c(1, 2, 3, 4, 5, 6, 7, 3)
)
cluster_author = as.factor(mapping$author[match(cluster_author, mapping$patient)])
for (s in samples) {
    srt_lst[[s]]@meta.data$cluster_author = cluster_author[srt@meta.data$Sample == s]
}
SpatialDimPlot(srt_lst$`HCC-1N`, group.by = "cluster_author", image.alpha = 0, stroke = NA) + 
    NoLegend() + 
    scale_fill_manual(values = author_col[unique(srt_lst$`HCC-1N`@meta.data$cluster_author)])
SpatialDimPlot(srt_lst$`HCC-1L`, group.by = "cluster_author", image.alpha = 0, stroke = NA) + 
    NoLegend() + 
    scale_fill_manual(values = author_col[unique(srt_lst$`HCC-1L`@meta.data$cluster_author)])
SpatialDimPlot(srt_lst$`HCC-1T`, group.by = "cluster_author", image.alpha = 0, stroke = NA) + 
    NoLegend() + 
    scale_fill_manual(values = author_col[unique(srt_lst$`HCC-1T`@meta.data$cluster_author)])
for (s in samples) {
    saveRDS(srt_lst[[s]], paste0(s, "/", s, "_filtered_srt.rds"))
}



# ---- HCC-2 ----
p = "HCC-2"
res = 0.3
res_var = paste0("RNA_snn_res.", res)

srt = readRDS(paste0("../patient/", p, "_processed_srt.rds"))
srt@meta.data$Tissue_Section = factor(srt@meta.data$Tissue_Section, 
                                      levels = c("N", "L", "T", "P"))
srt@meta.data$Sample = paste0(srt@meta.data$Patient, srt@meta.data$Tissue_Section)
DimPlot(srt, split.by = "Tissue_Section", group = res_var)

samples = grep(p, f, value = T)
srt_lst = lapply(paste0(samples, "/", samples, "_filtered_srt.rds"), 
                 readRDS)
names(srt_lst) = samples
for (s in samples) {
    srt_lst[[s]]@meta.data$cluster_patient = srt@meta.data[[res_var]][srt@meta.data$Sample == s]
}

for (s in samples) {
    cluster_author = srt_lst[[s]]@meta.data$cluster_patient
    
    srt_lst[[s]]@meta.data$cluster_author = cluster_author
}
SpatialDimPlot(srt_lst$`HCC-2N`, group.by = "cluster_patient") + NoLegend()
SpatialDimPlot(srt_lst$`HCC-2L`, group.by = "cluster_patient") + NoLegend()
SpatialDimPlot(srt_lst$`HCC-2T`, group.by = "cluster_patient") + NoLegend()
SpatialDimPlot(srt_lst$`HCC-2P`, group.by = "cluster_patient") + NoLegend()

cluster_author = srt@meta.data[[res_var]]
mapping = data.frame(
    patient = c(0, 1, 2, 3, 4, 5, 6, 7, 8),
    author  = c(1, 5, 4, 2, 3, 2, 7, 8, 6)
)
cluster_author = as.factor(mapping$author[match(cluster_author, mapping$patient)])
for (s in samples) {
    srt_lst[[s]]@meta.data$cluster_author = cluster_author[srt@meta.data$Sample == s]
}
SpatialDimPlot(srt_lst$`HCC-2N`, group.by = "cluster_author", image.alpha = 0, stroke = NA) + 
    NoLegend() + 
    scale_fill_manual(values = author_col[unique(srt_lst$`HCC-2N`@meta.data$cluster_author)])
SpatialDimPlot(srt_lst$`HCC-2L`, group.by = "cluster_author", image.alpha = 0, stroke = NA) + 
    NoLegend() + 
    scale_fill_manual(values = author_col[unique(srt_lst$`HCC-2L`@meta.data$cluster_author)])
SpatialDimPlot(srt_lst$`HCC-2T`, group.by = "cluster_author", image.alpha = 0, stroke = NA) + 
    NoLegend() + 
    scale_fill_manual(values = author_col[unique(srt_lst$`HCC-2T`@meta.data$cluster_author)])
SpatialDimPlot(srt_lst$`HCC-2P`, group.by = "cluster_author", image.alpha = 0, stroke = NA) + 
    NoLegend() + 
    scale_fill_manual(values = author_col[unique(srt_lst$`HCC-2P`@meta.data$cluster_author)])
for (s in samples) {
    saveRDS(srt_lst[[s]], paste0(s, "/", s, "_filtered_srt.rds"))
}



# ---- HCC-3 ----
p = "HCC-3"
res = 0.3
res_var = paste0("RNA_snn_res.", res)

srt = readRDS(paste0("../patient/", p, "_processed_srt.rds"))
srt@meta.data$Tissue_Section = factor(srt@meta.data$Tissue_Section, 
                                      levels = c("N", "L", "T"))
srt@meta.data$Sample = paste0(srt@meta.data$Patient, srt@meta.data$Tissue_Section)
DimPlot(srt, split.by = "Tissue_Section", group = res_var)

samples = grep(p, f, value = T)
srt_lst = lapply(paste0(samples, "/", samples, "_filtered_srt.rds"), 
                 readRDS)
names(srt_lst) = samples
for (s in samples) {
    srt_lst[[s]]@meta.data$cluster_patient = srt@meta.data[[res_var]][srt@meta.data$Sample == s]
}

for (s in samples) {
    cluster_author = srt_lst[[s]]@meta.data$cluster_patient
    
    srt_lst[[s]]@meta.data$cluster_author = cluster_author
}
SpatialDimPlot(srt_lst$`HCC-3N`, group.by = "cluster_patient") + NoLegend()
SpatialDimPlot(srt_lst$`HCC-3L`, group.by = "cluster_patient") + NoLegend()
SpatialDimPlot(srt_lst$`HCC-3T`, group.by = "cluster_patient") + NoLegend()

cluster_author = srt@meta.data[[res_var]]
mapping = data.frame(
    patient = c(0, 1, 2, 3, 4, 5, 6, 7),
    author  = c(1, 2, 4, 3, 6, 5, 7, 8)
)
cluster_author = as.factor(mapping$author[match(cluster_author, mapping$patient)])
for (s in samples) {
    srt_lst[[s]]@meta.data$cluster_author = cluster_author[srt@meta.data$Sample == s]
}
SpatialDimPlot(srt_lst$`HCC-3N`, group.by = "cluster_author", image.alpha = 0, stroke = NA) + 
    NoLegend() + 
    scale_fill_manual(values = author_col[unique(srt_lst$`HCC-3N`@meta.data$cluster_author)])
SpatialDimPlot(srt_lst$`HCC-3L`, group.by = "cluster_author", image.alpha = 0, stroke = NA) + 
    NoLegend() + 
    scale_fill_manual(values = author_col[unique(srt_lst$`HCC-3L`@meta.data$cluster_author)])
SpatialDimPlot(srt_lst$`HCC-3T`, group.by = "cluster_author", image.alpha = 0, stroke = NA) + 
    NoLegend() + 
    scale_fill_manual(values = author_col[unique(srt_lst$`HCC-3T`@meta.data$cluster_author)])
for (s in samples) {
    saveRDS(srt_lst[[s]], paste0(s, "/", s, "_filtered_srt.rds"))
}


# ---- HCC-4 ----
p = "HCC-4"
res = 0.3
res_var = paste0("RNA_snn_res.", res)

srt = readRDS(paste0("../patient/", p, "_processed_srt.rds"))
srt@meta.data$Tissue_Section = factor(srt@meta.data$Tissue_Section, 
                                      levels = c("N", "L", "T"))
srt@meta.data$Sample = paste0(srt@meta.data$Patient, srt@meta.data$Tissue_Section)
DimPlot(srt, split.by = "Tissue_Section", group = res_var)

samples = grep(p, f, value = T)
srt_lst = lapply(paste0(samples, "/", samples, "_filtered_srt.rds"), 
                 readRDS)
names(srt_lst) = samples
for (s in samples) {
    srt_lst[[s]]@meta.data$cluster_patient = srt@meta.data[[res_var]][srt@meta.data$Sample == s]
}

for (s in samples) {
    cluster_author = srt_lst[[s]]@meta.data$cluster_patient
    
    srt_lst[[s]]@meta.data$cluster_author = cluster_author
}
SpatialDimPlot(srt_lst$`HCC-4N`, group.by = "cluster_patient") + NoLegend()
SpatialDimPlot(srt_lst$`HCC-4L`, group.by = "cluster_patient") + NoLegend()
SpatialDimPlot(srt_lst$`HCC-4T`, group.by = "cluster_patient") + NoLegend()

cluster_author = srt@meta.data[[res_var]]
mapping = data.frame(
    patient = c(0, 1, 2, 3, 4),
    author  = c(1, 2, 3, 4, 5)
)
cluster_author = as.factor(mapping$author[match(cluster_author, mapping$patient)])
for (s in samples) {
    srt_lst[[s]]@meta.data$cluster_author = cluster_author[srt@meta.data$Sample == s]
}
SpatialDimPlot(srt_lst$`HCC-4N`, group.by = "cluster_author", image.alpha = 0, stroke = NA) + 
    NoLegend() + 
    scale_fill_manual(values = author_col[unique(srt_lst$`HCC-4N`@meta.data$cluster_author)])
SpatialDimPlot(srt_lst$`HCC-4L`, group.by = "cluster_author", image.alpha = 0, stroke = NA) + 
    NoLegend() + 
    scale_fill_manual(values = author_col[unique(srt_lst$`HCC-4L`@meta.data$cluster_author)])
SpatialDimPlot(srt_lst$`HCC-4T`, group.by = "cluster_author", image.alpha = 0, stroke = NA) + 
    NoLegend() + 
    scale_fill_manual(values = author_col[unique(srt_lst$`HCC-4T`@meta.data$cluster_author)])
for (s in samples) {
    saveRDS(srt_lst[[s]], paste0(s, "/", s, "_filtered_srt.rds"))
}

# ---- all samples ----
srt_lst = lapply(paste0(f, "/", f, "_filtered_srt.rds"), 
                 readRDS)
names(srt_lst) = f
srt = Reduce(function(x, y) merge(x, y), x = srt_lst)
cnt = srt@assays$Spatial@counts
meta = srt@meta.data
meta$Sample = paste0(meta$Patient, meta$Tissue_Section)
srt = CreateSeuratObject(counts = cnt, meta.data = meta)
srt = srt %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    harmony::RunHarmony(group.by.vars = c("Sample")) %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = seq(0.2, 0.8, 0.1)) %>%
    RunUMAP(reduction = "harmony", dims = 1:30)
saveRDS(srt, paste0("../all_samples_correct_sample_srt.rds"))

srt = readRDS("../all_samples_correct_patient_section_srt.rds")
DimPlot(srt, reduction = "umap", group.by = "Patient", split.by = "Tissue_Section")