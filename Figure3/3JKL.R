library(Seurat)
library(tidyverse)
library(ggforce)


setwd("/data/gpfs/projects/punim1612/MackayLab/Mackaylab_James_Runx/tony/")

hfhc=Seurat::Read10X(data.dir = 'koda_data/cd8_liver_hfhc/') %>% CreateSeuratObject()
res=Seurat::Read10X(data.dir = 'koda_data/cd8_liver_res/') %>% CreateSeuratObject()
nd=Seurat::Read10X(data.dir = 'koda_data/cd8_liver_nd/') %>% CreateSeuratObject()



x <- merge(hfhc, list(res, nd), 
           add.cell.ids = c("hfhc", "res", "nd"))     
x@meta.data$treatment <- rownames(x@meta.data) %>% str_split_i(pattern = '_', 1)


x[["percent_mito"]] <- PercentageFeatureSet(x, pattern = "^mt")
VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
x <- JoinLayers(x)

x <- x[,x$percent_mito < 10]

x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
x <- FindVariableFeatures(x)
VariableFeatures(x) <- VariableFeatures(x)[-grep(x = VariableFeatures(x), pattern = '(^Ig[klh][vj])|(^Tr[ab][vj])')]
x <- ScaleData(x, features = rownames(x))
x <- RunPCA(x, features = VariableFeatures(x))
x <- FindNeighbors(x, dims = 1:20)
x <- RunUMAP(x, dims = 1:20)
x <- FindClusters(x)


DimPlot(x, label=T)
y <- FindMarkers(x, ident.1 = c(11, 16), min.pct=0.3)


x = x[,!x$seurat_clusters %in% c(11,16)]

x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
x <- FindVariableFeatures(x)
VariableFeatures(x) <- VariableFeatures(x)[-grep(x = VariableFeatures(x), pattern = '(^Ig[klh][vj])|(^Tr[ab][vj])')]
x <- ScaleData(x, features = rownames(x))
x <- RunPCA(x, features = VariableFeatures(x))
x <- FindNeighbors(x, dims = 1:20)
x <- RunUMAP(x, dims = 1:20)
x <- FindClusters(x)

DimPlot(x, label=T)
y <- FindMarkers(x, ident.1 = c(15), min.pct=0.3)


x = x[,!x$seurat_clusters %in% c(15)]

x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
x <- FindVariableFeatures(x)
VariableFeatures(x) <- VariableFeatures(x)[-grep(x = VariableFeatures(x), pattern = '(^Ig[klh][vj])|(^Tr[ab][vj])')]
x <- ScaleData(x, features = rownames(x))
x <- RunPCA(x, features = VariableFeatures(x))
x <- FindNeighbors(x, dims = 1:20)
x <- RunUMAP(x, dims = 1:20)
x <- FindClusters(x, resolution = 0.8)

FeaturePlot(x, c('Klf2','Cxcr6', 'Itgae', 'Pdcd1'), order=F, 
            ncol=2, max.cutoff = 'q99') &
    scale_color_distiller(palette='Reds', direction = 1) & coord_fixed(1) -> p
p
ggsave(p, file='umaps_gene_expr_unordered.pdf', scale = 2)


#####
trm_sig = read.csv('/data/projects/punim1612/signatures/Mackay_sigs/mackay2016_coreSig.txt', sep='\t') %>% filter(Direction == 'up') %>% pull(Symbol)
tex_sig = read.csv('giles_2022', sep='\t') %>% filter(cluster == 'Exh-Term' & avg_log2FC > 0 & p_val_adj < 0.05) %>% pull(gene)

x = AddModuleScore(x, features = list(trm_sig), name = 'trm_2016_sig')
x = AddModuleScore(x, features = list(tex_sig), name = 'tex_2022_sig')


FeaturePlot(x, c('trm_2016_sig1','tex_2022_sig1', 'Tox'), order=T, 
            ncol=2, max.cutoff = 'q99', min.cutoff = 'q1') &
    scale_color_distiller(palette='Reds', direction = 1) & coord_fixed(1) -> p
p
ggsave(p, file='umaps_gene_expr2_upordered.pdf', scale = 2)

x@meta.data$celltype = case_when(x$seurat_clusters %in% c(1,3,7) ~ 'Tn',
          x$seurat_clusters %in% c(9, 5, 15, 13,2,12) ~ 'Tm',
          x$seurat_clusters %in% c(14) ~ 'Cyc',
          x$seurat_clusters %in% c(0,11,10,4,6) ~ 'Tex',
          x$seurat_clusters %in% c(8) ~ 'CD103neg_Trm',
          x$seurat_clusters %in% c(16) ~ 'CD103pos_Trm',
          T~'')


DimPlot(x, label=T, group.by = 'celltype') & coord_fixed(1) -> p
p


y = FindMarkers(x, 16, 8, ) %>% rownames_to_column('gene')

features = c('Itgae',"Il17a", "Il17b", "Il17c", "Il17d", "Il17f", "Il25",
             "Tgfb1", "Tgfb2", "Tgfb3",
             "Csf2",
             "Il6",
             "Il10",
             "Il4",
             "Gzma", "Gzmb", "Gzmc", "Gzmd", "Gzme",
             "Gzmf", "Gzmg", "Gzmk", "Gzmm", "Gzmn")
ggplot(y, aes(x = avg_log2FC, y=-log10(p_val_adj))) + geom_point() + 
    ggrepel::geom_text_repel(data = y %>% filter(gene %in% features), aes(label=gene), color = '#767676') +
    theme_bw() -> p

p %>% ggsave(filename = 'koda_cd103posVsneg_volcano.pdf', scale=3)
DotPlot(x[,x$celltype %in% c('CD103neg_Trm', 'CD103pos_Trm')], 
        features = features, group.by = 'celltype') -> p


p %>% ggsave(filename = 'koda_cd103posVsneg_dotplot.pdf', height = 3, width = 15)
y %>% write.csv('koda_cd103posVsneg_differential.csv')
ggsave(p, file='umaps_gene_clus.pdf', scale = 1)


##### Split #####
x = trms
temp = cbind(x@meta.data, x@reductions$umap@cell.embeddings)


ggplot(temp[,c("umap_1", "umap_2")], aes(x=umap_1, y=umap_2)) + 
    geom_point(color='#bdbdbd') + 
    geom_point(data=temp, color='red') + 
    facet_wrap(~treatment) + theme_bw() -> p0

ggplot(temp[,c("umap_1", "umap_2")], aes(x=umap_1, y=umap_2)) + 
    geom_point(color='#bdbdbd') + 
    geom_point(data=temp, aes(color=celltype)) + 
    facet_wrap(~treatment) + theme_bw() -> p1

ggplot(temp,aes(x=celltype, fill=treatment)) + 
    geom_bar(position='fill', color='black', linewidth = 0.3) + 
    theme_bw() + 
    theme(aspect.ratio = 2/1, 
          axis.text.x = element_text(angle=45, hjust=1)) + 
    scale_y_continuous(expand = c(0,0,0,0.1)) + labs(y='Frequency', x='') -> p2

temp %>% 
    group_by(celltype, treatment) %>% 
    summarise(count = n()) %>% 
    group_by(celltype) %>% mutate(freq = count/sum(count)) -> df_2

ggplot(temp,aes(x=treatment, fill=seurat_clusters)) + 
    geom_bar(position='fill', color='black', linewidth = 0.3) + 
    theme_bw() + 
    theme(aspect.ratio = 2/1, 
          axis.text.x = element_text(angle=45, hjust=1)) + 
    scale_y_continuous(expand = c(0,0,0,0.1)) + labs(y='Frequency', x='') -> p3

temp %>% 
    group_by(treatment, seurat_clusters) %>% 
    summarise(count = n()) %>% 
    group_by(treatment) %>% mutate(freq = count/sum(count)) -> df_3


ggplot(temp,aes(x=treatment, fill=celltype)) + 
    geom_bar(position='fill', color='black', linewidth = 0.3) + 
    theme_bw() + 
    theme(aspect.ratio = 2/1, 
          axis.text.x = element_text(angle=45, hjust=1)) + 
    scale_y_continuous(expand = c(0,0,0,0.1)) + labs(y='Frequency', x='') -> p4

temp %>% 
    group_by(celltype, treatment) %>% 
    summarise(count = n()) %>% 
    group_by(treatment) %>% mutate(freq = count/sum(count)) -> df_4



reprocess = function(seu){
    seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
    seu <- FindVariableFeatures(seu)
    VariableFeatures(seu) <- VariableFeatures(seu)[-grep(x = VariableFeatures(seu), pattern = '(^Ig[klh][vj])|(^Tr[ab][vj])')]
    seu <- ScaleData(seu, features = rownames(seu))
    seu <- RunPCA(seu, features = VariableFeatures(seu))
    seu <- FindNeighbors(seu, dims = 1:20)
    seu <- RunUMAP(seu, dims = 1:20)
    seu <- FindClusters(seu)
}

treatments = c('hfhc','nd','res')

objs = list()

for(treatment in treatments){
    
    objs[[treatment]] = reprocess(seu = x[,x$treatment == treatment])
}

objs %>% 
    lapply(function(x){
        DimPlot(x,group.by='celltype', label=T, label.box = T) + 
            labs(title = unique(x$treatment)) + NoLegend() + coord_fixed()
        }) %>% 
    patchwork::wrap_plots() -> p5
 

ggsave(p0, filename='plots2/umap_split_treatment.pdf', scale=3)
ggsave(p1, filename='plots2/umap_split_treatment_color_celtype.pdf', scale=3)
ggsave(p2, filename='plots2/bar_celltype_x_treatment_fill.pdf', scale=3)
ggsave(p3, filename='plots2/bar_treatment_x_cluster_fill.pdf', scale=3)
ggsave(p4, filename='plots2/bar_treatment_x_celltype_fill.pdf', scale=3)
ggsave(p5, filename='plots/umap_reprocess_split_treatment.pdf', scale=5)

write.csv(df_2, file = 'plots/freqtable_celltype_treatment.csv')
write.csv(df_3, file = 'plots/freqtable_cluster_treatment.csv')



#### END PREPROCESSING ####




liver_runx3_sig <- c('Nr4a1', 'Irf4', 'Fos', 'Fosl2', 'Ifitm2', 'Cdkn1a', 'Hpgds', 
                     'Tnfrsf9', 'Tmem123', 'Cx3cr1', 'Cdh1', 'Slamf7', 'Litaf', 
                     'Il2ra', 'Ccr8', 'Vps37b', 'Runx3', 'Gzmc','Gzmb','Cmah', 'Fgf13',
                     'Adgrg1', 'Tnf', 'Fn1', 'Xcl1', 'Itgae', 'Cd101', 'Rgs2', 'Havcr2', 
                     'Fut11')

liver_tgfb_sig <- c('Chn2', 'Itgae', 'Inpp4b', 'Qpct', 'Litaf', 'Ahr', 'Junb')


runx_sig = read.csv('sigs/Fonseca_R3Sig_up_down_human_mapped.csv')
liver_runx3_sig2 = runx_sig %>% filter(Direction == 'up' ) %>% pull(gene_mouse)
liver_runx3_sig3 = runx_sig %>% filter(Direction == 'down' ) %>% pull(gene_mouse)



tgfb_sig <- read.csv('sigs/Christo_PubTGFsig_up_down_human_mapped.csv')

tgfb_sig_up = tgfb_sig %>% filter(Direction == 'UP') %>% pull(gene_mouse) #%>%  str_to_title()
tgfb_sig_down = tgfb_sig %>% filter(Direction == 'DOWN') %>% pull(gene_mouse) #%>%  str_to_title()

#x <- AddModuleScore(x, features = list(liver_runx3_sig), name='runx3_sig')
#x <- AddModuleScore(x, features = list(liver_tgfb_sig), name='tgfb_sig')

FeaturePlot(x, c('runx3_sig1', 'tgfb_sig1'), max.cutoff = 'q99', pt.size = 1) & scale_color_viridis_c()

trms <- x[,x$celltype %>% str_detect(pattern= 'Trm')]

trms <- NormalizeData(trms, normalization.method = "LogNormalize", scale.factor = 10000)
trms <- FindVariableFeatures(trms)
VariableFeatures(trms) <- VariableFeatures(trms)[-grep(x = VariableFeatures(trms), pattern = '(^Ig[klh][vj])|(^Tr[ab][vj])')]
trms <- ScaleData(trms, features = rownames(trms))
trms <- RunPCA(trms, features = VariableFeatures(trms))
trms <- FindNeighbors(trms, dims = 1:20)
trms <- RunUMAP(trms, dims = 1:20)
trms <- FindClusters(trms, resolution = 0.1)
(DimPlot(trms, pt.size = 1) + coord_fixed()) %>% ggsave(filename='umap_basic.pdf')

DotPlot(trms, features = c('Runx3', 'Itgae', 'Pdcd1', 'Tbx21', 'Il2rb', 'Itgal')) + 
    scale_color_distiller(palette = 'Reds', direction = 1,breaks = c(-1,0,1)) + 
    labs(x='', y='Cluster') + 
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          legend.ticks = element_line(color='black')) -> p
p

trms@meta.data %>% group_by(seurat_clusters, treatment) %>% 
    summarise(count=n()) %>% group_by(seurat_clusters) %>% mutate(freq = count/sum(count)) -> freqs


# Reshape to contingency table
cont_table <- freqs %>%
    select(seurat_clusters, treatment, count) %>%
    pivot_wider(names_from = treatment, values_from = count, values_fill = 0)

# Chi-squared test
chisq.test(cont_table[,-1])


ggsave(p, filename ='plots/dotplot.pdf', scale=2)


trms <- AddModuleScore(trms, features = list(liver_runx3_sig), name='runx3_sig')
trms <- AddModuleScore(trms, features = list(liver_runx3_sig2), name='new_runx3_sig')
trms <- AddModuleScore(trms, features = list(liver_runx3_sig3), name='down_runx3_sig')
trms <- AddModuleScore(trms, features = list(liver_tgfb_sig), name='tgfb_sig')

trms <- AddModuleScore(trms, features = list(tgfb_sig_up), name='up_tgfb_sig')
trms <- AddModuleScore(trms, features = list(tgfb_sig_down), name='down_tgfb_sig')


FeaturePlot(trms, c('runx3_sig1', 'new_runx3_sig1'), max.cutoff = 'q99', pt.size = 1) & 
    scale_color_viridis_c()
trms$merged_runx_sig = trms$new_runx3_sig1 + (-1 * trms$down_runx3_sig1)

VlnPlot(trms, c('new_runx3_sig1','up_tgfb_sig1'), group.by = 'celltype', pt.size = 0) %>% ggsave(filename='violin_basic_2.pdf')

write.csv(x = cbind(trms@meta.data, trms@reductions$umap@cell.embeddings), file = 'trm_meta_new.csv')


x = read.csv('trm_meta_new.csv')


