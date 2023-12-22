options(Seurat.object.assay.version = "v3") 
 
#' @title add_clonotype 
#' @author Jiahao Wang 
#' @export
add_clonotype <- function(tcr_prefix, seurat_obj, type="t"){ 
    tcr <- read.csv(paste0(tcr_prefix, "/filtered_contig_annotations.csv")) 
    tcr <- tcr[!duplicated(tcr$barcode), ] 
    tcr <- tcr[,c("barcode", "raw_clonotype_id")] 
    names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id" 
 
    clono <- read.csv(paste0(tcr_prefix, "/clonotypes.csv")) 
    tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")]) 
    names(tcr)[names(tcr) == "cdr3s_aa"] <- "cdr3s_aa" 
 
    tcr <- tcr[, c(2,1,3)] 
    rownames(tcr) <- tcr[,1] 
    tcr[,1] <- NULL 
    colnames(tcr) <- paste(type, colnames(tcr), sep="_") 
    clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr) 
    return(clono_seurat) 
} 
 
#' @title extractVDJ 
#' @author Jiahao Wang 
#' @export
extractVDJ <- function(vdjOutDir, type = NULL){ 
 
    if(is.null(type)){ 
        stop("Please specific 'type', t/b.") 
    } 
    contigAnno = read.csv(paste0(vdjOutDir, "/filtered_contig_annotations.csv")) 
    contigAnno = contigAnno[!duplicated(contigAnno$barcode), ] 
    contigAnno = contigAnno[,c("barcode", "raw_clonotype_id")] 
    names(contigAnno)[names(contigAnno) == "raw_clonotype_id"] <- "clonotype_id" 
 
    clonoTypes = read.csv(paste0(vdjOutDir, "/clonotypes.csv")) 
    contigAnno = merge(contigAnno, clonoTypes[, c("clonotype_id", "cdr3s_aa")]) 
 
    contigAnno = contigAnno[, c(2, 1, 3)] 
    rownames(contigAnno) = contigAnno[, 1] 
    contigAnno[, 1] = NULL 
    colnames(contigAnno) <- paste0(type, "_", colnames(contigAnno)) 
    return(contigAnno) 
} 
 
#' @title readr 
#' @author Jiahao Wang 
#' @export
readr <- function(file_path){ 
    res = get(load2(file_path)) 
    return(res) 
} 
 
#' @title readRData 
#' @author Jiahao Wang 
#' @export
readRData <- function(file_path, idx = NULL, rm = FALSE, verbose = FALSE){ 
    obj_name = load(file_path) 
 
    if(is.null(idx) & length(obj_name) > 1){ 
            cat0("This RData contains multiple objects: ", paste0(obj_name, collapse = ", "), "") 
            cat0("Default return the first object, use 'idx' to select manully.")            
    } 
 
    if(is.null(idx)){ 
        idx = 1 
    } 
 
    if(length(idx) > 1){ 
        if(verbose){ 
            cat0("Loadding object: ", paste0(obj_name[idx], collapse = " "))     
        } 
        res = sapply(idx, function(x) get(obj_name[x]), simplify = FALSE) 
        names(res) = obj_name[idx] 
    } else if(idx == "all"){ 
        idx = 1:length(obj_name) 
        res = sapply(idx, function(x) get(obj_name[x]), simplify = FALSE) 
        names(res) = obj_name[idx] 
    } else {         
        if(verbose){ 
            cat0("Loadding object: ", paste0(obj_name[idx], collapse = " "))     
        } 
        res = get(obj_name[idx]) 
    } 
    return(res) 
} 
 
# to do: cluster可指定多个 
#' @title sc_subset_obj 
#' @author Jiahao Wang 
#' @export
sc_subset_obj <- function(obj, cluster, group = "seurat_clusters", verbose = TRUE){ 
    if(!group %in% colnames(obj@meta.data)){ 
        info = paste0("Given group '", group, "' can not be found in meta data, please check it!") 
        report(info, "E") 
        stop() 
    } 
 
    if(!cluster %in% obj@meta.data[, group]){ 
        info = paste0("Given cluster '", cluster, "' can not be found in specific group '", group, "', please check it!") 
        report(info, "E") 
        stop()        
    } 
    xobj = obj[, obj@meta.data[, group] == cluster] 
    if(verbose){ 
        report(paste0("Subset ", ncol(xobj), "/", ncol(obj), " cells from cluster '", cluster, "' in group '", group, "'")) 
    }     
 
    return(xobj) 
} 
 
#' @title sc_get_count 
#' @author Jiahao Wang 
#' @export
sc_get_count <- function(obj, cluster = NULL, group = "seurat_clusters", verbose = TRUE){ 
    if(!is.null(cluster)){ 
        cluster %in% obj@meta.data[, group] 
        obj = sc_subset_obj(obj, cluster, group, verbose = verbose) 
    } 
    count_tb = obj@assays$RNA@counts 
    return(count_tb) 
} 
#' @title sc_VolcanoPlot 
#' @author Jiahao Wang 
#' @export
sc_VolcanoPlot <- function(diff,  
                           log2FC = log2(1.5),  
                           padj = 0.05,  
                           label.symbols = NULL,  
                           label.max = 30, 
                           cols = c("#F93463", "#27A2EF"),  
                           title = ""){ 
    if( all( !c("log2FoldChange", "padj", "symbol") %in% colnames(diff) )){ 
    stop("Colnames must include: log2FoldChange, padj, symbol") 
    } 
    rownames(diff) = diff$symbol 
 
    # (1) define up and down 
    diff$threshold="ns"; 
    if(length(which(diff$log2FoldChange > log2FC & diff$padj < padj)) != 0){ 
        diff[which(diff$log2FoldChange > log2FC & diff$padj <padj),]$threshold="up" 
    } 
 
    if(length(which(diff$log2FoldChange < (-log2FC) & diff$padj < padj)) != 0){ 
        diff[which(diff$log2FoldChange < (-log2FC) & diff$padj < padj),]$threshold="down"; 
    } 
 
    diff$threshold=factor(diff$threshold, levels=c('down','ns','up')) 
    #head(diff) 
    # 
    tb2=table(diff$threshold)# ; print(tb2) 
    library(ggplot2) 
    # (2) plot 
    g1 = ggplot(data=diff, aes(x=log2FoldChange, y=-log10(padj), color=threshold)) + 
    geom_point(alpha=0.8, size=1.5) + 
    geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+ 
    geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+ 
    labs(title= ifelse(""==title, "", paste(title)))+ 
    xlab(bquote(Log[2]*FoldChange))+ 
    ylab(bquote(-Log[10]*italic(P.adj)) )+ 
    theme_classic(base_size = 14) + 
    theme(legend.box = "horizontal", 
          legend.position="top", 
          legend.spacing.x = unit(0, 'pt'), 
          legend.text = element_text( margin = margin(r = 20) ), 
          legend.margin=margin(b= -10, unit = "pt"), 
          plot.title = element_text(hjust = 0.5, size=10) 
          ) + 
    coord_flip() + setText(20) 
 
    no_up = length(which(diff$log2FoldChange > log2FC & diff$padj < padj)) == 0 
    no_dn = length(which(diff$log2FoldChange < (-log2FC) & diff$padj < padj)) == 0 
 
    if(no_up){ 
        g1 <- g1 + scale_color_manual('',labels=c(paste0("down(",tb2[[1]],')'),'ns'), 
                           values=c(cols[1], "grey") ) 
    } 
 
    if(no_dn){ 
        g1 <- g1 + scale_color_manual(labels=c('ns', 
                                           paste0("up(",tb2[[3]],')' )), 
                               values=c("grey", cols[2]) ) 
    } 
 
    if(!no_up & !no_dn){ 
        g1 <- g1 + scale_color_manual('',labels=c(paste0("down(",tb2[[1]],')'),'ns', 
                                       paste0("up(",tb2[[3]],')' )), 
                           values=c(cols[1], "grey", cols[2]) ) 
    }         
 
    g1 <- g1 + guides(color=guide_legend(override.aes = list(size=3, alpha=1))) 
 
    # (3)label genes 
    if(is.null(label.symbols)){ 
    diff.sig=diff[which(diff$threshold != "ns" ), ] 
    len=nrow(diff.sig) 
    if(len<label.max){ 
      label.symbols=rownames(diff.sig) 
    }else{ 
      diff.sig=diff.sig[order(diff.sig$log2FoldChange), ] 
      diff.sig= rbind(diff.sig[1:(label.max/2),], diff.sig[(len-label.max/2):len,]) 
      label.symbols=rownames(diff.sig) 
    } 
    } 
    dd_text = diff[label.symbols, ] 
    # print((dd_text)) 
    # add text 
    library(ggrepel) 
    p <- g1 + geom_text_repel(data=dd_text, 
                       aes(x=log2FoldChange, y=-log10(padj), label=row.names(dd_text)), 
                         #size=2.5,  
                         colour="black",alpha=1, max.overlaps = 50) 
    return(p) 
} 
 
 
#' @title sc_merge_qc 
#' @author Jiahao Wang 
#' @export
sc_merge_qc <- function(input_dir, type = "count"){ 
    if(type %!in% c("count", "multi_5_Only", "multi_5_PE")){ 
        cat("Only support type: count|multi_5_Only|multi_5_PE!\n") 
        stop() 
    } 
 
    key_list <- list( 
        count = c( 
            "Estimated Number of Cells", 
            "Median Genes per Cell", 
            "Number of Reads", 
            "Valid Barcodes", 
            "Sequencing Saturation", 
            "Reads Mapped to Genome", 
            "Reads Mapped Confidently to Genome", 
            "Fraction Reads in Cells", 
            "Total Genes Detected", 
            "Median UMI Counts per Cell"), 
        multi_5_Only = c(40, 28, 42, 6, 4, 41, 36, 2, 45, 23, 26, 64, 72, 13, 17, 48, 56), 
        multi_5_PE   = c(41, 44, 43, 6, 4, 42, 37, 2, 46, 23, 26, 65, 73, 13, 17, 49, 57) 
    ) 
 
    tags = lf("x", input_dir) 
    qc_table = data.frame() 
    if(type == "count"){ 
        cat0(type, "\n") 
        for(tag in tags){ 
            xqc_table = read.csv(paste0(input_dir, "/", tag, "/outs/metrics_summary.csv"), check.names = FALSE, header = FALSE) 
            print(xqc_table) 
            rownames(xqc_table) = xqc_table[, 1] 
            xqc_table = cbind(Sample = tag, xqc_table[key_list[[type]], ]) 
            qc_table = rbind(qc_table, xqc_table) 
        } 
    } else { 
        for(tag in tags){ 
            xqc_table = read.csv(paste0(input_dir, "/", tag, "/outs/per_sample_outs/", tag, "/metrics_summary.csv"), check.names = FALSE) 
            xqc_table = xqc_table[key_list[[type]], c(2, 5:6)] 
            xqc_table = cbind(Sample = tag, xqc_table) 
            qc_table = rbind(qc_table, xqc_table) 
        } 
    } 
    return(qc_table) 
} 
 
 
#' @title sc_dotplot_general 
#' @author Jiahao Wang 
#' @export
sc_dotplot_general <- function(obj, group, file){ 
    celltypes = c("B", "Endo", "Epis", "Fibro", "Myeloid", "T_I_NK") 
    markers = unlist(sapply(celltypes, selectMarkers)) 
    names(markers) = NULL 
     
    obj = ScaleData(obj, features = markers) 
    p <- DotPlot(obj, group.by = group, features = markers) + coord_flip() 
    p <- p + setText(18) 
    p <- p + theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1,vjust = 0.5)) 
    p <- p + labs(x = NULL, y = NULL) 
    p <- p + scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#bebcbe', '#E71908')) 
    p <- p + guides(color = guide_colorbar(barwidth = unit(4.5, "cm"), title = "Average Expression")) 
    p <- p + theme(legend.position = "top") 
    p <- p + theme(axis.text.x = element_text(size = 15 * 0.6, hjust = 0.5)) 
    save.graph(p, file = file, 10, 7, PDF = FALSE) 
} 
 
#' @title cluster_start_with_one 
#' @author Jiahao Wang 
#' @export
cluster_start_with_one <- function(obj, cluster_name = "seurat_clusters"){ 
    if(sum(obj[[cluster_name]][, 1] == 0) == 0){ 
        cat("Maybe you have run this function before!\n") 
    } else { 
        cluster = as.numeric(as.character(obj[[cluster_name]][, 1])) + 1 
        cluster = factor(cluster, levels = c(1:length(unique(cluster)))) 
        obj[[cluster_name]] = cluster 
    } 
    return(obj)         
} 
 
 
#' @title prepareCounts 
#' @author Jiahao Wang 
#' @export
prepareCounts <- function(outs){ 
    counts = Seurat::Read10X(outs) 
    features = data.table::fread(paste0(outs, "/features.tsv.gz"), data.table = FALSE, header = FALSE) 
    counts = suppressWarnings(cbind(features[, -3], counts)) 
    colnames(counts)[1:2] = c("GeneID", "GeneName") 
    write("c", counts, file = paste0(outs, "/../counts.txt")) 
    write("c", counts[1:100, 1:102], file = paste0(outs, "/../counts_100_x_100.txt")) 
    invisible(counts) 
} 
 
#' @title FindClusters2 
#' @author Jiahao Wang 
#' @export
FindClusters2 <- function(obj, cluster.range = c(35, 40), by = 0.1, res = 1, verbose = FALSE){ 
    if(verbose) 
        cat("Find suitable resolution, start with", res, "\n") 
 
    if(length(cluster.range) == 1){ 
        cluster.range = c(cluster.range, cluster.range) 
    } 
 
    sigmoid <- function(x) 1 / (1 + exp(-x)) 
    x = -log(10/res - 1) 
    plusCounter = minusCounter = 0 
    nCluster = 1 
    while(nCluster < cluster.range[1] | nCluster > cluster.range[2]){ 
        if(verbose) 
            cat("resolution", res, "... ") 
        obj <- FindClusters(obj, resolution = res, verbose = FALSE) 
        nCluster = length(unique(obj$seurat_clusters)) 
 
        if(nCluster < cluster.range[1]){ 
            x = x + by 
            plusCounter = plusCounter + 1 
        } else if (nCluster > cluster.range[2]){ 
            x = x - by 
            minusCounter = minusCounter + 1 
        } else{             
            break 
        } 
        res = signif(sigmoid(x) * 10, 3) 
        if(plusCounter & minusCounter){ 
            cat("\n") 
            stop("Specific cluster ranger was skipped! Try expanding the cluster range or reducing the resolution step size.") 
        } 
 
        if(verbose) 
            cat(nCluster, "clusters.\n") 
    } 
 
    obj = cluster_start_with_one(obj) 
    obj@misc[["best.resolution"]] = res 
    if(verbose) 
        cat("Final resolution:", res, "with", nCluster, "clusters.\n") 
    return(obj) 
} 
 
#' @title numToMax 
#' @author Jiahao Wang 
#' @export
numToMax <- function(x){ 
    x = max(x) 
    ceiling(x / 10 ^ floor(log10(x) - 1) + 4) * 10 ^ floor(log10(x) - 1) 
} 
 
#' @title RemoveDoublets 
#' @author Jiahao Wang 
#' @export
RemoveDoublets <- function( 
    obj, 
    doublet.rate, 
    sample_name, 
    pN = 0.25, 
    pc.num = 1:30, 
    use.SCT = FALSE, 
    remove = TRUE 
  ){ 
 
    # ref: 
    # - https://www.jianshu.com/p/6770c6a05287 
    # - https://www.jianshu.com/p/0cb401b1ebe6 
    # - https://zhuanlan.zhihu.com/p/469335625 
 
    ## 寻找最优pK值 
    sweep.res.list <- paramSweep_v3(obj, PCs = pc.num, sct = F) 
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)   
    bcmvn <- find.pK(sweep.stats) 
    pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() 
 
    ## 排除不能检出的同源doublets，优化期望的doublets数量 
    homotypic.prop <- modelHomotypic(obj$seurat_clusters)   # 最好提供celltype 
    nExp_poi <- round(doublet.rate*ncol(obj))  
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) 
    seu.scored <- doubletFinder_v3(obj, PCs = pc.num, pN = 0.25, pK = pK_bcmvn,  
                                   nExp = nExp_poi.adj, reuse.pANN = F, sct = F) 
    # pick out doublets 
    cname <-colnames(seu.scored[[]]) 
    DF <- cname[grep('^DF',cname)] 
    seu.scored[["doublet"]] <- as.numeric(seu.scored[[DF]]=="Doublet") 
     
    table(seu.scored@meta.data$doublet) 
    # remove doublets 
    if(remove == "TRUE"){ 
        seu.removed <- subset(seu.scored, subset = doublet != 1) 
    } else{ 
        seu.removed <- seu.scored 
    } 
 
    p1 <- DimPlot(seu.scored, group.by = DF) 
    res.list <- list("plot" = p1, "obj" = seu.removed) 
    return(res.list) 
} 
 
# subclusterMakers <- list( 
#     "T/I/NK cell"      = c("CD3D", "CD3", "CD8A", "CD8B", "CD4", "TRGC1", "TRGC2", "TRDC", "KRT86", "KIT",  
#                            "LST1", "FGFBP2", "TYROBP", "GNLY", "FCGR3A", "KLRF1", "MKI67", "STMN1"), 
#     "Epithelial cell"  = c(), 
#     "B cell"           = c('IL1B', 'PTGS2', 'VEGFA', 'FCN1', 'S100A8', 'S100A9', 'S100A12', 'CD1C',  
#                            'FCER1A', 'CLEC10A', 'LILRA4', 'CLIC3', 'CXCR3', 'CCL19', 'LAMP3', 'CCR7',  
#                            'CALHM6', 'VAMP5', 'LDHA', 'CHI3L1', 'SPP1', 'FBP1', 'MMP9', 'C1QC', 'FOLR2',  
#                            'SELENOP', 'C1QB', 'TMSB4X', 'YBX1', 'DNASE1L3', 'PCLAF', 'TYMS', 'UBE2C'), 
#     "Myeloid cell"           = c('CD19', 'MS4A1', 'CD79A', 'CD79B', 'FCER2', 'TCL1A', 'IGHD', 'SELL', 'NR4A1',  
#                            'NR4A2', 'CD69', 'CXCR4', 'TOMM7', 'LTB', 'CD48', 'BANK1', 'RGS13', 'LMO2',  
#                            'NEIL1', 'BCL6', 'XBP1', 'MZB1', 'IGHA1', 'IGHA2', 'IGHG1', 'IGHG2', 'STMN1',  
#                            'MKI67', 'UBE2C', 'HMGB2', 'PCLAF'), 
#     "Endothelial cell" = c('CXCL12', 'SULF1', 'SEMA3G', 'GJA5', 'IGFBP3', 'PLVAP', 'PLPP3', 'COL4A1',  
#                            'COL4A2', 'COL15A1', 'IL32', 'SYNE2', 'ARGLU1', 'WSB1', 'N4BP2L2', 'ACKR1',  
#                            'SOCS3', 'GADD45B', 'HLA−DRA', 'ZFP36', 'CCL21', 'LYVE1', 'PR', 'OX1', 'TFF3',  
#                            'AKAP12', 'CPE', 'MADCAM1', 'ADGRG6', 'CLDN5'), 
#     "Fibroblast"       = c('MYH11', 'ACTA2', 'MCAM', 'CAV1', 'CXCL12', 'CFD', 'MFAP5', 'IGFBP6', 'DCN',  
#                            'RSPO3', 'WNT2B', 'CCL2', 'POSTN', 'WNT5A', 'COL3A1', 'TMEM158', 'CTHRC1',  
#                            'COL6A3', 'FAP') 
#                         ) 
 
###### color for use 
mycolors2 = c('#E4E4E4','#FFF7EC','#FEE8C8','#FDD49E','#FDBB84','#FC8D59','#EF6548','#D7301F','#B30000','#7F0000') 
 
col.cluster2.1 <- c('#1C79B7','#F38329','#2DA248','#DC403E','#976BA6','#8D574C','#D07DB0','#BDBF3B', '#27BDD0','#B0C7E5','#F9BA79','#A0CC8A','#F09EA1','#C7B0D1','#D4BCB8','#F7CCDD', 
                   '#DBDD8D','#A3D7E5','#B04E4F','#A38A58','#ED5351','#0C8945','#0F71A7','#A82764', '#F8DAE4','#7A4B1F','#5E6BAE','#8AC997','#DAC9AC','#0F5148','#A0B5DC','#9F858E', 
                   '#5C181D','#7B8380','#E8DDDA','#264220','#5AB747','#5169AE','#4B3D58','#CD428A', '#62615A','#B82129','#66762E') 
col.cluster <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', 
                 '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', 
                 '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') 
 
llc_cols <- c('#b3d3e8','#7dcade','#368dab', '#259b39', '#65b72f', '#a2ce89', '#e4361e', '#93646d', '#c0a26e', '#edb000', '#eb6938', '#ed996d', '#ed7372', '#e969a2', '#d3b7d0', '#825796', '#cd0081', '#e3e39a', '#7491ca', '#b85919') 
 
#### yellow-blue color 
yellow_blue_color <- c('#450458', '#48196B', '#482878', '#472F7D', '#443A83', '#3E4C8A', '#375A8C', '#34608D', '#2C718E','#287C8E', '#20928C', '#1E9A8A', '#28AE80', '#3FBC72', '#60CA60', '#90D743', '#B2DD2D', '#CDE11D', '#E7E419','#EFE51C') 
 
stage_cols <- c('#F3BA2F', '#E5916F', '#068085', '#A72F5C', '#FBE4B8', '#E42988', '#D8ACCF', '#21863A', '#F2CF1E', '#DADADA', '#1591C0', '#844482', '#293464', '#F4B2BC', '#AFDCE0', '#897664', '#E5481A', '#9FD08C', '#748595', '#EAE443') 
 
#' @title creatSeuratObj 
#' @author Jiahao Wang 
#' @export
creatSeuratObj  <- function(countFile, metaFile){ 
 
    counts = readCounts(countFile) 
    meta = readMeta(metaFile) 
    obnj = CreateSeuratObject(counts = counts, meta.data = meta) 
} 
 
 
#' @title FeaturePlot2 
#' @author Jiahao Wang 
#' @export
FeaturePlot2  <- function(obj, features, axis = TRUE, miniAxis = FALSE, legend = TRUE){ 
    loadp(ggplot2) 
    umapM = Embeddings(obj, reduction = "umap") 
    plotL = list() 
    umapM = Embeddings(obj, reduction = "umap") 
    pdata = data.frame(UMAP_1 = umapM[, 1], UMAP_2 = umapM[, 2]) 
    for(feature in features){ 
        pdata$Expression = log2(obj@assays$RNA@counts[feature ,] + 1) 
        p <- ggplot(pdata, aes(UMAP_1, UMAP_2)) 
        p <- p + geom_point(aes(color = Expression), size = 0.3) + setText(20) 
        p <- p + scale_color_continuous(low = "lightgrey", high = "#DE1F1F") 
        p <- p + labs(title = feature) 
        p <- p + guides(color = guide_colorbar(barwidth = unit(0.6, "cm"), barheight = unit(3, "cm"))) 
 
        if(!axis){ 
            p <- p + theme_blank() 
        } 
 
        if(!legend){ 
            p <- p + theme(legend.position = "none") 
        } 
 
        if(miniAxis){ 
            x.vector = pdata[, 1] 
            y.vector = pdata[, 2] 
            segTable = data.frame(x = c(min(x.vector) - 1.5 - 0.03, min(x.vector) - 1.5), 
                                  y = c(min(y.vector) - 1.5, min(y.vector) - 1.5), 
                                  vx = c(8, 0), 
                                  vy = c(0, 8)) 
            p <- p + geom_segment(data = segTable, mapping = aes(x = x, y = y, xend = x + vx, yend = y + vy), 
                                                                 arrow = arrow(length=unit(0.25, "cm")),  
                                                                 size = 0.5) 
            p <- p + annotate("text", x = min(x.vector) + 1.1, y = min(y.vector) - 2.3, label = "UMAP_1", size = 4, hjust = 0, parse = TRUE) 
            p <- p + annotate("text", x = min(x.vector) - 2.5, y = min(y.vector) + 1.4, label = "UMAP_2", size = 4, hjust = 0, parse = TRUE, angle = 90) 
        } 
        plotL[[feature]] <- p 
    } 
 
    loadp(gridExtra) 
    res <- grid.arrange(grobs = plotL, nrow = 1)  
    return(res) 
} 
 
#' @title getMultipletRate 
#' @author Jiahao Wang 
#' @export
getMultipletRate  <- function(obj){ 
    # https://kb.10xgenomics.com/hc/en-us/articles/360059124751-Why-is-the-multiplet-rate-different-for-the-Next-GEM-Single-Cell-3-LT-v3-1-assay-compared-to-other-single-cell-applications- 
    library_size = seq(1:10) * 1000 
    multiplet_rate = c(0.8, 1.5, 2.3, 3, 3.8, 4.6, 5.3, 6.1, 6.8, 8.0) / 100 
    multiplet_rate_table = data.frame(library_size,multiplet_rate) 
 
    library_size_obj = round(dim(obj)[2] / 1000) 
    library_size_obj[library_size_obj > 10] = 10 
    multiplet_rate_obj = multiplet_rate_table[library_size_obj, 2] 
    return(multiplet_rate_obj) 
} 
 
#' @title mergeSingleCellranger 
#' @author Jiahao Wang 
#' @export
mergeSingleCellranger <- function(dir, min.cells = 1, min.features = 1, method = "count"){ 
    sample = lf("x", dir) 
    dirEach = switch(method, 
                count = paste0(dir, "/", sample, "/outs/filtered_feature_bc_matrix"), 
                multi = paste0(dir, "/", sample, "/outs/per_sample_outs/", sample, "/count/sample_filtered_feature_bc_matrix") 
                ) 
    # check 
    notExist = dirEach[!file.exists(paste0(dirEach, "/matrix.mtx.gz"))] 
    if(length(notExist) > 0){ 
        message("Directory below does not exist or incomplete:") 
        invisible(sapply(notExist, function(x) message("\t", x))) 
        stop() 
    }     
    names(dirEach) = gsub("_", "-", sample) 
    counts = Read10X(data.dir = dirEach) 
    obj = CreateSeuratObject(counts, min.cells = min.cells, min.features = min.features) 
    return(obj) 
} 
 
#' @title mkdir 
#' @author Jiahao Wang 
#' @export
mkdir <- function(...){ 
    folders = as.character(list(...)) 
    for(folder in folders){ 
        if(!dir.exists(folder)) 
            dir.create(folder, recursive = TRUE) 
    } 
} 
 
#' @title multiPara 
#' @author Jiahao Wang 
#' @export
multiPara <- function(func, # 要执行的函数  
                    ..., # func的动态参数 
                    paraL = NULL, # func的静态参数 
                    bind = "c", # 结果的合并方式 
                    cores = NULL # 使用的核数 
                    ){ 
 
    suppressPackageStartupMessages(library(doParallel)) 
    if(is.null(cores)){ 
        cores <- detectCores(logical = FALSE) 
    } 
 
    cl <- makeCluster(cores) 
    registerDoParallel(cl) 
 
    argL <- list(...) 
    result <- suppressWarnings(foreach(i = seq(length(argL[[1]])), .combine = bind) %dopar%  
        do.call(func, c(lapply(argL, `[`, i), paraL))) 
    stopCluster(cl) 
 
    invisible(result) 
} 
 
 
 
#' @title processSigMarker 
#' @author Jiahao Wang 
#' @export
processSigMarker <- function(object, i, n = 50, write = F) { 
                                                            df <- FindMarkers(object, ident.1 = i, only.pos = T) 
                                                            df.sig <- subset(df, df$p_val < 0.05) 
                                                            df.sig <- df.sig[with(df.sig, order(p_val)),] 
                                                            print(dim(df.sig)) 
                                     
                                                            if (write == T) 
                                                            { 
                                                                out <- cbind(gene_name = rownames(df.sig), df.sig) 
                                                                filname = paste0(as.character(i), '_',as.character(length(out[,1])), 'sig_markers.txt') 
                                                                write.table(out, file = filname, row.names = F, col.names = T, quote = F, sep='\t') 
                                                            } 
 
                                                            topsig <- rownames(df.sig[with(df.sig, order(p_val)),][1:n,]) 
                                                            return(topsig) 
                                                            } 
 
#' @title QC_MT 
#' @author Jiahao Wang 
#' @export
QC_MT <- function(object, n = 30) 
{ 
    object <- PercentageFeatureSet(object, pattern = "^MT-", col.name = "percent.mt") 
    temp <- table(cut(object$percent.mt, breaks = seq(0, 100, 5))) 
    temp <- as.data.frame(temp) 
    print(ggplot(temp, aes(x = Var1, y = Freq)) + geom_bar(position = "stack", stat = "identity", width = 0.6) + gg_qc) 
 
    print(n) 
    print(dim(object)[2]) 
    print(dim(subset(object, percent.mt < n))[2]) 
    print(dim(subset(object, percent.mt < n))[2]/dim(object)[2]) 
 
    return(object) 
} 
 
#' @title QC_nFnC 
#' @author Jiahao Wang 
#' @export
QC_nFnC <- function(object, n = 0.8) 
{ 
    object$nFeature.nCount <- object$nFeature_RNA/object$nCount_RNA 
    temp <- table(cut(object$nFeature.nCount, breaks = seq(0,1,0.05))) 
    temp <- as.data.frame(temp) 
    print(ggplot(temp, aes(x = Var1, y = Freq)) + geom_bar(position = "stack", stat = "identity", width = 0.6) + gg_qc) 
 
    print(n) 
    print(dim(object)[2]) 
    print(dim(subset(object, nFeature.nCount < n))[2]) 
    print(dim(subset(object, nFeature.nCount < n))[2]/dim(object)[2]) 
     
    return(object) 
} 
 
#' @title QC_Scat 
#' @author Jiahao Wang 
#' @export
QC_Scat <- function(object, hyl =200, hyh = 7000, vx = 4e+4) 
{ 
    print(FeatureScatter(object, feature1 = 'nCount_RNA', feature2 = "nFeature_RNA", pt.size = 0.01) + gg_style + theme_bw() + 
          geom_hline(yintercept = hyl, color = "darkred", linetype = "dashed") +  
          geom_hline(yintercept = hyh, color = "darkred", linetype = "dashed") +  
          geom_vline(xintercept = vx, color = "darkblue", linetype = "dashed")) 
 
    print(hyl);print(hyh) 
    print(dim(object)[2]) 
    print(dim(subset(object, nFeature_RNA > hyl & nFeature_RNA < hyh))[2]) 
    cat(dim(subset(object, nFeature_RNA > hyl & nFeature_RNA < hyh))[2]/dim(object)[2]) 
 
    ### nCount hist  
    temp <- table(cut(object$nCount_RNA, breaks = seq(0, 2e+5, 1e+4))) 
          temp <- as.data.frame(temp) 
    print(ggplot(temp, aes(x = Var1, y = Freq)) + geom_bar(position = "stack", stat = "identity", width = 0.6) + gg_qc) 
    print(vx) 
    print(dim(subset(object, nCount_RNA < vx))[2]) 
    print(dim(subset(object, nCount_RNA < vx))[2]/dim(object)[2]) 
     
    return(object) 
} 
 
#' @title QC_workflow 
#' @author Jiahao Wang 
#' @export
QC_workflow <- function(object, mt_n = 30, nFnC_n = 0.8, hyl =200, hyh = 7000, vx = 4e+4, method = 'check') 
{ 
    if (method == 'check') 
    { 
        object <- QC_MT(object, mt_n) 
        object <- QC_nFnC(object, nFnC_n) 
        object <- QC_Scat(object, hyl, hyh, vx) 
 
        return(object) 
    } 
     
    if (method == 'final') 
    { 
        object <- subset(object, percent.mt<mt_n & nFeature_RNA > hyl & nFeature_RNA < hyh & nCount_RNA < vx & nFeature.nCount < nFnC_n) 
        print(object) 
        return(object)   
    } 
} 
 
#' @title readCounts 
#' @author Jiahao Wang 
#' @export
readCounts <- function(file){ 
 
    library(data.table) 
    counts = fread(file, check.names = FALSE) 
    row.names = as.character(data.frame(counts[, 1])[, 1]) 
    counts[, 1] = NULL 
    counts = data.frame(counts, check.names = FALSE) 
    rownames(counts) = row.names 
    return(counts) 
} 
 
#' @title readfile2seurat 
#' @author Jiahao Wang 
#' @export
readfile2seurat <- function(file_path) # by fox 
 
{ 
    loadp(data.table, Seurat) 
    data <- fread(file_path) 
    rownames(data) <- data$V1 
    data$V1 <- NULL 
    data <- CreateSeuratObject(counts = data, min.cell = 0) 
    return(data) 
} 
 
#' @title selectMarkers 
#' @author Jiahao Wang 
#' @export
selectMarkers <- function(celltype = NULL){ 
    markers <- list( 
        B              = c("MS4A1", "MZB1", "CD79A", "SDC1"), 
        T              = c("CD3D", "CD3E"), 
        Endo           = c("VWF", "PECAM1", "ACKR1", "CLDN5", "ENG"), 
        Epis           = c("EPCAM", "CD24", "CEACAM5", "KRT8"), 
        Fibro          = c("COL3A1", "DCN", "ACTA2", "MYH11"), 
        Mast           = c("TPSAB1", "CPA3", "KIT"), # 归入Myeloid 
        Myeloid        = c("FCGR3A", "CD14", "CD68", "S100A8", "FCN1", "TPSAB1", "CPA3", "FLT3", "LAMP3"), 
        Neutrophils    = c("FCGR3A", "FCGR3B", "CST3", "LYZ", "CSF3R"), 
        DC             = c("FLT3", "LILRA4", "XCR1", "CLEC9A", "CD1C", "CD1E", "LAMP3", "FCER1A"), # 归入Myeloid 
        Smooth_muscle  = c("ACTA2", "TAGLN", "MUSTN1", "MYH11"), 
        T_I_NK         = c("CD3D", "CD3E", "CD3G", "CD2"), 
        MAIT           = c("SLC4A10", "TRAV1-2", "KLRB1"), 
        gdT            = c("TRDV1", "TRDV2", "TRDC", "KIR2DL4"), 
        ILC            = c("KRT86", "KIT", "LST1"), 
        CD8            = c("CD8A", "CD8B"), 
        CD4            = c("CD4"), 
        NK             = c("TYROBP", "GNLY", "FCGR3A", "NCAM1"), 
        Prolif         = c("MKI67", "STMN1", "UBE2C") 
        ) 
 
    if(is.null(celltype)){ 
        cat("Availabled celltype: \n\t") 
        cat(paste(names(markers)), "\n") 
    } else { 
        if(! celltype %in% names(markers)){ 
            cat("Not availabled celltype! Only support below: \n\t") 
            cat(paste(names(markers)), "\n") 
        } else {         
            return(markers[[celltype]]) 
        } 
    } 
} 
 
#' @title startSC 
#' @author Jiahao Wang 
#' @export
startSC <- function(lib = FALSE){ 
 
    if(lib) loadp(Seurat, ggplot2, dbplyr, RColorBrewer) 
    print(paste('Start', as.character(date()), sep = '    ')) 
} 
 
 
#' @title test_workflow 
#' @author Jiahao Wang 
#' @export
test_workflow <- function(object = NULL, nPC.method = 'sigPC', N = 40, markers.1st2check = general_marker, assay = 'RNA') 
                            { 
                                #### help 
                                if(class(object) != 'Seurat') 
                                { 
                                    cat("  Help for Seurat Workflow 1\n 
                                           test_workflow(object, nPC.method = 'sigPC', N = 40, markers.1st2check = general_marker, assay = 'RNA')\n  
                                           nPC.method:\n  
                                           sigPC: automaticaly choosing significant PCs by ScoreJackStraw\n 
                                           topPC: mannually choosing Top N PCs\n 
                                           markers.1st2check: general_marker by defalt 
                                           \n 
                                           example: lst <- test_workflow(seurat_object)\n 
                                                    seurat_object <- lst[[1]]\n 
                                                    lst\n 
                                           example2:lst <- test_workflow(seurat_object, 'topPC', 30)\n 
                                                    seurat_object <- lst[[1]]\n 
                                                    lst\n\n") 
 
                                }  
                                #### Main 
                                if(class(object) == 'Seurat') 
                                { 
                                    print(paste0('Started to process seurat testing workflow with ', nPC.method, ' method'))     
 
                                    object <- NormalizeData(object) 
                                    object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, verbose = FALSE) 
                                    object <- ScaleData(object, features = rownames(object), verbose = F) 
                                    object <- RunPCA(object) 
 
                                    #### PCA Plots 
                                    DimHeatmap(object, dims = 1:15, cells = 100, balanced = TRUE) 
                                    pe <- ElbowPlot(object, ndims = 50)     
 
                                    #### significant PCs selected by jackstraw score : nPC      
                                    object <- JackStraw(object, dims = 50) 
                                    object <- ScoreJackStraw(object, dims = 1:50) 
                                    pj <- JackStrawPlot(object, dims = 1:50) 
                                    JackResult <- object@reductions$pca@jackstraw@overall.p.values 
 
                                    #### nPC.method 
                                    if (nPC.method == 'sigPC') 
                                    { 
                                        nPC <- JackResult[,'PC'][which(JackResult[,'Score'] < 0.05)] 
                                    } else if (nPC.method == 'topPC') 
                                    { 
                                        nPC <- N 
                                    } 
 
                                    #### Run TSNE UMAP 
                                    object <- RunTSNE(object, reduction = "pca", dims = nPC, seed.use = 7, verbose = F) 
                                    object <- RunUMAP(object, reduction = "pca", dims = nPC, verbose = F) 
                                     
                                    #### orig.ident 
                                    pumap <- DimPlot(object, group.by = 'orig.ident', cols = 'grey', pt.size = 0.1) + gg_style 
                                    ptsne <- DimPlot(object, reduction = 'tsne', group.by = 'orig.ident', cols = 'grey', pt.size=0.1) + gg_style 
 
                                    #### Run SNN 
                                    object <- FindNeighbors(object, reduction = "pca", dims = nPC, verbose = F) 
                                    object <- FindClusters(object, resolution = seq(0.1, 1.2, by=0.1), verbose = F) 
                                     
                                    if (assay == 'RNA') 
                                    { 
                                        p0.7umap <- DimPlot(object, group.by = 'RNA_snn_res.0.7', label=T, pt.size=0.1) + gg_style 
                                        p0.7tsne <- DimPlot(object, reduction = 'tsne', group.by = 'RNA_snn_res.0.7', label=T, pt.size=0.01) + gg_style 
                                    } 
 
                                    pf <- FeaturePlot(object, reduction = 'tsne', features = markers.1st2check,  pt.size = 0.1, cols = mycolors2) 
                                     
                                    ### Return Out list 
                                    out_list <- list('object' = object, 'elbow.plot' = pe, 'jack.plot' = pj, 'tsne.plot'= ptsne, 'umap.plot' = pumap,  
                                                     'umap0.7' = p0.7umap, 'tsne0.7' = p0.7tsne, 'general.feature.plot' = pf, 'nPC' = nPC) 
                                     
                                    print('seurat testing workflow completed') 
                                     
                                    return(out_list) 
                                } 
                            } 
