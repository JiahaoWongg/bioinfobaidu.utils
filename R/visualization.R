#' @title show_my_pet 
#' @author Jiahao Wang 
#' @export
show_my_pet <- function(pet_name = "lying_cat", list_my_pet = FALSE){ 
    image_dir = "/public1/home/scg5695/code/myscript/yyds/ASCII_image" 
	if(list_my_pet){ 
		my_pets = gsub(".txt", "", lf("x", image_dir)) 
		cat(paste(my_pets, collapse = "\n"), "\n") 
	} else { 
	    file_path = paste0(image_dir, "/", pet_name, ".txt") 
	    if(!file.exists(file_path)){ 
	        cat(paste0("Can't find your pet named '", pet_name, "', maybe lost?\n")) 
	    } else { 
			cat(paste0("\n", paste(readLines(file_path), collapse = "\n"), "\n\n")) 
	    } 
	} 
} 
 
#' @title myDonkey 
#' @author Jiahao Wang 
#' @export
myDonkey <- function(icansay = NULL){ 
	cat("\n")  
	cat( "        \\  ^__^     ", paste0(icansay, "~haw~"), "           \n") 
	cat( "         \\ (oo)\\________         \n") 
	cat( "           (__)\\        )\\/\\    \n") 
	cat( "               ||-----w |          \n") 
	cat( "               ||      ||          \n\n") 
} 
 
Blue = "#27A2EF" 
Red = "#F93463" 
 
#' @title barplot2 
#' @author Jiahao Wang 
#' @export
barplot2  <- function(data, 
				  	  outName = NULL,  
					  PDF = FALSE,  
					  PNG = TRUE,  
					  w = 6,  
					  h = 6, 
					  n = 20, 
					  title = "",  
					  padjust = TRUE,  
					  addLine = FALSE,  
					  fill = "#0089CB" 
					  ){ 
 
	if(is.null(out)) 
		stop("Please specific parameter 'out' to save pdf file.\n") 
 
	if(n > nrow(data)) 
		n = nrow(data) 
	data = data.frame(data)[1:n, ] 
 
	if(padjust){ 
		data = data[order(data$p.adjust), ] 
		data$Description = factor(data$Description, levels = rev(data$Description)) 
		p <- ggplot(data) + setTheme() + setText(50) 
		p <- p + geom_bar(aes(x = Description, y = -log2(p.adjust)), stat = 'identity', width = 0.8, fill = fill, alpha = 1) 
	} else { 
		data = data[order(data$pvalue), ] 
		data$Description = factor(data$Description, levels = rev(data$Description)) 
		p <- ggplot(data) + setTheme() + setText(50) 
		p <- p + geom_bar(aes(x = Description, y = -log2(pvalue)), stat = 'identity', width = 0.8, fill = fill, alpha = 1) 
 
	} 
 
	p <- p + labs(x = "", title = title) + coord_flip() + theme(plot.title = element_text(hjust = 0)) 
	p <- p + geom_hline(aes(yintercept = -log2(0.05)), colour = "black", linetype = "dashed", size = 1.5) 
	p <- p + theme(axis.text.y = element_text(hjust = 1)) + theme(axis.text.y  = element_text(color = "black")) 
 
	if(addCount){ 
		p <- p + geom_line(aes(x = Description, y = Count, group = 1), size = 2, color = "#BC0909") 
		p <- p + geom_point(aes(x = Description, y = Count), size = 5) 
	} 
	save.graph(p, file = outName, w, h, PDF = PDF, PNG = PNG) 
} 
 
#' @title barplotGO 
#' @author Jiahao Wang 
#' @export
barplotGO  <- function(data, size = 30){ 
 
	data = data[order(data$p.adjust), ] 
	data$Description = factor(data$Description, levels = rev(data$Description)) 
	data$GeneRatioValue = as.numeric(subString(data$GeneRatio, 1, "/")) / as.numeric(subString(data$GeneRatio, 2, "/")) 
	data = data[order(data$GeneRatioValue, decreasing = TRUE), ] 
 
	p <- ggplot(data) + setTheme() + setText(size) 
	p <- p + geom_bar(aes(x = Description, y = -log10(p.adjust)), stat = 'identity', width = 0.8, fill = "#DD1C77", alpha = 1) 
	p <- p + scale_colour_gradient(low = "#DD1C77", high = "#3182bd", limit = c(min(data$p.adjust), max(data$p.adjust)), guide = guide_colourbar(reverse = TRUE)) 
	p <- p + labs(x = "GeneRatio", y = "", title = "") 
 
	p <- p + labs(x = "", title = "title") + coord_flip() + theme(plot.title = element_text(hjust = 0)) 
	p <- p + geom_hline(aes(yintercept = -log10(0.05)), colour = "black", linetype = "dashed", size = 1.5) 
	p <- p + theme(axis.text.y = element_text(hjust = 1)) + theme(axis.text.y  = element_text(color = "black")) 
 
	p <- p + facet_grid(ONTOLOGY ~ ., scale="free") 
	p <- p + theme(strip.text.y = element_text(size = 15, angle = 0)) 
	p <- p + theme(strip.background.y = element_rect(fill = "white")) 
	p <- p + theme(legend.key.size = unit(40, "pt")) 
	p <- p + theme(legend.title=element_text(size=unit(25, "pt"))) 
	p <- p + theme(legend.text=element_text(size=unit(20, "pt"))) 
	return(p) 
} 
 
#' @title bedGraphToGR 
#' @author Jiahao Wang 
#' @export
bedGraphToGR  <- function(Tx){ 
 
	chr = as.character(Tx[,1]) 
	ST = as.numeric(Tx[,2]) + 1 
	ED = as.numeric(Tx[,3]) 
	methyCount = as.numeric(Tx[,5]) 
	readCount = methyCount + as.numeric(Tx[,6]) 
	beta = signif(methyCount / readCount, 4) 
	gr = GRanges(seqnames = Rle(chr),  ranges = IRanges(ST,ED),  strand = "*", methyCount, readCount, beta) 
	return(gr) 
} 
 
 
#' @title dotplot2 
#' @author Jiahao Wang 
#' @export
dotplot2  <- function(data,  
					 w 		   = 6,  
					 h 		   = 6,  
					 n 		   = 20,  
					 colorHigh = "#DD1C77",  
					 colorLow  = "#3182bd", 
					 sortBy    = "p.adjust",  
					 title     = "",  
					 padjust   = TRUE, 
					 outName 	   = NULL, 
					 PDF = FALSE, 
					 PNG = TRUE 
					 ){ 
 
	if(is.null(outName)) 
		stop("Please specific parameter 'outName' to save pdf file.\n") 
 
	if(n > nrow(data)) 
		n = nrow(data) 
	data = data.frame(data)[1:n, ] 
 
	data = data[order(data$Count, decreasing = TRUE), ] 
	if(sortBy == "Count"){ 
		data$Description = factor(data$Description, levels = rev(data$Description)) 
	}else if(sortBy == "p.adjust"){ 
		if(padjust){ 
			data = data[order(data$p.adjust), ] 
			data$Description = factor(data$Description, levels = rev(data$Description)) 
		}else{ 
			data = data[order(data$pvalue), ] 
			data$Description = factor(data$Description, levels = rev(data$Description)) 
		}			 
	} 
 
	data$GeneRatioValue = as.numeric(subString(data$GeneRatio, 1, "/")) / as.numeric(subString(data$GeneRatio, 2, "/")) 
	data = data[order(data$GeneRatioValue, decreasing = TRUE), ] 
 
	p <- ggplot(data, aes(x = GeneRatioValue, y = Description)) + theme_bw() 
	if(padjust){ 
		p <- p + geom_point(aes(color = p.adjust, size = Count)) 
		p <- p + scale_colour_gradient(low = colorHigh, high = colorLow, limit = c(min(data$p.adjust), max(data$p.adjust)), guide = guide_colourbar(reverse = TRUE)) 
	} else{ 
		p <- p + geom_point(aes(color = pvalue, size = Count))		 
		p <- p + scale_colour_gradient(low = colorHigh, high = colorLow, limit = c(min(data$pvalue), max(data$pvalue)), guide = guide_colourbar(reverse = TRUE)) 
	} 
 
	p <- p + labs(x = "GeneRatio", y = "", title = title) + theme(plot.title = element_text(hjust = 0.5)) 
	p <- p + theme(axis.text.y  = element_text(color = "black")) 
 
	if(!is.null(outName)){ 
		if(is.null(w)){ 
		    max_nchar = max(nchar(as.character(data$Description))) 
		    w =  max_nchar * 0.13 + 2.2 
		} 
		save.graph(p, file = outName, w, h, PDF = PDF, PNG = PNG) 
	} 
	return(p) 
} 
 
#' @title dotplot3 
#' @author Jiahao Wang 
#' @export
dotplot3  <- function(egoT, top = 20, padjust = TRUE, title = "GO enrichment"){ 
 
	library("gridExtra") 
 
	data = egoT 
	onts = unique(data$ONTOLOGY) 
 
	plots = list() 
	for(i in 1:length(onts)){ 
 
		ont = onts[i] 
		xdata = data[data$ONTOLOGY == ont, ] 
		xdata = xdata[order(xdata$Count, decreasing = TRUE), ] 
		xdata$Description = factor(xdata$Description, levels = rev(xdata$Description)) 
		xdata$GeneRatioValue = as.numeric(subString(xdata$GeneRatio, 1, "/")) / as.numeric(subString(xdata$GeneRatio, 2, "/")) 
 
		p <- ggplot(xdata, aes(x = GeneRatioValue, y = Description)) + theme_bw() 
		if(padjust){ 
			p <- p + geom_point(aes(color = p.adjust, size = Count)) + setText(16) 
			p <- p + scale_colour_gradient(low = "red", high = "blue", limit = c(min(data$p.adjust), max(data$p.adjust)), guide = guide_colourbar(reverse = TRUE)) 
		} else{ 
			p <- p + geom_point(aes(color = pvalue, size = Count))		 
			p <- p + scale_colour_gradient(low = "red", high = "blue", limit = c(min(data$pvalue), max(data$pvalue)), guide = guide_colourbar(reverse = TRUE)) 
		} 
 
		p <- p + labs(x = "GeneRatio", y = "", title = paste("GO:", ont)) 
		p <- p + theme(axis.text.y  = element_text(color = "black")) 
		p <- p + theme(axis.text.y = element_text(hjust = 1)) 
		p <- p + theme(plot.title = element_text(hjust = 0.5)) 
		plots[[i]] = p 
	} 
 
	p_res <- grid.arrange(grobs = plots, ncol = 3)  
	return(p_res) 
 
} 
 
#' @title dotplotGO 
#' @author Jiahao Wang 
#' @export
dotplotGO  <- function(data, size = 28, padj = TRUE,  
					   outName = NULL, w = NULL, h = 8, 
					   PDF = FALSE, PNG = TRUE){ 
   
    if(padj){ 
        data = data[data$p.adjust < 0.05, ] 
        data = data[order(data$p.adjust), ] 
    } else { 
        data = data[data$pvalue < 0.05, ] 
        data = data[order(data$pvalue), ]         
    } 
 
    # only show top 10 pathways 
    data_BP = data[data$ONTOLOGY == "BP", ] 
    data_BP = data_BP[1:ifelse(nrow(data_BP) < 10, nrow(data_BP), 10), ] 
 
    data_CC = data[data$ONTOLOGY == "CC", ] 
    data_CC = data_CC[1:ifelse(nrow(data_CC) < 10, nrow(data_CC), 10), ] 
 
    data_MF = data[data$ONTOLOGY == "MF", ] 
    data_MF = data_MF[1:ifelse(nrow(data_MF) < 10, nrow(data_MF), 10), ] 
    data = rbind(data_BP, data_CC, data_MF) 
 
    data$Description = factor(data$Description, levels = rev(data$Description)) 
    data$GeneRatioValue = as.numeric(subString(data$GeneRatio, 1, "/")) / as.numeric(subString(data$GeneRatio, 2, "/")) 
    data = data[order(data$GeneRatioValue, decreasing = TRUE), ] 
 
    # p <- ggplot(data, aes(x = 1, y = Description)) + theme_bw() + setTheme() + setText(25)\ 
    loadp(ggplot2, ggthemes) 
    p <- ggplot(data, aes(x = GeneRatioValue, y = Description)) + setText(25, graph.theme = "bw") 
    if(padj){ 
        p <- p + geom_point(aes(color = p.adjust, size = Count)) 
        p <- p + scale_colour_gradient(low = "#DD1C77", high = "#3182bd", limit = c(min(data$p.adjust), max(data$p.adjust)), guide = guide_colourbar(reverse = TRUE)) 
    } else { 
        p <- p + geom_point(aes(color = pvalue, size = Count)) 
        p <- p + scale_colour_gradient(low = "#DD1C77", high = "#3182bd", limit = c(min(data$pvalue), max(data$pvalue)), guide = guide_colourbar(reverse = TRUE)) 
    } 
    p <- p + labs(x = "Gene Ratio", y = "", title = "") 
    p <- p + theme(axis.text.y  = element_text(color = "black")) + theme(axis.text.y = element_text(hjust = 1)) 
    p <- p + theme(axis.text.x  = element_text(angle = 45, hjust = 1)) 
    # p <- p + theme(axis.ticks.x = element_blank()) 
    p <- p + facet_grid(ONTOLOGY ~ ., scale="free") 
    p <- p + theme(strip.text.y = element_text(size = size * 0.5, angle = 0)) 
    p <- p + theme(strip.background.y = element_rect(fill = "white")) 
 
    if(!is.null(outName)){ 
        if(is.null(w)){ 
            max_nchar = max(nchar(as.character(data$Description))) 
            w =  max_nchar * 0.13 + 2.2 
        } 
        save.graph(p, file = outName, w, h, PDF = PDF, PNG = PNG) 
    } 
    return(p) 
} 
 
#' @title dotplotKEGG 
#' @author Jiahao Wang 
#' @export
dotplotKEGG  <- function(data, size = 28, padj = TRUE, out = NULL, w = NULL, h = 8){ 
   
    if(padj){ 
        data = data[data$p.adjust < 0.05, ] 
        data = data[order(data$p.adjust), ] 
    } else { 
        data = data[data$pvalue < 0.05, ] 
        data = data[order(data$pvalue), ]         
    } 
 
    # only show top 30 pathways 
    data = na.omit(data[1:30, ]) 
 
    data$Description = factor(data$Description, levels = rev(data$Description)) 
    data$GeneRatioValue = as.numeric(subString(data$GeneRatio, 1, "/")) / as.numeric(subString(data$GeneRatio, 2, "/")) 
    data = data[order(data$GeneRatioValue, decreasing = TRUE), ] 
 
    # p <- ggplot(data, aes(x = 1, y = Description)) + theme_bw() + setTheme() + setText(25)\ 
    loadp(ggthemes) 
    p <- ggplot(data, aes(x = GeneRatioValue, y = Description)) + setText(25, graph.theme = "bw") 
    if(padj){ 
        p <- p + geom_point(aes(color = p.adjust, size = Count)) 
        p <- p + scale_colour_gradient(low = "#DD1C77", high = "#3182bd", limit = c(min(data$p.adjust), max(data$p.adjust)), guide = guide_colourbar(reverse = TRUE)) 
    } else { 
        p <- p + geom_point(aes(color = pvalue, size = Count)) 
        p <- p + scale_colour_gradient(low = "#DD1C77", high = "#3182bd", limit = c(min(data$pvalue), max(data$pvalue)), guide = guide_colourbar(reverse = TRUE)) 
    } 
    p <- p + labs(x = "Gene Ratio", y = "", title = "") 
    p <- p + theme(axis.text.y  = element_text(color = "black")) + theme(axis.text.y = element_text(hjust = 1)) 
    p <- p + theme(axis.text.x  = element_text(angle = 45, hjust = 1)) 
    # p <- p + theme(axis.ticks.x = element_blank()) 
 
    if(!is.null(out)){ 
        if(is.null(w)){ 
            max_nchar = max(nchar(as.character(data$Description))) 
            w =  max_nchar * 0.1109 + 1.5147 
        }             
        pdf(out, w, h) 
            print(p) 
        dev.off() 
        } 
    return(p) 
} 
 
#' @title gghelp 
#' @author Jiahao Wang 
#' @export
gghelp <- function(){ 
 
	cat('plot.title   = element_text(color = "black", size = title.size, face = face, hjust = 0.5)\n') 
	cat('plot.title   = element_text(color = "black", size = title.size, face = face, hjust = 0.5)\n') 
	cat('plot.title   = element_text(color = "black", size = title.size, face = face, hjust = 0.5)\n') 
	cat('plot.title   = element_text(color = "black", size = title.size, face = face, hjust = 0.5)\n') 
	cat('plot.title   = element_text(color = "black", size = title.size, face = face, hjust = 0.5)\n') 
} 
 
#' @title pitchPlot 
#' @author Jiahao Wang 
#' @export
pitchPlot <- function(){ 
 
	if(FALSE){ 
		library("gridExtra") 
		pdf(paste0("PDF/TCGA_mutate_driver_", treat, ".pdf"), 15, 9) 
			grid.arrange(grobs = plots, ncol = 8)  
		dev.off() 
	} 
} 
 
#' @title plotEGO 
#' @author Jiahao Wang 
#' @export
plotEGO <- function(ego, title = ""){ 
	p <- dotplot(ego, showCategory = nrow(ego)) + labs(title = title) 
	p <- p + theme(plot.title = element_text(size = 20, face = "bold")) 
	return(p) 
} 
 
#' @title plotEmpty 
#' @author Jiahao Wang 
#' @export
plotEmpty <- function(){ 
	p <- ggplot() + geom_point(aes(1, 1), colour = "white") +    
		 theme(axis.ticks = element_blank(),  
        	   panel.background = element_blank(),  
        	   axis.line = element_blank(),  
        	   axis.text.x = element_blank(),  
        	   axis.text.y = element_blank(),  
        	   axis.title.x = element_blank(),  
        	   axis.title.y = element_blank()) 
	return(p) 
} 
 
#' @title plotVolcano 
#' @author Jiahao Wang 
#' @export
plotVolcano <- function(diffM,  
						DIFFcutoff = 2,  
						FDRcutoff = 0.05,  
						title = "DGE_volcano",  
						xtitle = "log2(FoldChange)",  
						xlim = 5,  
						ylim = 10, 
						padjust = FALSE, 
						size = 0.1, 
						pos = 8.5, 
						label = TRUE 
						){ 
 
	library("ggplot2") 
	if(padjust){ 
		p <- ggplot(diffM, aes(x = log2FC, y = -log10(padj), color = DGE)) + xlim(-xlim, xlim) + ylim(0, ylim) 
		p <- p + theme_bw() + labs(title = title, x = xtitle, y = "-log10(FDR)") 
	}else{ 
		p <- ggplot(diffM, aes(x = log2FC, y = -log10(pvalue), color = DGE)) + xlim(-xlim, xlim) + ylim(0, ylim) 
		p <- p + theme_bw() + labs(title = title, x = xtitle, y = "-log10(P-value)") 
	} 
 
	p <- p + geom_point(size = size) + setText(20) + theme(plot.title = element_text(hjust = 0.5)) 
	p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
	p <- p + scale_color_manual(values = c("UP" = "red", "DN" = "dodgerblue4", "NC" = "grey")) 
	p <- p + theme(legend.position = "none") 
 
	p <- p + geom_vline(xintercept = c(log2(DIFFcutoff), -log2(DIFFcutoff)), lty = 4, col = "black", lwd = 0.6) 
	p <- p + geom_hline(yintercept = -log10(FDRcutoff),lty = 4, col= "black", lwd= 0.6) 
 
	N_DN = sum(DEGtable$DGE == "DN") 
	p <- p + annotate("text", x = -xlim*2/5, y = pos, label = paste("DN:", N_DN), size = 4, hjust = 1, color = "dodgerblue4", parse = TRUE) 
	N_UP = sum(DEGtable$DGE == "UP") 
	p <- p + annotate("text", x = xlim*2/5, y = pos, label = paste("UP:", N_UP), size = 4, hjust = 0, color = "red", parse = TRUE) 
 
	if(label){ 
		p <- p + geom_text_repel(data = diffM, aes(log2FC, -log10(pvalue), label = lable), 
		                		 size = 3, box.padding = unit(0.5, "lines"), 
		                		 point.padding = unit(0.8, "lines"),  
		                		 segment.color = "black",  
		                		 show.legend = FALSE, max.overlaps = 50) 
	} 
	return(p) 
} 
 
#' @title region_go 
#' @author Jiahao Wang 
#' @export
region_go <- function(bed, PDF){ 
 
    # check input 
    keyColumns = c("SYMBOL") 
    if(!all(keyColumns %in% colnames(bed))){ 
        message("Lost key columns: ", paste(setdiff(keyColumns, colnames(bed)), collapse = " ")) 
        stop() 
    } 
 
    genes = unique(bed[, "SYMBOL"]) 
    ego = runEGO(genes, 1, 1) 
    p <- dotplotGO(ego, out = PDF, padj = FALSE) 
    return(p) 
} 
 
#' @title region_kegg 
#' @author Jiahao Wang 
#' @export
region_kegg <- function(bed, PDF){ 
 
    # check input 
    keyColumns = c("SYMBOL") 
    if(!all(keyColumns %in% colnames(bed))){ 
        message("Lost key columns: ", paste(setdiff(keyColumns, colnames(bed)), collapse = " ")) 
        stop() 
    } 
 
    genes = unique(bed[, "SYMBOL"]) 
    ekegg = runKEGG(genes, 1, 1) 
    p <- dotplotKEGG(ekegg, out = PDF, padj = FALSE) 
    return(p) 
} 
 
#' @title region_plot_chromosome_chordal 
#' @author Jiahao Wang 
#' @export
region_plot_chromosome_chordal <- function(bed, PDF = NULL, Rplot = TRUE){ 
 
	loadp(circlize) 
 
    # check input 
    keyColumns = c("chr", "ST", "ED", "delta", "diff_type") 
    if(!all(keyColumns %in% colnames(bed))){ 
        message("Lost key columns: ", paste(setdiff(keyColumns, colnames(bed)), collapse = " ")) 
        stop() 
    } else { 
        if(!all(unique(bed[, "diff_type"]) %in% c("Hyper", "Hypo"))){ 
            message("Wrong diff_type, should be Hyper or Hypo.") 
            stop() 
        } 
    } 
 
    circlize_plot <- function(bed){ 
 
        bed = bed[bed$chr %in% paste0("chr", 1:22), ] 
        input_data2 = bed 
        input_data2[2:3] <- lapply(input_data2[2:3], as.numeric) 
        CGI_hyper = input_data2[input_data2$diff_type == 'Hyper',] 
        CGI_hypo = input_data2[input_data2$diff_type == 'Hypo',] 
        bed_list = list(CGI_hyper[c("chr", "ST", "ED")], CGI_hypo[c("chr", "ST", "ED")]) 
        bed_list2 = list(CGI_hyper[, c("chr", "ST", "ED", "delta")], CGI_hypo[, c("chr", "ST", "ED", "delta")]) 
 
        circos.clear()   
        human_cytoband = read.cytoband(species = "hg19")$df 
        human_cytoband2 = human_cytoband[!human_cytoband$V1 %in% c('chrX','chrY'),] 
        circos.par("start.degree" = 78, "gap.degree" = rep(c(2, 2), 12)) ##以1开头 
        circos.initializeWithIdeogram(plotType = NULL) 
        circos.track(ylim = c(0, 1), panel.fun = function(x, y) { 
        circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(3),  
                   ifelse(CELL_META$sector.index %in% c('chrX','chrY'),'',CELL_META$sector.index),facing = "clockwise",font=2,adj =c(0.4,0), 
                   cex = 1.1, niceFacing = TRUE) 
        }, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0) , bg.border = NA) 
        circos.genomicIdeogram(human_cytoband2,track.height = 0.06) 
 
        circos.genomicTrack(bed_list2,track.height = 0.16, 
                         panel.fun = function(region, value, ...) { 
                           circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                              border = ifelse(value[[1]] > 0, "red", "green"), ...) 
                           circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040") 
                           circos.yaxis(side = "left", at = seq(-1, 1, by = 0.5),#tick.length =0.1, 
                                        sector.index = get.all.sector.index()[1], labels.cex = 0.6) 
                         } , bg.border = NA) 
        circos.text(0,0, 'Mean \n methylation',sector.index = get.all.sector.index()[1],facing = 'clockwise',adj =c(0.5, -1.22), cex = 0.8, font = 2) 
        set_track_gap(mm_h(0.01)) 
        circos.genomicDensity(CGI_hyper, col = c("#FF000080"),track.height = 0.15,bg.border = NA)  ##画基因的密度 
        circos.text(0,0, '  Hyper \n Density',sector.index = get.all.sector.index()[1],facing = 'clockwise',adj =c(0.5,-0.8),cex = 0.8,font=2)   
        circos.genomicDensity(CGI_hypo, col = c("#0000FF80"), track.height = 0.15,bg.border = NA) 
        circos.text(0,0, '  Hypo \n Density',sector.index = get.all.sector.index()[1],facing = 'clockwise',adj =c(0.5,-0.8),cex = 0.8,font=2)  
    } 
 
    if(Rplot){ 
        suppressMessages(circlize_plot(bed)) 
    } 
 
    if(!is.null(PDF)){ 
        pdf(PDF, 10, 10) 
            suppressMessages(circlize_plot(bed)) 
        dev.off()        
    } 
    # invisible(circlize_plot(bed)) 
} 
 
#' @title region_plot_chromsome_bar 
#' @author Jiahao Wang 
#' @export
region_plot_chromsome_bar <- function(bed, PDF = NULL, Rplot = TRUE){ 
 
    # check input 
    keyColumns = c("chr", "ST", "ED", "diff_type") 
    if(!all(keyColumns %in% colnames(bed))){ 
        message("Lost key columns: ", paste(setdiff(keyColumns, colnames(bed)), collapse = " ")) 
        stop() 
    } else { 
        if(!all(unique(bed[, "diff_type"]) %in% c("Hyper", "Hypo"))){ 
            message("Wrong diff_type, should be Hyper or Hypo.") 
            stop() 
        } 
    } 
 
    loadp(ggplot2) 
    countT = reshape2::melt(table(bed$chr, bed$diff_type)) 
    colnames(countT) = c("chr", "diff_type", "count") 
    countT$chr = factor(countT$chr, levels = paste0("chr", 1:22)) 
    Max = max(countT$count) + max(countT$count) * 0.15; Max 
    p <- ggplot(countT, aes(chr, count, fill = diff_type)) 
    p <- p + geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.6) 
    p <- p + setTheme() 
    p <- p + setText(20, x.text.angle = 45, y.title.vjust = 1.5) 
    p <- p + setLimit(y.expand = c(0, 0), y.lim = c(0, Max)) 
    p <- p + labs(x = "", y = "Count of DMRs", title = "") 
    p <- p + scale_fill_manual(values = c("Hyper" = "#CE0977", "Hypo" = "#0457BD")) 
    p <- p + theme(legend.position = "top") 
    p <- p + theme(legend.title = element_blank()) 
    p <- p + theme(legend.text = element_text(size = 15)) 
 
    if(Rplot){ 
        print(p) 
    } 
 
    if(!is.null(PDF)){ 
        pdf(PDF, 10, 5) 
            print(p) 
        dev.off()        
    } 
    invisible(p) 
} 
 
#' @title region_plot_genomic_bar 
#' @author Jiahao Wang 
#' @export
region_plot_genomic_bar <- function(bed, PDF = NULL, Rplot = TRUE){ 
 
    # check input 
    keyColumns = c("chr", "ST", "ED", "annotation", "diff_type") 
    if(!all(keyColumns %in% colnames(bed))){ 
        message("Lost key columns: ", paste(setdiff(keyColumns, colnames(bed)), collapse = " ")) 
        stop() 
    } else { 
        if(!all(unique(bed[, "diff_type"]) %in% c("Hyper", "Hypo"))){ 
            message("Wrong diff_type, should be Hyper or Hypo.") 
            stop() 
        } 
    } 
 
    countT = reshape2::melt(table(bed$annotation, bed$diff_type)) 
    colnames(countT) = c("feature", "diff_type", "count") 
    # countT$count = sample(10:10000, 12) 
    Max = max(countT$count) + max(countT$count) * 0.15 
    p <- ggplot(countT, aes(x = feature, y = count, fill = feature, group = diff_type)) 
    p <- p + geom_bar(stat = 'identity', width = 0.6, colour = "black", position = position_dodge(width = 0.8)) 
    p <- p + labs(x = "", y = "Count of DMRs", title = "Genomic Count of DMRs") 
    p <- p + setTheme() 
    p <- p + setText(20, x.text.angle = 45, y.title.vjust = 1.5) 
    p <- p + setLimit(y.expand = c(0, 0), y.lim = c(0, Max)) 
    p <- p + theme(legend.position = "none") 
    p <- p + scale_fill_manual(values = rev(col.cluster2.1[1:(nrow(countT)/2)])) 
    p <- p + geom_text(aes(x = feature, y = count + Max * 0.05, label = count), size = 4.3, color = "black", position = position_dodge(width = 0.8)) 
    p <- p + facet_grid(. ~ diff_type) 
    p <- p + theme(strip.text.x = element_text(size = 20 * 0.8, color = "black")) 
    p <- p + theme(strip.background.x = element_rect(fill = "white", colour = "white")) 
    p <- p + theme(strip.switch.pad.grid = unit(0.5, "cm")) 
 
    if(Rplot){ 
        print(p) 
    } 
 
    if(!is.null(PDF)){ 
        pdf(PDF, 10, 5) 
            print(p) 
        dev.off()        
    } 
    invisible(p) 
} 
 
#' @title region_plot_genomic_stack 
#' @author Jiahao Wang 
#' @export
region_plot_genomic_stack <- function(bed, PDF = NULL, Rplot = TRUE){ 
 
    # check input 
    keyColumns = c("chr", "ST", "ED", "annotation", "diff_type") 
    if(!all(keyColumns %in% colnames(bed))){ 
        message("Lost key columns: ", paste(setdiff(keyColumns, colnames(bed)), collapse = " ")) 
        stop() 
    } else { 
        if(!all(unique(bed[, "diff_type"]) %in% c("Hyper", "Hypo"))){ 
            message("Wrong diff_type, should be Hyper or Hypo.") 
            stop() 
        } 
    } 
 
    loadp(ggplot2) 
    countT = reshape2::melt(table(bed$annotation, bed$diff_type)) 
    colnames(countT) = c("feature", "diff_type", "count") 
    countT$freq = countT$count / rep(aggregate(countT$count, by = list(countT$diff_type), sum)[, 2], each = nrow(countT) / 2) * 100 
    countT$feature = factor(countT$feature, levels = rev(c('Promoter', 'Exon', 'Intron', 'Intergenic', 'UTR', 'Downstream'))) 
 
    p <- ggplot(countT, aes(x = diff_type, weight = freq, fill = feature)) 
    p <- p + geom_bar( position = "stack") + coord_flip() + setTheme(grid.major = T) + setText(20) 
    p <- p + labs(x = "", y = "Percentage(%)", title = "DMRs Genomic Distribution") 
    p <- p + scale_fill_manual(values = rev(col.cluster2.1[1:(nrow(countT)/2)]), name = "Features", breaks = rev(levels(countT$feature)))    # p <- p + scale_fill_manual(values = col.cluster2.1[1:(nrow(countT)/2)], name = "Features", breaks = levels(countT$feature)) 
     
    if(Rplot){ 
        print(p) 
    } 
 
    if(!is.null(PDF)){ 
        pdf(PDF, 10, 5) 
            print(p) 
        dev.off()        
    } 
    invisible(p) 
} 
 
#' @title region_plot_length_bar 
#' @author Jiahao Wang 
#' @export
region_plot_length_bar <- function(bed, PDF = NULL, Rplot = TRUE){ 
 
    keyColumns = c("Width") 
    if(!all(keyColumns %in% colnames(bed))){ 
        message("Lost key columns: ", paste(setdiff(keyColumns, colnames(bed)), collapse = " ")) 
        stop() 
    } 
     
    Width = bed$Width 
    Width = Width[Width < 2000] 
 
    df <- data.frame(length = Width) 
    density <- density(df$length, adjust = 1/5) #调整值使直方图更光滑 
    df.density <- data.frame(x=density$x, y=density$y) 
 
    p <- ggplot(df, aes(x=length)) + theme_classic() + setText(20) 
    p <- p + geom_histogram(aes(y=..density..),  binwidth = 10, colour="black", fill="lightblue", alpha=0.5) 
    p <- p 
    p <- p + labs(title = "DMRs Length Distribution", x = "Length", y = "Density") 
 
    if(Rplot){ 
        print(p) 
    } 
 
    if(!is.null(PDF)){ 
        pdf(PDF, 10, 5) 
            print(p) 
        dev.off()        
    } 
 
    invisible(p) 
} 
 
#' @title save.graph 
#' @author Jiahao Wang 
#' @export
save.graph <- function(p, file, w = 4, h = 4, mar = NULL, PDF = TRUE, PNG = TRUE){ 
    if(!is.null(mar)){ 
        p <- p + theme(plot.margin = margin(mar[1], mar[2], mar[3], mar[4])) 
    } 
     
    if(PDF){ 
	    ggsave(paste0(file, ".pdf"), p, width = w, height = h)	 
    } 
 
    if(PNG){ 
	    ggsave(paste0(file, ".png"), p, width = w, height = h) 
    } 
} 
 
#' @title setBlank 
#' @author Jiahao Wang 
#' @export
setBlank <- function(x.title  = 1, 
                     x.text   = 1, 
                     x.ticks  = 1, 
                     y.title  = 1, 
                     y.text   = 1, 
                     y.ticks  = 1 
                    ){ 
     
    p = theme() 
    if(x.title) 
        p <- p + theme(axis.title.x = element_blank()) 
    if(x.text) 
        p <- p + theme(axis.text.x = element_blank()) 
    if(x.text) 
        p <- p + theme(axis.ticks.x = element_blank()) 
     
    if(y.title) 
        p <- p + theme(axis.title.y = element_blank()) 
    if(y.text) 
        p <- p + theme(axis.text.y = element_blank()) 
    if(y.ticks) 
        p <- p + theme(axis.ticks.y = element_blank()) 
     
    return(p) 
} 
 
#' @title setLegend 
#' @author Jiahao Wang 
#' @export
setLegend <- function(direction   = "vertical",  
					  position    = "right", # "none",   
					  title       = NA, #  
					  order       = NA, # c() 
					  labels      = NA, 
					  text.size   = 10, 
					  text.col    = "black", 
					  text.face   = "plain", 
					  text.angle  = 0, 
					  text.hjust  = 0, 
					  text.vjust  = 0, 
					  title.size  = 10, 
					  title.col   = "black", 
					  title.face  = "plain", 
					  title.angle = 0, 
					  title.hjust = 0, 
					  title.vjust = 0, 
					  key.size    = 2, 
					  key.color   = NA, 
					  key.fill    = NA, 
					  margin      = 10 
		 			 ){ 
 
  	# basic 
	p <- theme(legend.direction = direction,  
			   legend.position  = position, 
		  	   legend.text      = element_text(size = text.size, color = text.col, face = text.face, angle = text.angle, hjust = text.hjust, vjust = text.vjust), 
		  	   legend.key       = element_rect(fill = key.fill, color = key.color), 
		  	   legend.key.size  = unit(key.size, "cm")) 
 
	# magin 
  	if(length(margin) == 1){ 
	  	p <- p + theme(legend.margin = margin(rep(margin, 4), "pt")) 
  	} else{ 
	  	p <- p + theme(legend.margin = margin(margin[1], margin[2], margin[3], margin[4], "pt")) 
  	} 
 
  	# title 
  	if(is.na(title)){ 
	  	p <- p + theme(legend.title = element_text(size = title.size, color = title.col, face = title.face, angle = title.angle, hjust = title.hjust, vjust = title.vjust)) 
	} else if(isFALSE(title)){ 
  		p <- p + theme(legend.title = element_blank()) 
  	} else if(is.character(title)){ 
	  	p <- p + guides(color = guide_legend(title = title)) 
  	} 
 
	# order, labels and palette 
  	if(is.na(palette)){ 
		if(!is.na(order) & is.na(labels)){ 
			p <- p + scale_color_discrete(breaks = order) 
		} else if(is.na(order) & !is.na(labels)){ 
			p <- p + scale_color_discrete(labels = labels) 
		} else if(!is.na(order) & !is.na(labels)){ 
			p <- p + scale_color_discrete(breaks = order, labels = labels) 
		} 
  	} else{ 
  		if(!is.na(order) & is.na(labels)){ 
			p <- p + scale_color_manual(values = palette, breaks = order) 
		} else if(is.na(order) & !is.na(labels)){ 
			p <- p + scale_color_manual(values = palette, labels = labels) 
		} else if(!is.na(order) & !is.na(labels)){ 
			p <- p + scale_color_manual(values = palette, breaks = order, labels = labels) 
		} 
  	} 
 
	return(p) 
} 
 
#' @title setLimit 
#' @author Jiahao Wang 
#' @export
setLimit <- function(x.lim     = 0,  
					 y.lim     = 0, 
					 x.expand = 0, 
					 y.expand = 0 
					 ){ 
 
	if(length(x.lim) != 1){ 
		if(length(x.lim) == 2){ 
			if(length(x.expand) == 1){ 
				p <- scale_x_continuous(limits = c(x.lim[1], x.lim[2])) 
			} else{ 
				p <- scale_x_continuous(limits = c(x.lim[1], x.lim[2]), expand = x.expand) 
			} 
		} else{ 
			if(length(x.expand) == 1){ 
				p <- scale_x_continuous(limits = c(x.lim[1], x.lim[2]), breaks = seq(x.lim[1], x.lim[2], x.lim[3])) 
			} else{ 
				p <- scale_x_continuous(limits = c(x.lim[1], x.lim[2]), breaks = seq(x.lim[1], x.lim[2], x.lim[3]), expand = x.expand) 
			} 
		} 
		if(length(y.lim) != 1) 
			if(length(y.lim) == 2){ 
				if(length(y.expand) == 1){ 
					p <- p + scale_y_continuous(limits = c(y.lim[1], y.lim[2])) 
				} else{ 
					p <- p + scale_y_continuous(limits = c(y.lim[1], y.lim[2]), expand = y.expand) 
				} 
 
			} else{ 
				if(length(y.expand) == 1){ 
					p <- p + scale_y_continuous(limits = c(y.lim[1], y.lim[2]), breaks = seq(y.lim[1], y.lim[2], y.lim[3])) 
				} else{ 
					p <- p + scale_y_continuous(limits = c(y.lim[1], y.lim[2]), breaks = seq(y.lim[1], y.lim[2], y.lim[3]), expand = y.expand) 
				} 
 
			} 
	} else{ 
		if(length(y.lim) == 2){ 
			if(length(y.expand) == 1){ 
				p <- scale_y_continuous(limits = c(y.lim[1], y.lim[2])) 
			} else{ 
				p <- scale_y_continuous(limits = c(y.lim[1], y.lim[2]), expand = y.expand) 
			} 
 
		} else{ 
			if(length(y.expand) == 1){ 
				p <- scale_y_continuous(limits = c(y.lim[1], y.lim[2]), breaks = seq(y.lim[1], y.lim[2], y.lim[3])) 
			} else{ 
				p <- scale_y_continuous(limits = c(y.lim[1], y.lim[2]), breaks = seq(y.lim[1], y.lim[2], y.lim[3]), expand = y.expand) 
			} 
 
		} 
	} 
 
    return(p) 
} 
 
#' @title setText 
#' @author Jiahao Wang 
#' @export
setText <- function(title.size   = 10, 
                    title.face = "bold", 
                    title.hjust = 0.5, 
                    title.vjust = 0, 
                    x.title.size = NA, 
                    x.title.face = NA, 
                    x.title.vjust = 0.5, 
                    x.text.size  = NA, 
                    x.text.angle = 0, 
                    x.text.hjust = 0.5, 
                    x.text.vjust = 1, 
                    y.title.size = NA, 
                    y.title.vjust = 0, 
                    y.text.size  = NA, 
                    face         = "plain", 
                    y.title.face = NA, 
                    y.text.angle = 0, 
                    y.text.hjust = 1, 
                    graph.theme = "classic" 
){ 
	 
    checkNA <- function(x, y, ratio = NA){ 
         
        if(is.na(x)) 
            if(is.na(ratio)){ 
                x = y 
            } else{ 
                x = y * ratio 
            } 
        return(x) 
    } 
     
    x.title.size = checkNA(x.title.size, title.size, 0.8) 
    y.title.size = checkNA(y.title.size, title.size, 0.8) 
    x.text.size  = checkNA(x.text.size, title.size, 0.6) 
    y.text.size  = checkNA(y.text.size, title.size, 0.6) 
    x.title.face = checkNA(x.title.face, face) 
    y.title.face = checkNA(y.title.face, face) 
    if(x.text.angle){ 
	    if(x.text.hjust != 45){ 
	        x.text.hjust = 1 
	    } 
    } 
    if(y.text.angle) 
        y.text.hjust = 1 
    if(x.text.angle == 90){ 
    	x.text.vjust = 0.5 
    } 
 
    mytheme = switch(graph.theme, 
    		"bw" = theme_bw(), 
    		"classic"  = theme_classic() 
    		) 
    p <- mytheme 
    p <- p + theme(plot.title   = element_text(color = "black", size = title.size, face = title.face, hjust = title.hjust, vjust = title.vjust), 
               axis.title.x = element_text(color = "black", size = x.title.size, face = x.title.face, vjust = x.title.vjust), 
               axis.title.y = element_text(color = "black", size = y.title.size, face = y.title.face, vjust = y.title.vjust), 
               axis.text.x  = element_text(color = "black", size = x.text.size, angle = x.text.angle, hjust = x.text.hjust, vjust = x.text.vjust), 
               axis.text.y  = element_text(color = "black", size = y.text.size, angle = y.text.angle, hjust = y.text.hjust), 
               plot.margin  = margin(t = 5.5, r = 5.5, b = 5.5, l = 15, unit = "pt"), 
               legend.title = element_text(size = unit(x.title.size, "pt")), 
               legend.key.size = unit(x.text.size, "pt"), 
               legend.text = element_text(size = unit(x.text.size, "pt"))) 
    if(graph.theme == "blank"){ 
    	p <- p + theme(axis.line = element_blank(),  
    				   axis.text.x = element_blank(), 
    				   axis.text.y = element_blank(),  
    				   axis.ticks = element_blank(), 
    				   axis.title.x = element_blank(), 
    				   axis.title.y = element_blank(), 
    				   panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), 
    				   panel.grid.minor = element_blank(), plot.background = element_blank()) 
    } 
 
    return(p) 
} 
 
#' @title setTheme 
#' @author Jiahao Wang 
#' @export
setTheme <- function(theme.bw   = TRUE,  
					 grid.major = FALSE,  
					 grid.minor = FALSE 
					 ){ 
	p <- theme() 
	if(theme.bw) 
		p <- p + theme_bw()  
 
	if(!grid.major) 
		p <- p + theme(panel.grid.major = element_blank()) 
 
	if(!grid.minor) 
		p <- p + theme(panel.grid.minor = element_blank()) 
	 
	return(p) 
} 
 
#' @title show_cols 
#' @author Jiahao Wang 
#' @export
show_cols <- function(col, plot = FALSE){ 
 
    loadp(ggplot2) 
    p <- ggplot(data.frame(col = factor(col, levels = rev(col))), aes(col)) 
    p <- p + geom_bar(aes(fill = rev(col)), width = 1) 
    p <- p + scale_fill_manual(values = col) 
    p <- p + scale_y_continuous(expand = c(0, 0)) 
    p <- p + theme(legend.position = "none") + labs(x = NULL, y = NULL) 
    p <- p + theme(axis.text.x = element_blank(), panel.background = element_blank(), axis.ticks = element_blank()) 
    p <- p + theme(axis.text.y = element_text(face = "bold", size = 15)) 
    p <- p + coord_flip() 
    if(plot){ 
        pdf("show_cols.pdf", 3, 5) 
            p 
        dev.off() 
        cat("File save to 'show_cols.pdf.'\n") 
    } 
    return(p) 
} 
 
#' @title theme_blank 
#' @author Jiahao Wang 
#' @export
theme_blank <- function(){ 
    p <- theme(axis.line = element_blank(),  
               axis.text.x = element_blank(), 
               axis.text.y = element_blank(),  
               axis.ticks = element_blank(), 
               axis.title.x = element_blank(), 
               axis.title.y = element_blank(), 
               panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(), plot.background = element_blank()) 
    return(p) 
} 
 
#' @title venn 
#' @author Jiahao Wang 
#' @export
venn <- function(..., Rplot = TRUE, PDF = NULL, PPT = NULL,  
                 cex.label = 8, cex.num = 8,  
                 alpha.fill = 0.7, col.fill = c("#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0"),  
                 alpha.stroke = 1, size.stroke = 1, 
                 col.stroke = "black", 
                 show_percentage = FALSE, 
                 returnP = FALSE 
                ){ 
     
    data_list = list(...) 
    names(data_list) = getName(...) 
     
    # build venn table 
    Len = length(data_list) 
    vennM = matrix("", Len, Len, dimnames = list(getName(...), getName(...))) 
    for(i in 1:Len){ 
        for(j in i:Len){ 
            vennM[i, j] = length(intersect(data_list[[i]], data_list[[j]])) 
        } 
    } 
    resM = apply(vennM, 1, as.numeric) 
    dimnames(resM) = dimnames(vennM) 
     
    if(Len > 4){ 
        message("More than 4 objects, will not plot venn.") 
        return(resM) 
    } 
 
    # plot venn 
    require(ggvenn) 
    p <- suppressWarnings(ggvenn( 
        data = data_list, 
        show_percentage = show_percentage, 
        set_name_size = cex.label, 
        text_size = cex.num, 
        stroke_alpha = alpha.stroke, 
        stroke_size = size.stroke, 
        fill_alpha = alpha.fill, 
        fill_color = col.fill, 
        stroke_linetype = "longdash", 
    )) 
     
    if(Rplot){ 
        print(p) 
    } 
     
    if(!is.null(PDF)){ 
        pdf(PDF, 9, 9) 
            print(p) 
        dev.off() 
    } 
 
    if(!is.null(PPT)){ 
        require(export) 
        graph2ppt(p, file = PPT, width = 9, height = 9) 
    } 
     
    if(returnP){ 
	    invisible(p) 
    } else { 
	    invisible(resM) 
    } 
} 
 
 
