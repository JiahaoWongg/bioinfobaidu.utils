#' @title `%/%` 
#' @author Jiahao Wang 
#' @export
`%/%` <- function(x, y) as.numeric(x) / as.numeric(y) 
 
#' @title write.qc.bsmap 
#' @author Jiahao Wang 
#' @export
write.qc.bsmap <- function(fileIn, fileOut = NULL){ 
	lines = readLines(fileIn) 
	file_R1         = extractFields(lines, "Input read file #1", " ", 5) 
	file_R2          = extractFields(lines, "Input read file #2", " ", 5) 
	total_read_pairs = extractFields(lines, "total read pairs", " ", 10) 
	if(total_read_pairs == "pairs:") 
		total_read_pairs = extractFields(lines, "total read pairs", " ", 11) 
 
	aligned_pairs = extractFields(lines, "aligned pairs", " ", 3) 
	aligned_pairs_ratio = signif(aligned_pairs %/% total_read_pairs, 6) 
	unique_pairs = extractFields(lines, "unique pairs", " ", 7) 
	unique_pairs_ratio = signif(unique_pairs %/% total_read_pairs, 6) 
	suppressed_non_unique_pairs = extractFields(lines, "suppressed non-unique pairs", " ", 12) 
	suppressed_non_unique_pairs_ratio = signif(suppressed_non_unique_pairs %/% total_read_pairs, 6) 
 
	unpaired_read_1 = extractFields(lines, "unpaired read #1", " ", 4) 
	unpaired_read_1_ratio = signif(unpaired_read_1 %/% total_read_pairs, 6) 
	unpaired_unique_read_1 = extractFields(lines, "unpaired read #1", " ", 8) 
	unpaired_unique_read_1_ratio = signif(unpaired_unique_read_1 %/% total_read_pairs, 6) 
	suppressed_non_unique_read_1 = extractFields(lines, "unpaired read #1", " ", 13) 
	suppressed_non_unique_read_1_ratio = signif(suppressed_non_unique_read_1 %/% total_read_pairs, 6) 
 
	unpaired_read_2 = extractFields(lines, "unpaired read #2", " ", 4) 
	unpaired_read_2_ratio = signif(unpaired_read_2 %/% total_read_pairs, 6) 
	unpaired_unique_read_2 = extractFields(lines, "unpaired read #2", " ", 8) 
	unpaired_unique_read_2_ratio = signif(unpaired_unique_read_2 %/% total_read_pairs, 6) 
	suppressed_non_unique_read_2 = extractFields(lines, "unpaired read #2", " ", 13) 
	suppressed_non_unique_read_2_ratio = signif(suppressed_non_unique_read_2 %/% total_read_pairs, 6) 
 
	resT = data.frame(file_R1, file_R2, total_read_pairs, 
					  aligned_pairs, aligned_pairs_ratio, unique_pairs, unique_pairs_ratio,  
					  suppressed_non_unique_pairs, suppressed_non_unique_pairs_ratio, 
					  unpaired_read_1, unpaired_read_1_ratio, unpaired_unique_read_1, unpaired_unique_read_1_ratio, 
					  suppressed_non_unique_read_1, suppressed_non_unique_read_1_ratio, 
					  unpaired_read_2, unpaired_read_2_ratio, unpaired_unique_read_2, unpaired_unique_read_2_ratio, 
					  suppressed_non_unique_read_2, suppressed_non_unique_read_2_ratio) 
	Features = colnames(resT) 
	Value = as.character(resT[1, ]) 
	resT = data.frame(Features, Value) 
 
	if(is.null(fileOut)){ 
		fileOut = gsub("report.txt", "stat.txt", fileIn) 
	} 
	write("c", resT, fileOut) 
	invisible(resT) 
} 
 
#' @title extractFields 
#' @author Jiahao Wang 
#' @export
extractFields <- function(lines, anchor, sep, idx, rev = FALSE){ 
	line = lines[grepl(anchor, lines)] 
	res = subString(line, idx, sep, rev = rev) 
	return(res) 
} 
 
 
 
#' @title write.qc.trim_galore 
#' @author Jiahao Wang 
#' @export
write.qc.trim_galore <- function(fileIn, fileOut = NULL){ 
	lines = readLines(fileIn) 
	file_name = extractFields(lines, "Input filename", ": ", 2) 
	reads_raw = gsub(",", "", extractFields(lines, "Total reads processed", " ", 1, rev = TRUE)) 
	reads_with_adapter = gsub(",", "", extractFields(lines, "Reads with adapters", " ", 2, rev = TRUE)) 
	reads_with_adapter_ratio = signif(reads_with_adapter %/% reads_raw, 6) 
 
	reads_write = gsub(",", "", extractFields(lines, "Reads written", " ", 2, rev = TRUE)) 
	reads_write_ratio = signif(reads_write %/% reads_raw, 6) 
 
	bp_processed = gsub(",", "", subString(extractFields(lines, "Total basepairs processed", " ", 2, rev = TRUE), 1, " ")) 
	bp_qc_remove = gsub(",", "", extractFields(lines, "Quality-trimmed", " ", 3, rev = TRUE)) 
	bp_qc_remove_ratio = signif(bp_qc_remove %/% bp_processed, 6) 
	bp_write = gsub(",", "", extractFields(lines, "Total written", " ", 3, rev = TRUE)) 
	bp_write_ratio = signif(bp_write %/% bp_processed, 6) 
 
	resT = data.frame(file_name, reads_raw, reads_with_adapter, reads_with_adapter_ratio, 
										reads_write, reads_write_ratio, bp_processed, bp_qc_remove, 
										bp_qc_remove_ratio, bp_write, bp_write_ratio) 
	Features = colnames(resT) 
	Value = as.character(resT[1, ]) 
	resT = data.frame(Features, Value) 
 
	if(is.null(fileOut)){ 
		fileOut = gsub("report.txt", "stat.txt", fileIn) 
	} 
	write("c", resT, fileOut) 
	invisible(resT) 
} 
 
#' @title aggregateOverlaps 
#' @author Jiahao Wang 
#' @export
aggregateOverlaps <- function(gr, Intervals, ignore.strand = TRUE){ 
 
	suppressMessages(library("GenomicRanges")) 
	gr = gr[countOverlaps(gr, Intervals, type = "any", ignore.strand = ignore.strand) == 1] 
	x  = findOverlaps(gr, Intervals, type = "any", ignore.strand = ignore.strand) 
	cM = data.frame(mcols(gr)) 
	if(nrow(cM) == 0){ 
		stop("No common seqname or range in the two GRanges objects.") 
	} else{ 
		agT = aggregate(cM[, 1:2], by = list(subjectHits(x)), sum) 
		mgr = Intervals[as.integer(agT[,1])] 
		mcols(mgr) = agT[,2:3] 
		mgr$beta = signif(agT[, 2] / agT[, 3], 6) 
		if(exists("tag")) mgr$SID = tag 
	} 
	return(mgr) 
} 
 
#' @title aggregateOverlaps2 
#' @author Jiahao Wang 
#' @export
aggregateOverlaps2 <- function(gr, Intervals){ 
 
	suppressMessages(library("GenomicRanges")) 
	gr = gr[countOverlaps(gr, Intervals) == 1] 
	x  = findOverlaps(gr, Intervals) 
	cM = data.frame(mcols(gr)) 
	if(nrow(cM) == 0){ 
		stop("No common seqname or range in the two GRanges objects.") 
	} else{ 
		agT = aggregate(cM[, 1:2], by = list(subjectHits(x)), sum) 
		mgr = Intervals[as.integer(agT[,1])] 
		mcols(mgr) = agT[,2:3] 
		mgr$beta = signif(agT[, 2] / agT[, 3], 4) 
		if(exists("tag")) mgr$SID = tag 
	} 
	return(mgr) 
} 
 
#' @title arrayIntervalMethylation 
#' @author Jiahao Wang 
#' @export
arrayIntervalMethylation <- function(MeM, IntervalGR){ # by boss 
 
	source("/public1/home/scg5712/SHARE3/code/yyds/Methylation_workflow.R") 
	if(!exists("RnBeads_450K_hg19_GR")) 
		load("/public1/home/scg5712/SHARE3/iGenome/hg19/RData/RnBeads_450K_hg19_Probes_GR.RData") # RnBeads_450K_hg19_GR 
	cpgGR = RnBeads_450K_hg19_GR 
 
	# check seqlevels 
	if(length(intersect(seqlevels(IntervalGR), seqlevels(cpgGR))) == 0) 
		stop("No shared seqlevels found.\n") 
 
	# only keep CpGs that locate in intervals 
	Idx1 = countOverlaps(cpgGR, IntervalGR, type = "any", ignore.strand = TRUE) > 0 
	Idx2 = names(cpgGR) %in% rownames(MeM) 
	Idx  = Idx1 & Idx2 
 
	if(sum(Idx) == 0) 
		stop("No CpGs found in pre-defined regions.\n") 
	 
	cpgGR = cpgGR[Idx, ] 
 
	# keep shared CpGs 
	cNames = intersect(rownames(MeM), names(cpgGR)) 
	cpgGR  = cpgGR[cNames] 
	MeM    = MeM[cNames, ] 
	 
	# Sample-by-sample processing to control memory footprint 
	mergedM = matrix(NA, length(IntervalGR), ncol(MeM)) 
	x = findOverlaps(cpgGR, IntervalGR, type = "any", ignore.strand = TRUE) 
	aT  = aggregate(MeM, by = list(factor(subjectHits(x))), FUN="mean0") 
	aID = as.character(aT[,1]) 
	aMM = as.matrix(aT[,-1]) 
	mergedM[as.integer(aID), ] = aMM 
	dimnames(mergedM) = list(as.character(IntervalGR), colnames(MeM)) 
	 
	## return 
	return(mergedM) 
 
} 
 
#' @title ArrayIntervalMethylation 
#' @author Jiahao Wang 
#' @export
ArrayIntervalMethylation <- function(MM, IntervalGR, cpgGR){ 
 
	# check seqlevels 
	if(length(intersect(seqlevels(IntervalGR), seqlevels(cpgGR))) == 0) 
		stop("No shared seqlevels found.\n") 
 
	# only keep CpGs that locate in intervals 
	Idx1 = countOverlaps(cpgGR, IntervalGR, type = "any", ignore.strand = TRUE) > 0 
	Idx2 = names(cpgGR) %in% rownames(MM) 
	Idx  = Idx1 & Idx2 
 
	if(sum(Idx) == 0) 
		stop("No CpGs found in pre-defined regions.\n") 
	 
	cpgGR = cpgGR[Idx, ] 
 
	# keep shared CpGs 
	cNames = intersect(rownames(MM), names(cpgGR)) 
	cpgGR  = cpgGR[cNames] 
	MM    = MM[cNames, ] 
	 
	# Sample-by-sample processing to control memory footprint 
	mergedM = matrix(NA, length(IntervalGR), ncol(MM)) 
	x = findOverlaps(cpgGR, IntervalGR, type = "any", ignore.strand = TRUE) 
	aT  = aggregate(MM, by = list(factor(subjectHits(x))), FUN="mean0") 
	aID = as.character(aT[,1]) 
	aMM = as.matrix(aT[,-1]) 
	mergedM[as.integer(aID), ] = aMM 
	dimnames(mergedM) = list(as.character(IntervalGR), colnames(MM)) 
	 
	## return 
	return(mergedM) 
} 
 
#' @title buildGR 
#' @author Jiahao Wang 
#' @export
buildGR <- function(Chr, ST, ED, strand = "*", meta = NA){ 
 
	suppressMessages(library("GenomicRanges")) 
	Chr = as.character(Chr) 
	ST = as.numeric(ST) 
	ED = as.numeric(ED) 
	gr = GRanges(seqnames = Rle(Chr),  ranges = IRanges(ST, ED), strand = strand) 
	if(class(meta) == "data.frame"){ 
		mcols(gr) = meta	 
	} else if(class(meta) != "logical"){ 
		stop("'meta' must be a dataframe.") 
	} 
	return(gr) 
} 
 
#' @title CGItoGene 
#' @author Jiahao Wang 
#' @export
CGItoGene <- function(CGIs, DISTcutoff = 0){ 
	# check 
	if(!grepl("\\chr", CGIs[1])[1]){ 
	        stop("CGIs format: chr1:135124-135563:+") 
    } 
 
	load("/sibcb2/bioinformatics2/wangjiahao/Data/RData/CGIanno.RData") 
	Genes = as.character(CGIanno$Gene[CGIanno$CGI %in% CGIs]) 
	idx = Genes != "None" 
	GenesNoneRemoved = Genes[idx] 
	# if(length(which(!idx))) 
		# cat("Warning:", length(which(!idx)), "CGIs no corresponding genes, removed.\n") 
	resGene = unique(GenesNoneRemoved) 
	if(DISTcutoff){ 
		subAnno = CGIanno[CGIanno$CGI %in% CGIs,] 
		subAnnoNoNone = subAnno[Genes != "None",] 
		idx = as.character(subAnnoNoNone$Dist) <= DISTcutoff 
		Genes_too_far_removed = as.character(subAnnoNoNone$Gene[idx]) 
		# if(length(which(!idx))) 
			# cat("Warning:", length(which(!idx)), "mapped genes too far away, removed.\n") 
		resGene = unique(Genes_too_far_removed) 
	} 
	# cat("Result:", length(resGene), "genes acquired.\n") 
	return(resGene) 
} 
 
 
#' @title convertEnrichT 
#' @author Jiahao Wang 
#' @export
convertEnrichT <- function(type, enrichT){ 
	Tx = data.frame(enrichT) 
	if(type == "go"){ 
		Tx$geneID = as.character(unlist(sapply(Tx$geneID, function(x) paste0(unique(strsplit(x, "/")[[1]]), collapse = "/")))) 
	} else{ 
		Tx$Symbol = as.character(unlist(sapply(Tx$geneID, function(x) paste0(unique(sort(convertId(strsplit(x, "/")[[1]], "ENTREZID", "SYMBOL"))), collapse = "/")))) 
	} 
	Tx 
} 
 
#' @title convertId 
#' @author Jiahao Wang 
#' @export
convertId <- function(genes = NULL, from, to, OrgDb = "org.Hs.eg.db"){ 
 
	if(is.null(genes)){ 
		cat("\nAvailable gene types:\n") 
		cat("\n\tACCNUM, ALIAS, ENSEMBL, ENSEMBLPROT, ENSEMBLTRANS, ENTREZID, ENZYME\n") 
		cat("\n\tEVIDENCE, EVIDENCEALL, GENENAME, GO, GOALL, IPI, MAP, OMIM, ONTOLOGY\n") 
		cat("\n\tONTOLOGYALL, PATH, PFAM, PMID, PROSITE, REFSEQ, SYMBOL, UCSCKG, UNIGENE, UNIPROT\n\n") 
	}else{ 
		suppressPackageStartupMessages(library("clusterProfiler")) 
		res = suppressMessages(suppressWarnings( 
			  bitr(genes, fromType = from, toType = to, OrgDb = OrgDb))) 
		if(length(to) == 1) res = unique(res[, 2]) 
		return(res) 
	} 
} 
 
 
#' @title findGO 
#' @author Jiahao Wang 
#' @export
findGO <- function(pattern, method = "key"){ 
 
	if(!exists("GO_DATA")) 
		load("/sibcb2/bioinformatics2/wangjiahao/Data/RData/GO_DATA.RData") 
 
	if(method == "key"){ 
		pathways = cbind(GO_DATA$PATHID2NAME[grep(pattern, GO_DATA$PATHID2NAME)]) 
	} else if(method == "gene"){ 
		pathways = cbind(GO_DATA$PATHID2NAME[GO_DATA$EXTID2PATHID[[pattern]]]) 
	} 
 
	colnames(pathways) = "pathway" 
 
	if(length(pathways) == 0){ 
		cat("No results!\n") 
	} else{ 
		return(pathways) 
	} 
} 
 
 
#' @title getGO 
#' @author Jiahao Wang 
#' @export
getGO <- function(ID){ 
 
	if(!exists("GO_DATA")) 
		load("/sibcb2/bioinformatics2/wangjiahao/Data/RData/GO_DATA.RData") 
 
	allNAME = names(GO_DATA$PATHID2EXTID) 
	if(ID %in% allNAME){ 
		geneSet = GO_DATA$PATHID2EXTID[ID] 
		names(geneSet) = GO_DATA$PATHID2NAME[ID] 
		return(geneSet)		 
	} else{ 
		cat("No results!\n") 
	} 
} 
 
 
#' @title getKEGG 
#' @author Jiahao Wang 
#' @export
getKEGG <- function(ID){ 
 
	library("KEGGREST") 
 
	gsList = list() 
	for(xID in ID){ 
 
		gsInfo = keggGet(xID)[[1]] 
		if(!is.null(gsInfo$GENE)){ 
			geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])	 
			xgeneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])			 
			NAME = sapply(strsplit(gsInfo$NAME, " - "), function(x) x[1]) 
			names(xgeneSet) = NAME 
			gsList[NAME] = xgeneSet  
		} else{ 
			cat(" ", xID, "No corresponding gene set in specific database.\n") 
		} 
	} 
	return(gsList) 
} 
 
#' @title getOverlaps 
#' @author Jiahao Wang 
#' @export
getOverlaps <- function(gr, Intervals){ 
 
	x = countOverlaps(gr, Intervals, type = "any", ignore.strand = TRUE) 
	if(!length(x)){	 
		stop("No common seqname or range in the two GRanges objects.") 
	} 
	res = gr[x == 1] 
	return(res) 
} 
 
#' @title grListToGR 
#' @author Jiahao Wang 
#' @export
grListToGR <- function(grList){ 
 
	return(unlist(GRangesList(grList))) 
} 
 
#' @title grToBed 
#' @author Jiahao Wang 
#' @export
grToBed <- function(GR){ 
 
	suppressMessages(library("GenomicRanges")) 
	seqs   = as.character(GR) 
	chr    = sapply(strsplit(seqs, "[-:]"), function(z) z[1]) 
	ST     = start(GR) 
	ED     = end(GR) 
	Strand = strand(GR) 
	metas  = mcols(GR) 
	bed    = data.frame(chr, ST, ED, Strand, metas)	 
	return(bed) 
} 
 
#' @title grToFasta 
#' @author Jiahao Wang 
#' @export
grToFasta <- function(gr, file){ 
 
	hg19.fa = "/sibcb2/bioinformatics2/wangjiahao/iGenome/Bismark/hg19/hg19.fa" 
	pos = as.character(gr) 
	write.table(pos, "tmp.pos", sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE) 
	system(paste("samtools faidx -n 50", hg19.fa, "-r tmp.pos >", file)) 
	system(paste("rm tmp.pos & samtools faidx", file)) 
} 
 
#' @title keggtToGeneT 
#' @author Jiahao Wang 
#' @export
keggtToGeneT <- function(ekegg){ 
 
	geneList = strsplit(ekegg$Symbol, "/") 
	genes = sort(unique(unlist(geneList))) 
	resT = data.frame() 
	for(Gene in genes){ 
 
		Pathway = paste(ekegg$Description[sapply(geneList, function(x) Gene %in% x)], collapse = " / ")  
		xresT = data.frame(Gene, Pathway) 
		resT  = rbind(resT, xresT)  
	} 
	return(resT) 
} 
 
#' @title mergeGrStrand 
#' @author Jiahao Wang 
#' @export
mergeGrStrand <- function(gr){ 
 
	Tx = data.frame(gr) 
	CGI = paste0(seqnames(gr), ":", start(gr), "-", end(gr)) 
	mTx = aggregate(mcols(gr)[, 1:2], by = list(factor(CGI, levels = unique(CGI))), sum)[, -1] 
	mgr = toGR(unique(CGI)) 
	beta = signif(mTx[, 1] / mTx[, 2], 4) 
	mcols(mgr) = data.frame(mTx, beta, SID = unique(gr$SID)) 
	return(mgr) 
} 
 
#' @title PromoToGene 
#' @author Jiahao Wang 
#' @export
PromoToGene <- function(Colnames){ 
	gene = sapply(strsplit(Colnames, ";"), function(x) x[[2]]) 
	return(unlist(gene)) 
} 
 
#' @title removeCols 
#' @author Jiahao Wang 
#' @export
removeCols <- function(Matrix, cutoff = 0.85){ 
	 
	removeID = c() 
	for(i in 1:ncol(Matrix)){ 
 
		tag = colnames(Matrix)[i] 
		x = Matrix[, i] 
		ratio = sum(!is.na(x)) / length(x) 
		if(ratio < cutoff)  
			removeID = c(removeID, tag) 
	} 
	res = Matrix[, -match(removeID, colnames(Matrix))] 
 
	# report result 
	N = length(removeID) 
	Ratio = paste0(signif(N / ncol(Matrix), 3) * 100, "%")  
	cat("\t", paste0(N, "(", Ratio, ") variables removed.\n")) 
	return(res) 
} 
 
#' @title removeColsAllNa 
#' @author Jiahao Wang 
#' @export
removeColsAllNa <- function(x) x[, apply(x, 2, function(y) any(!is.na(y)))] 
 
#' @title runDEG 
#' @author Jiahao Wang 
#' @export
runDEG <- function(expM,  
				   s1,  
				   s2,  
				   FDRcutoff = 0.05,  
				   DIFFcutoff = 2,  
				   min = 3,  
				   padjust = TRUE 
				   ){ 
	Rownames = rownames(expM) 
	meanI = meanII = log2FC = pvalue = padj = size = rep(NA, nrow(expM)) 
 
	for(i in 1:nrow(expM)){ 
 
		x = as.numeric(na.omit(as.numeric(expM[i, s1]))) 
		y = as.numeric(na.omit(as.numeric(expM[i, s2]))) 
		if(length(x) < min | length(y) < min) 
			next 
 
		if((length(unique(x)) == 1 & length(unique(y)) == 1)) 
			if(unique(x) == unique(y)) 
				next 
 
		meanI[i]  = mean(x) 
		meanII[i] = mean(y) 
		log2FC[i] = log2(mean(x) / mean(y)) 
		pvalue[i] = suppressWarnings(wilcox.test(x, y)$p.value) 
		size[i]   = paste0(length(x), "/", length(y)) 
	} 
 
	Rownames = Rownames[!(is.na(meanI) | is.na(pvalue))] 
	DEGtable = na.omit(data.frame(meanI, meanII, log2FC, pvalue, size)) 
 
	DEGtable$padj = p.adjust(DEGtable$pvalue, method = "fdr", n = length(DEGtable$pvalue)) 
 
	DEG <- rep("NC", nrow(DEGtable)) 
	if(padjust){ 
		DEG[((DEGtable$padj) < FDRcutoff) & (DEGtable$log2FC > log2(DIFFcutoff))]  = "UP" 
		DEG[((DEGtable$padj) < FDRcutoff) & (DEGtable$log2FC < -log2(DIFFcutoff))] = "DN" 
	} else{ 
		DEG[((DEGtable$pvalue) < FDRcutoff) & (DEGtable$log2FC > log2(DIFFcutoff))]  = "UP" 
		DEG[((DEGtable$pvalue) < FDRcutoff) & (DEGtable$log2FC < -log2(DIFFcutoff))] = "DN" 
	} 
	DEGtable$DEG = DEG 
	rownames(DEGtable) = Rownames 
	return(DEGtable) 
} 
 
#' @title runDMR 
#' @author Jiahao Wang 
#' @export
runDMR <- function(betaM,  
				   s1,  
				   s2,  
				   FDRcutoff = 0.05,  
				   delta = 0.1,  
				   min = 3,  
				   padjust = TRUE 
				   ){ 
	Rownames = rownames(betaM) 
	meanI = meanII = delta = pvalue = padj = size = rep(NA, nrow(betaM)) 
 
	for(i in 1:nrow(betaM)){ 
 
		x = as.numeric(na.omit(as.numeric(betaM[i, s1]))) 
		y = as.numeric(na.omit(as.numeric(betaM[i, s2]))) 
		if(length(x) < min | length(y) < min) 
			next 
 
		if((length(unique(x)) == 1 & length(unique(y)) == 1)) 
			if(unique(x) == unique(y)) 
				next 
 
		meanI[i]  = mean(x) 
		meanII[i] = mean(y) 
		delta[i] = mean(x) - mean(y) 
		pvalue[i] = suppressWarnings(wilcox.test(x, y)$p.value) 
		size[i]   = paste0(length(x), "/", length(y)) 
	} 
 
	Rownames = Rownames[!(is.na(meanI) | is.na(pvalue))] 
	DMRtable = na.omit(data.frame(meanI, meanII, delta, pvalue, size)) 
 
	DMRtable$padj = p.adjust(DMRtable$pvalue, method = "fdr", n = length(DMRtable$pvalue)) 
 
	# DEG <- rep("NC", nrow(DMRtable)) 
	# if(padjust){ 
	# 	DEG[((DMRtable$padj) < FDRcutoff) & (DMRtable$delta > delta)]  = "Hyper" 
	# 	DEG[((DMRtable$padj) < FDRcutoff) & (DMRtable$delta < -delta)] = "Hypo" 
	# } else{ 
	# 	DEG[((DMRtable$pvalue) < FDRcutoff) & (DMRtable$delta > delta)]  = "Hyper" 
	# 	DEG[((DMRtable$pvalue) < FDRcutoff) & (DMRtable$delta < -delta)] = "Hypo" 
	# } 
	# DMRtable$DEG = DEG 
	rownames(DMRtable) = Rownames 
	return(DMRtable) 
} 
 
#' @title runEGO 
#' @author Jiahao Wang 
#' @export
runEGO <- function(sigG,  
				   pvalueCutoff = 0.05,  
				   qvalueCutoff = 0.05, 
				   universe = NULL,  
				   out = NULL, 
				   keyType = "ENSEMBL" , 
				   ont = "ALL", 
                   sig = TRUE, 
				   raw = FALSE 
				   ){ 
 
 
	if(sum(c("clusterProfiler", "org.Hs.eg.db", "msigdbr") %in% .packages()) != 3){ 
		message("Loading packages ...") 
		loadp(clusterProfiler, org.Hs.eg.db, msigdbr) 
	}		 
 
	TYPE = ifelse(subString(sigG[1], 1:4) == "ENSG", "ENSEMBL", ifelse(subString(sigG[1], 1) %in% LETTERS, "SYMBOL", "ENTREZID")) 
	if(TYPE != "ENSEMBL") sigG = convertId(sigG, TYPE, "ENSEMBL") 
 
	message("Running ONTOLOGY: ", ont, " ...") 
	suppressMessages( 
	ego <- enrichGO( 
		gene          = sigG, 
		keyType 	  = keyType, 
		OrgDb         = org.Hs.eg.db, 
		ont           = ont, 
		pAdjustMethod = "fdr", 
		pvalueCutoff  = pvalueCutoff, 
		qvalueCutoff  = qvalueCutoff, 
		readable      = TRUE 
	)) 
 
    if(sig){ 
        ego = ego[ego$pvalue < 0.05, ] 
    } 
 
	if(raw) 
		return(ego) 
 
	egoT = convertEnrichT("go", data.frame(ego)) 
	egoT$Count = sapply(strsplit(egoT$geneID, "/"), length) 
	if(!is.null(out)) 
		write.csv(egoT, file = out) 
 
	if(is.null(nrow(egoT))){ 
		message("No terms enriched") 
		egoT[1, ] = rep(NA, ncol(egoT)) 
	} else{ 
		message(paste(nrow(egoT), " terms enriched")) 
	} 
	return(egoT) 
} 
 
#' @title runEGO3 
#' @author Jiahao Wang 
#' @export
runEGO3 <- function(ENSEMBL){ 
 
	ego_BP = runEGO(ENSEMBL, ont = "BP") 
	ego_BP2 = simplifyEGO(ego_BP, 15)[1:10] 
 
	ego_CC = runEGO(ENSEMBL, ont = "CC") 
	ego_CC2 = simplifyEGO(ego_CC, 15)[1:10] 
 
	ego_MF = runEGO(ENSEMBL, ont = "MF") 
	ego_MF2 = simplifyEGO(ego_MF, 15)[1:10] 
 
	egoT2 = rbind(data.frame(ONTOLOGY = "BP", ego_BP2), data.frame(ONTOLOGY = "CC", ego_CC2), data.frame(ONTOLOGY = "MF", ego_MF2)) 
	data = egoT2[order(egoT2$Count, decreasing = TRUE), ] 
	data$Description = factor(data$Description, levels = rev(data$Description)) 
	data$GeneRatioValue = as.numeric(subString(data$GeneRatio, 1, "/")) / as.numeric(subString(data$GeneRatio, 2, "/")) 
	return(data) 
} 
 
#' @title runGSEA 
#' @author Jiahao Wang 
#' @export
runGSEA <- function(DEGtable  = NULL,  
					geneSet   = NULL,  
					rnkFile   = NULL,  
					chipFile  = NULL,  
					gmtFile   = NULL, 
					prefix    = "GSEA", 
					outFolder = "GSEA", 
					min 	  = 15, 
					max		  = 500,  
					verbo     = FALSE, 
					PDF       = FALSE 
					){ 
 
	if(prefix %in% list.files(outFolder)) 
		stop("\tPrefix already exists! please replace prefix.") 
 
	if((is.null(geneSet) & is.null(gmtFile)) | (!is.null(geneSet) & !is.null(gmtFile))) 
		stop("\tParameter 'geneSet' or 'gmtFile' must and only be supplied one.") 
 
	check = suppressMessages(suppressWarnings(require("xlsx"))) 
	if(!check) 
		warning("\tPackage 'xlsx' is not exist, results will save as txt format instead.") 
 
	tmpFolder = paste0(outFolder, "/tmp_", prefix) 
	dir.create(tmpFolder) 
 
	if(!is.null(DEGtable)){ 
 
		cat(" Prepare files ...\n") 
		Tx = read.table(DEGtable, row.names = 1, header = TRUE, sep = "\t") 
		Tx <- Tx[!is.na(Tx$pvalue), ] 
		Tx <- Tx[!is.na(Tx$log2FoldChange), ] 
 
		ID <- rownames(Tx) 
		ENSG <- sapply(strsplit(ID, " "), function(z) z[1]) 
		Symbol <- sapply(strsplit(ID, " "), function(z) z[2]) 
		pvalue <- as.numeric(Tx$pvalue) 
		pvalue[pvalue < 10^-300] <- 10^-300 
		zscore <- qnorm(pvalue/2, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)*sign(Tx$log2FoldChange) 
		oTable <- data.frame(ENSG, zscore) 
		rnkFile <- paste0(tmpFolder, "/", prefix, ".rnk") 
		write.table(oTable, file = rnkFile, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) 
 
		chip = na.omit(data.frame(ENSG, Symbol, Symbol)) 
		colnames(chip) = c("Probe Set ID", "Gene Symbol", "Gene Title") 
		chip = chip[!duplicated(as.character(chip$"Probe Set ID")), ] 
		chipFile = paste0(tmpFolder, "/", prefix, ".chip") 
		write.table(chip, file = chipFile, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE) 
 
		if(!is.null(geneSet)){ 
			gmtFile = paste0(tmpFolder, "/", prefix, ".gmt") 
			write.gmt(geneSet, gmtFile) 
		} 
	} 
 
	# run GSEA 
	cat(" Run GSEA ...\n") 
	cmd = "java -cp /sibcb2/bioinformatics/software/GSEA/gsea2-2.2.4.jar xtools.gsea.GseaPreranked" 
	cmd = paste(cmd, "-gmx", gmtFile) 
	cmd = paste(cmd, "-rnk", rnkFile) 
	cmd = paste(cmd, "-chip", chipFile) 
	cmd = paste(cmd, "-rpt_label", prefix) 
	cmd = paste(cmd, "-out", outFolder) 
	cmd = paste(cmd, "-scoring_scheme weighted -collapse true -mode Max_probe -norm None") 
	cmd = paste(cmd, "-nperm 1000 -include_only_symbols true") 
	cmd = paste(cmd, "-make_sets true -plot_top_x 20 -rnd_seed timestamp") 
	cmd = paste(cmd, "-set_min", min, "-set_max", max, "-zip_report false -gui false") 
	if(verbo){ 
		system(cmd) 
	}else{ 
		system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE) 
	}  
	system(paste("rm -r", tmpFolder)) 
 
	# get GSEA plots in PDF format 
	if(PDF){ 
		cat(" Get PDF format ...\n") 
		gseaRunFolder = list.files(outFolder, paste0(prefix, ".GseaPreranked"), full.names = TRUE) 
		if(length(gseaRunFolder) != 1) 
			stop("None or multiple GSEA run was found.\n") 
		cmd = paste("/sibcb2/bioinformatics/software/BcbioNG/anaconda/bin/gseapy replot -i", gseaRunFolder, "-o", gseaRunFolder) 
		system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE) 
	} 
 
	# summary 
	gseaRunFolder = list.files(outFolder, paste0(prefix, ".GseaPreranked"), full.names = TRUE) 
	newRunFolder  = gsub(basename(gseaRunFolder), prefix, gseaRunFolder) 
	system(paste("mv", gseaRunFolder, newRunFolder)) 
 
	neg = grep(".xls", list.files(newRunFolder, "gsea_report_for_na_neg", full.names = TRUE), value = TRUE) 
	pos = grep(".xls", list.files(newRunFolder, "gsea_report_for_na_pos", full.names = TRUE), value = TRUE) 
 
	posT = read.table(pos, header = TRUE, sep = "\t")[, 1:11] 
	negT = read.table(neg, header = TRUE, sep = "\t")[, 1:11] 
	resT = rbind(posT, negT)[, c(1, 4, 5, 7, 8)] 
	resT$State = c(rep("UP", nrow(posT)), rep("DN", nrow(negT))) 
 
	resT = resT[order(resT$FDR), ] 
	dimnames(resT) = list(1:nrow(resT), c("Pathyway", "Size", "NES", "Pvalue", "FDR", "State")) 
	resT[, 3:5] = apply(resT[, 3:5], 2, function(x) signif(x, 3)) 
 
	if(check){ 
		resL = split(resT, factor(resT$State, levels = c("UP", "DN"))) 
		if(!nrow(resL[["UP"]])) 
			resL[["UP"]] = matrix(NA, 1, 6, dimnames = list(1, colnames(resT))) 
		if(!nrow(resL[["DN"]])) 
			resL[["DN"]] = matrix(NA, 1, 6, dimnames = list(1, colnames(resT))) 
 
		write.xlsx(resL[["UP"]], file = paste0(newRunFolder, "/00_", prefix, "_GSEA_summary.xlsx"), sheetName = "UP", row.names = FALSE, append = TRUE) 
		write.xlsx(resL[["DN"]], file = paste0(newRunFolder, "/00_", prefix, "_GSEA_summary.xlsx"), sheetName = "DN", row.names = FALSE, append = TRUE) 
	} else{ 
		write.table(resT, file = paste0(newRunFolder, "/00_", prefix, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE) 
	} 
	invisible(resL) 
} 
 
#' @title runKEGG 
#' @author Jiahao Wang 
#' @export
runKEGG <- function(sigG,  
				    pvalueCutoff = 0.05,  
				    qvalueCutoff = 0.05, 
				    organism = "hsa", 
					out = NULL 
					){ 
 
	suppressPackageStartupMessages(library("clusterProfiler")) 
	suppressPackageStartupMessages(library("org.Hs.eg.db")) 
 
	TYPE = ifelse(subString(sigG[1], 1:4) == "ENSG", "ENSEMBL", ifelse(subString(sigG[1], 1) %in% LETTERS, "SYMBOL", "ENTREZID")) 
	if(TYPE != "ENTREZID") sigG = convertId(sigG, TYPE, "ENTREZID") 
    # if(TYPE != "ENTREZID") sigG = as.numeric(convertId(sigG, TYPE, "ENTREZID")) 
 
	suppressMessages( 
	ekegg <- enrichKEGG( 
		gene          = sigG, 
		organism      = organism, 
		keyType 	  = "kegg", 
		pAdjustMethod = "fdr",  
		pvalueCutoff  = pvalueCutoff, 
		qvalueCutoff  = qvalueCutoff, 
		use_internal_data = FALSE 
	)) 
 
	ekeggT = convertEnrichT("kegg", ekegg) 
	if(!is.null(out)) 
		write.csv(ekeggT, file = out) 
 
	if(is.null(nrow(ekeggT))){ 
		cat("No terms enriched\n") 
		ekeggT[1, ] = rep(NA, ncol(ekeggT)) 
	} else{ 
		cat(paste(nrow(ekeggT), "terms enriched\n")) 
	} 
	return(ekeggT) 
} 
 
#' @title runPara 
#' @author Jiahao Wang 
#' @export
runPara <- function(Vector, fun, combine = "c", cores = NULL){ 
 
	suppressPackageStartupMessages(library(doParallel)) 
 
	if(is.null(cores)){ 
		cores = detectCores() - 2 
	} else { 
		if(cores > detectCores()){ 
			cores = detectCores() - 2 
		} 
	} 
	# cat("runPara:", cores, "cores will be used.\n") 
 
	registerDoParallel(cores) 
	res <- foreach(v = Vector, .combine = combine) %dopar% fun(v) 
	invisible(res) 
} 
 
#' @title seqkitStat2xlsx 
#' @author Jiahao Wang 
#' @export
seqkitStat2xlsx <- function(fileIn, fileOut){ 
	loadp(xlsx) 
	Tx = read.table(fileIn, header = TRUE) 
	By = list(factor(subString(basename(Tx$file), 1, "_"))) 
	Tx_1 = aggregate(Tx[, 4:5], by = By, sum)[, -1] 
	Tx_2 = aggregate(Tx[, 6], by = By, min)[, -1] 
	Tx_3 = aggregate(Tx[, 7], by = By, mean)[, -1] 
	Tx_4 = aggregate(Tx[, 8], by = By, max)[, -1] 
 
	Tx_2 = aggregate(Tx[, c(6:8, 14:15)], by = By, mean)[, -1] 
	statT = cbind(unique(By[[1]]), Tx_1, Tx_2) 
	statT$sum_len = signif(Tx_1$sum_len / 10^9, 4) 
	colnames(statT) = c("sample_ID", "total_reads", "total_bases(G)", "min_len", "avg_len", "max_len", "Q20(%)", "Q30(%)") 
	rownames(statT) = NULL 
	table2xlsx(statT, row.names = FALSE, file = fileOut) 
} 
 
#' @title simplifyAnno 
#' @author Jiahao Wang 
#' @export
simplifyAnno <- function(anno){ 
 
    res = rep(NA, length(anno)) 
    features = c("Promoter", "Exon", "Intron", "Intergenic", "UTR", "Downstream") 
    for(feature in features) res[grep(feature, anno)] = feature 
    res = factor(res, levels = features) 
    return(res) 
} 
 
#' @title simplifyEGO 
#' @author Jiahao Wang 
#' @export
simplifyEGO <- function(ego, N = 20){ 
 
	if(ego@ontology == "GOALL"){ 
 
 
 
	} else{ 
		if(nrow(ego) <= N){ 
			check = FALSE 
		} else{ 
			check = TRUE 
		} 
		cutoff = 0.7 
		while(check){ 
			ego = simplify(ego, cutoff = cutoff) 
			check = !nrow(ego) <= N 
			cutoff = cutoff - 0.05 
		}	 
	} 
 
	cat(nrow(ego), "terms Remain.\n") 
	return(ego) 
} 
 
#' @title summaryGSEA 
#' @author Jiahao Wang 
#' @export
summaryGSEA <- function(gseaRunFolder){ 
 
	resL = list() 
	names = c("pos", "neg") 
	for(name in names){ 
		reportFile = list.files(gseaRunFolder, paste0("gsea_report_for_na_", name), full.name = TRUE) 
		reportFile = reportFile[grep("xls$", reportFile)] 
		reportTx = read.table(reportFile, header = TRUE, sep = "\t") 
 
		pathways = as.character(reportTx$NAME)[reportTx$GS.DETAILS != ""] 
		coreEnricheds = c() 
		for(pathway in pathways){ 
 
			file  = paste0(gseaRunFolder, "/", pathway, ".xls") 
			xresT = read.table(file, header = TRUE, sep = "\t") 
			coreEnriched = as.character(xresT$PROBE[xresT$CORE.ENRICHMENT == "Yes"]) 
			coreEnriched = paste(coreEnriched, collapse = "/") 
			coreEnricheds = c(coreEnricheds, coreEnriched) 
		} 
		outTx = reportTx[reportTx$GS.DETAILS != "", c(1, 4, 6, 7, 8)] 
		outTx$coreEnrichment = coreEnricheds 
		colnames(outTx) = c("pathway", "size", "NES", "pvalue", "FDR", "coreEnrichment") 
		resL[[name]] = outTx 
	} 
	return(resL) 
} 
 
#' @title summaryHomer 
#' @author Jiahao Wang 
#' @export
summaryHomer <- function(outFolder){ 
 
	homerFolder = paste0(outFolder, "/homerResults") 
	xFiles = list.files(homerFolder, ".motif$") 
	xFiles = xFiles[-grep("similar", xFiles)] 
	xFiles = xFiles[-grep("RV", xFiles)] 
	xFiles = xFiles[order(as.numeric(gsub("\\.", "", gsub("motif", "", xFiles))))] 
	texts  = sapply(paste0(homerFolder, "/", xFiles), readLines) 
	chunks = sapply(texts, function(x) strsplit(x[1], "[\t]")) 
 
	motif = sapply(chunks, function(x) subString(x[1], 2, ">")) 
	match = sapply(chunks, function(x) subString(subString(x[2], 2, "BestGuess:"),  1, "/")) 
	score = sapply(chunks, function(x) rev(strsplit(x[2], "[()]")[[1]])[1]) 
	count = sapply(chunks, function(x) subString(x[6], 3, "[T:()]")) 
	ratio = sapply(chunks, function(x) subString(x[6], 2, "[()]")) 
	p_value = sapply(chunks, function(x) subString(x[6], 2, "P:")) 
 
	xresT = data.frame(motif,  
					   match,  
					   score = as.numeric(score),  
					   count = as.numeric(count), 
					   ratio_perc = as.numeric(gsub("%", "", ratio)),  
					   p_value = as.numeric(p_value) 
					   ) 
	rownames(xresT) = gsub(".motif", "", basename(rownames(xresT))) 
	return(xresT) 
} 
 
#' @title summaryHomerKnown 
#' @author Jiahao Wang 
#' @export
summaryHomerKnown <- function(outFolder){ 
 
	knownFolder = paste0(outFolder, "/knownResults") 
	xFiles = list.files(knownFolder, ".motif$") 
	xFiles = xFiles[order(as.numeric(gsub("\\.motif", "", gsub("known", "", xFiles))))] 
	texts  = sapply(paste0(knownFolder, "/", xFiles), readLines) 
	chunks = sapply(texts, function(x) strsplit(x[1], "[\t]")) 
 
	motif = sapply(chunks, function(x) subString(x[1], 2, ">")) 
	TF    = sapply(chunks, function(x) subString(x[2], 1, "/")) 
	count = sapply(chunks, function(x) subString(x[6], 3, "[T:()]")) 
	ratio = sapply(chunks, function(x) subString(x[6], 2, "[()]")) 
	p_value = sapply(chunks, function(x) subString(x[6], 2, "P:")) 
 
	xresT = data.frame(motif,  
					   TF,  
					   count = as.numeric(count), 
					   ratio_perc = as.numeric(gsub("%", "", ratio)),  
					   p_value = as.numeric(p_value) 
					   ) 
	rownames(xresT) = gsub("\\.motif", "", basename(rownames(xresT))) 
	return(xresT) 
} 
 
#' @title table2xlsx 
#' @author Jiahao Wang 
#' @export
table2xlsx <- function(Tx, file, sheetName = "sheet1", row.names = TRUE, col.names = TRUE){ 
	loadp(xlsx) 
	wb <- createWorkbook(type = "xlsx") 
	#TABLE_COL_STYLE <- CellStyle(wb) + Font(wb) + Alignment(wrapText = TRUE, horizontal = "ALIGN_RIGHT") 
	TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold = TRUE) +  
		Alignment(wrapText = TRUE, horizontal = "ALIGN_CENTER") +  
		Border(color = "black", position=c("BOTTOM", "LEFT", "RIGHT", "TOP"),  
		       pen = c("BORDER_THIN", "BORDER_THIN", "BORDER_THIN", "BORDER_THIN")) 
	TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold = TRUE) +  
		Alignment(wrapText = TRUE, horizontal = "ALIGN_CENTER") +  
		Border(color = "black", position=c("BOTTOM", "LEFT", "RIGHT", "TOP"),  
		       pen = c("BORDER_THIN", "BORDER_THIN", "BORDER_THIN", "BORDER_THIN")) 
	sheet <- createSheet(wb, sheetName = sheetName) 
	addDataFrame(Tx, sheet, col.names = col.names, row.names = row.names,  
	           startRow = 1, startColumn = 1, #colStyle = TABLE_COL_STYLE, 
	           colnamesStyle = TABLE_COLNAMES_STYLE, rownamesStyle = TABLE_ROWNAMES_STYLE) 
	if(row.names){ 
		autoSizeColumn(sheet, colIndex = c(1:(ncol(Tx) + 1))) 
	} else { 
		autoSizeColumn(sheet, colIndex = c(1:ncol(Tx))) 
	} 
	saveWorkbook(wb, file) 
} 
 
#' @title TCGAtranslateID 
#' @author Jiahao Wang 
#' @export
TCGAtranslateID <- function(file_ids, legacy = FALSE){ 
 
	suppressPackageStartupMessages(library("GenomicDataCommons")) 
    info = files(legacy = legacy) %>% 
           filter( ~ file_id %in% file_ids) %>% 
           select('cases.samples.submitter_id') %>% 
    id_list = lapply(info$cases, function(x) x[[1]][[1]][[1]]) 
    barcodes_per_file = sapply(id_list,length) 
	res = data.frame(file_id = rep(ids(info),barcodes_per_file),  
					 submitter_id = unlist(id_list)) 
    return(res) 
} 
 
#' @title toGR 
#' @author Jiahao Wang 
#' @export
toGR <- function(obj, strandIdx = NA){ 
 
	suppressMessages(library("GenomicRanges")) 
	type = class(obj)[1] 
	if(type == "GRanges"){ 
		GR = obj 
	} else if(type == "list"){ 
		GR = unlist(GRangesList(grList)) 
	} else if(type == "character"){ 
 
		chr = sapply(strsplit(obj, "[-:]"), function(z) z[1]) 
		ST = as.numeric(sapply(strsplit(obj, "[-:]"), function(z) z[2])) 
		ED = as.numeric(sapply(strsplit(obj, "[-:]"), function(z) z[3])) 
		SD = as.character(sapply(strsplit(obj, "[-:]"), function(z) z[4])) 
		GR = buildGR(chr, ST, ED, SD) 
	} else if(type == "data.frame"){ 
 
		chr = obj[, 1] 
		ST  = obj[, 2] 
		ED  = obj[, 3] 
 
		if(ncol(obj) > 3){ 
			if(obj[1, 4] %in% c("+", "-", "*")){ 
				if(ncol(obj) > 4){ 
 
					meta = obj[, 5:ncol(obj)] 
					strand = obj[, 4] 
					GR = buildGR(chr, ST, ED, strand, meta)   
				} else{ 
					strand = obj[, 4] 
					GR = buildGR(chr, ST, ED, strand)   
				} 
			} else if(!is.na(strandIdx)){ 
 
				strand = obj[, strandIdx] 
				meta = obj[, -c(1:3, strandIdx)] 
				GR = buildGR(chr, ST, ED, strand, meta)   
		    } else{ 
				meta = obj[, -c(1:3), drop = FALSE] 
				GR = buildGR(chr, ST, ED, meta = meta)   
		    } 
		} else{ 
			GR = buildGR(chr, ST, ED) 
		} 
	} 
	return(GR) 
} 
 
#' @title tosymbol 
#' @author Jiahao Wang 
#' @export
tosymbol <- function(geneid){symbol <- bitr(unlist(strsplit(geneid,split='/')),fromType='ENTREZID',toType='SYMBOL',OrgDb="org.Hs.eg.db")[,2];paste(symbol,collapse='/')} 
