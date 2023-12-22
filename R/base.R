#' @title clear 
#' @author Jiahao Wang 
#' @export
clear <- function(N = 80){ 
	cat(rep("\n", N)) 
} 
 
#' @title cat0 
#' @author Jiahao Wang 
#' @export
cat0 <- function(..., newline = TRUE) { 
	args <- list(...) 
	for(arg in args) cat(arg, sep = "") 
	if(newline) cat("\n") 
} 
 
#' @title `%!in%` 
#' @author Jiahao Wang 
#' @export
`%!in%` <- function(x,y) !("%in%"(x,y)) 
 
#' @title send_to_dingding 
#' @author Jiahao Wang 
#' @export
send_to_dingding <- function(info){ 
    cmd = readLines("/public3/home/scg0045/SHARE/code/auto_rrbs/utils/dingding_cmd.txt") 
    cmd = gsub("COMMAND", info, cmd) 
    system(cmd) 
} 
 
#' @title writeXlsx 
#' @author Jiahao Wang 
#' @export
writeXlsx <- function(..., dfList = NULL, sheetNames = NULL, file = NULL, saveRowCol = NULL, append = FALSE){ 
    if(is.null(file)){ 
        stop("Please specify 'file' for saving!") 
    } 
     
    loadp(xlsx) 
    if(is.null(dfList)){ 
        dfList = list2(...) 
    } 
     
    loadp(xlsx) 
    Len = length(dfList) 
     
    if(!is.null(sheetNames)){ 
        names(dfList) = sheetNames 
    } else if(sum(duplicated(toupper(names(dfList))) > 0)){ 
        message("Duplicated object names after touppered! Use '", paste(paste0("sheet", 1:Len), collapse = " "), "' instead.") 
        names(dfList) = paste0("sheet", 1:Len) 
    } 
         
    if(is.null(saveRowCol)){ 
        for(i in 1:Len){ 
            write.xlsx(dfList[[i]], sheetName = names(dfList)[i], file = file, append = ifelse(i != 1 | append, TRUE, FALSE)) 
        } 
    } else {      
        for(i in 1:Len){ 
            write.xlsx(dfList[[i]], sheetName = names(dfList)[i], file = file, append = ifelse(i != 1 | append, TRUE, FALSE),  
                       row.names = saveRowCol[[1]][i], col.names = saveRowCol[[2]][i]) 
        }         
    } 
} 
 
#' @title cd 
#' @author Jiahao Wang 
#' @export
cd  <- function(folder){ 
	setwd(folder) 
} 
 
#' @title checkMd5 
#' @author Jiahao Wang 
#' @export
checkMd5  <- function(masterFolder){ 
	cwd = getwd() 
	setwd(masterFolder) 
	files = lf("rf", masterFolder, ".md5$") 
	for(file in files){ 
		setwd(dirname(file)) 
		system(paste("md5sum -c", file)) 
		setwd(cwd) 
	} 
} 
 
#' @title checkMd5File 
#' @author Jiahao Wang 
#' @export
checkMd5File  <- function(fileIn, fileOut){ 
	if(subString(fileIn, 1) != "/"){ 
		stop("md5 file must be sbsolute path!") 
	} 
 
	cd(dirname(fileIn)) 
	lines = readLines(fileIn) 
	res = runPara(lines, checkMd5Line, cores = 64) 
 
	file_FAILED = subString(grep("FAILED", res, value = TRUE), 1, ":") 
	if(length(file_FAILED) == 0){ 
		cat("\n\tAll", length(lines), "files OK!\n") 
	} else { 
		cat("\n\tFind", length(file_FAILED), "files FAILED:\n") 
		for(x in file_FAILED) cat("\t\t", x, "\n") 
 
	} 
	write.line(res, file = fileOut) 
} 
 
#' @title checkMd5Line 
#' @author Jiahao Wang 
#' @export
checkMd5Line  <- function(line){ 
 
	md5_check = subString(line, 1, " ") 
	filePath = subString(line, 3, " ") 
 
	md5 = subString(system(paste("md5sum", filePath), intern = TRUE), 1, " ") 
	if(md5 == md5_check){ 
		res = "OK" 
	} else { 
		res = "FAILED" 
	} 
	res = paste0(filePath, ": ", res) 
	cat(res, "\n") 
	return(res) 
} 
 
#' @title checkPackages 
#' @author Jiahao Wang 
#' @export
checkPackages  <- function(..., stop = FALSE, install = FALSE){ 
	pkgs = getName(...) 
	uninstall <- setdiff(pkgs, installed.packages()) 
	if(length(uninstall) > 0){ 
		if(stop){ 
			message("\nError: Packages '", paste0(uninstall, collapse = ", "), "' not found!") 
			stop() 
		} else { 
			message("\nWarning: Packages '", paste0(uninstall, collapse = ", "), "' not found, install them ...!") 
			if(install){ 
				installer(pkgs = uninstall) 
			} else { 
				message("       Run 'installer(", paste0(uninstall, collapse = ", "), ")' to install them.\n") 
			} 
		}	 
	} 
} 
 
#' @title checkQuota 
#' @author Jiahao Wang 
#' @export
checkQuota  <- function(limit = 2){ 
	all  = system("lfs quota -gh scg0045 /public3 | sed -n '3p' | awk '{print $3}' | sed 's/T//'", intern = TRUE) 
	used = system("lfs quota -gh scg0045 /public3 | sed -n '3p' | awk '{print $2}' | sed 's/T//'", intern = TRUE)  
	remain =  as.numeric(all) - as.numeric(used) 
	if(remain < limit) 
		report(paste0("Remaining quota is less than ", limit, "T."), "W") 
} 
 
#' @title conna.test 
#' @author Jiahao Wang 
#' @export
conna.test  <- function(value, group, cutoff = 3){ 
 
	Ns = length(unique(group)) 
	idx = 1 
	for(i in 1:(Ns-1)){ 
 
		if(class(group) == "factor"){ 
			group1 = levels(group)[i] 
		}else{ 
			group1 = unique(group)[i] 
		} 
 
		data1 = na.omit(value[group == group1]) 
		for(j in (i+1):Ns){ 
 
			if(class(group) == "factor"){ 
				group2 = levels(group)[j] 
			}else{ 
				group2 = unique(group)[j] 
			} 
 
			data2  = na.omit(value[group == group2]) 
			if(length(data1) < cutoff | length(data2) < cutoff) 
				next 
			pair   = paste0(group1, "-", group2) 
			medianI  = mean(data1) 
			medianII = mean(data2) 
			# medianI  = median(data1) 
			# medianII = median(data2) 
			diff   = medianI - medianII 
			pvalue = wilcox.test(data1, data2, exact = FALSE)$p.value 
			xresT  = data.frame(pair, medianI, medianII, diff, pvalue) 
			# rownames(xresT) = pair 
			if(idx == 1){ 
				resT = xresT 
			} else{ 
				resT = rbind(resT, xresT) 
			} 
			idx = idx + 1 
		} 
		# resT = data.frame(resT) 
	} 
	 
	if(idx == 1) 
		return() 
	resT$padj = signif(p.adjust(resT$pvalue, method = "fdr"), 4) 
	return(na.omit(resT)) 
} 
 
#' @title countNA 
#' @author Jiahao Wang 
#' @export
countNA  <- function(Matrix){ 
 
	x = table(is.na(Matrix)) 
	return(c(ALL = sum(x), x)) 
} 
 
#' @title cp 
#' @author Jiahao Wang 
#' @export
cp  <- function(from, to){ 
	sapply(1:length(from), function(x) system(paste("cp", from[x], to[x]))) 
} 
 
 
#' @title ensureExists 
#' @author Jiahao Wang 
#' @export
ensureExists  <- function(folder){ 
	while(TRUE){ 
		if(dir.exists(folder)) 
			break 
		else 
			Sys.sleep(10) 
	} 
} 
 
#' @title ensureJobEnd 
#' @author Jiahao Wang 
#' @export
ensureJobEnd  <- function(Job = "Job"){ 
 
	check = TRUE 
	while(check){ 
		check = suppressWarnings(length(system(paste("qstat | grep", Job), intern = TRUE))) 
		if(check) Sys.sleep(1) 
	} 
} 
 
 
#' @title fillNA 
#' @author Jiahao Wang 
#' @export
fillNA  <- function(Matrix, by = 2, method = "median"){ 
 
	replaceNA <- function(x, method){ 
		if(method == "median"){ 
			x[is.na(x)] = median(x, na.rm = TRUE)  
		}else{ 
			x[is.na(x)] = mean(x, na.rm = TRUE)  
		} 
		return(x) 
	} 
 
	if(by == 2){ 
		res = apply(Matrix, 2, function(x) replaceNA(x, method))  
	} else{ 
		res = apply(Matrix, 1, function(x) replaceNA(x, method))  
	} 
 
	return(as.data.frame(res)) 
} 
 
#' @title getName 
#' @author Jiahao Wang 
#' @export
getName  <- function(...){ 
    as.character(substitute(list(...)))[-1] 
} 
 
#' @title getNline 
#' @author Jiahao Wang 
#' @export
getNline  <- function(file, zip = FALSE){ 
 
	if(!zip){ 
		size = as.numeric(system(paste("zcat", file, "| wc -l"), intern = TRUE)) 
	}else{ 
		size = as.numeric(system(paste("cat", file, "| wc -l"), intern = TRUE)) 
	} 
	return(size) 
} 
 
#' @title getSize 
#' @author Jiahao Wang 
#' @export
getSize <- function(file, unit = "m"){ 
 
	if(unit == "k"){ 
		size = as.numeric(strsplit(system(paste("du -sh -k", file), intern = TRUE), "\t")[[1]][1]) 
	} else{ 
		size = as.numeric(strsplit(system(paste("du -sh -m", file), intern = TRUE), "\t")[[1]][1]) 
	} 
	return(size) 
} 
 
#' @title grepr 
#' @author Jiahao Wang 
#' @export
grepr <- function(string, pattern, ignore.case = FALSE){ 
 
	res = matchPattern2(string, pattern, ignore.case = ignore.case)[, 1] 
	return(res) 
} 
 
#' @title greprl 
#' @author Jiahao Wang 
#' @export
greprl <- function(string, pattern, ignore.case = FALSE){ 
 
	res = rep(FALSE, nchar(string)) 
	idx = matchPattern2(string, pattern, ignore.case = ignore.case)[, 1] 
	res[idx] = TRUE 
	return(res) 
} 
 
#' @title hd 
#' @author Jiahao Wang 
#' @export
hd <- function(obj, x = 5, y = NULL){ 
 
	if(class(obj)[1] == "function") 
		return(head(obj)) 
 
	check_len <- function(idx, max){ 
		if(length(idx) > 1){ 
			if(max(idx) > max){ 
				idx_res = idx[idx <= max] 
			}else{ 
				idx_res = idx 
			} 
		}else{ 
			if(idx > max){ 
				idx_res = 1:max 
			}else{ 
				idx_res = 1:idx 
			} 
		} 
		return(idx_res) 
	} 
 
	if(is.null(y)) 
		y = x 
	 
	dims = is.null(dim(obj)) 
	if(!dims){ 
 
		cat("dim:", nrow(obj), "×", ncol(obj), "\n") 
 
		idx_x = check_len(x, nrow(obj)) 
		idx_y = check_len(y, ncol(obj)) 
		res = data.frame(obj[idx_x, idx_y]) 
		if(!is.null(rownames(obj))) 
			rownames(res) = rownames(obj)[idx_x] 
		if(!is.null(colnames(obj))) 
			colnames(res) = colnames(obj)[idx_y] 
	} else{ 
		cat("dim:", length(obj), "\n") 
		idx_x = check_len(x, length(obj)) 
		res = obj[idx_x] 
		if(!is.null(names(obj))) 
			names(res) = names(obj)[idx_x] 
	} 
	return(res) 
} 
 
#' @title installer 
#' @author Jiahao Wang 
#' @export
installer <- function(..., pkg = NULL, force = TRUE){ 
 
	if(is.null(pkg)){ 
	    pkg = as.character(substitute(list(...)))[-1] 
	} 
 
    CRAN = names(available.packages()[, 1]) 
    BIOC = BiocManager::available() 
    installed = names(installed.packages()[, 1]) 
    installed_target = intersect(pkg, installed) 
 
    if(length(installed_target) > 0){ 
        message("Package(s) '", paste(installed_target, collapse = " "), "' has beed installed, omit, use 'force = TRUE' to reinstall it.") 
    } 
         
    pkg_to_install = setdiff(pkg, installed) 
    if(length(pkg_to_install) == 0){ 
        message("No Package(s) to installed.") 
    } else { 
    	notFound = c() 
        for(single_pkg in pkg_to_install){ 
            message("\n\n") 
            message("******** INSTALLING ", single_pkg, " ********") 
            message("\n\n") 
            if(single_pkg %in% CRAN){ 
                install.packages(pkg) 
            } else if (single_pkg %in% BIOC) { 
                require("BiocManager") 
                BiocManager::install(pkg) 
            } else { 
            	notFound = c(notFound, single_pkg) 
            } 
        } 
	    message("Error: Package(s) '", paste(notFound, collapse = " "), "' not found in CRAN or Bioconductor.")    	 
	    message("Maybe they can installed in github: search 'github pkg_name'.")    	 
	    message("Or you give the error package names, check it!")    	 
    } 
} 
 
#' @title lf 
#' @author Jiahao Wang 
#' @export
lf <- function(type = "x", folder, pattern = "", remove = NULL, include.dirs = FALSE, ignore.case = FALSE){ 
	if(!type %in% c("rf", "f", "r", "x")) 
		stop("\nOnly support below types: 
			\t- rf: recursive = TRUE, full.names = TRUE 
			\t-  f: recursive = TRUE, full.names = FALSE 
			\t-  r: recursive = FALSE, full.names = TRUE 
			\t-  x: recursive = FALSE, full.names = FALSE\n") 
	res = switch(type, 
			rf = list.files(folder, pattern, full.names = TRUE, recursive = TRUE, include.dirs = include.dirs, ignore.case = ignore.case), 
			f  = list.files(folder, pattern, full.names = TRUE, recursive = FALSE, include.dirs = include.dirs, ignore.case = ignore.case), 
			r  = list.files(folder, pattern, full.names = FALSE, recursive = TRUE, include.dirs = include.dirs, ignore.case = ignore.case), 
			x  = list.files(folder, pattern, full.names = FALSE, recursive = FALSE, include.dirs = include.dirs, ignore.case = ignore.case) 
	) 
	if(!is.null(remove)){ 
		for(x in remove){ 
			res = res[!grepl(x, res)] 
		} 
	} 
	return(res) 
} 
 
#' @title list2 
#' @author Jiahao Wang 
#' @export
list2 <- function(...){ 
    resL = list(...) 
    names(resL) = getName(...) 
    return(resL) 
} 
 
#' @title load2 
#' @author Jiahao Wang 
#' @export
load2 <- function(fileIn){ 
	load(fileIn, verbose = TRUE, envir = parent.frame()) 
} 
 
#' @title loadp 
#' @author Jiahao Wang 
#' @export
loadp <- function(...){ 
    pkgs = as.character(substitute(list(...)))[-1] 
    suppressMessages(for(pkg in pkgs) require(pkg, character.only = TRUE)) 
} 
 
#' @title mean0 
#' @author Jiahao Wang 
#' @export
mean0 <- function(z) mean(z, na.rm = TRUE) 
 
#' @title numToMax 
#' @author Jiahao Wang 
#' @export
numToMax <- function(x){ 
    x = max(x) 
    ceiling(x / 10 ^ floor(log10(x) - 1) + 2) * 10 ^ floor(log10(x) - 1) 
} 
 
#' @title mv 
#' @author Jiahao Wang 
#' @export
mv <- function(from, to){ 
    sapply(1:length(from), function(x) system(paste("mv", from[x], to[x]))) 
} 
 
#' @title read 
#' @author Jiahao Wang 
#' @export
read <- function(type, file, header = FALSE, sep = "\t"){ 
	if(!type %in% c("rc", "r", "c", "x")) 
		stop("\nOnly support below types: 
			\t- rc: save rownames and colnames 
			\t-  r: only save  
			\t-  c: only save colnames 
			\t-  x: save neither rownames nor colnames\n") 
	switch(type, 
		rc = read.table(data, file, header = TRUE, col.names = TRUE, sep = "\t", quote = FALSE), 
		r  = read.table(data, file, header = TRUE, col.names = FALSE, sep = "\t", quote = FALSE), 
		c  = read.table(data, file, header = FALSE, col.names = TRUE, sep = "\t", quote = FALSE), 
		x  = read.table(data, file, header = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) 
	) 
} 
 
#' @title read.bedgraph 
#' @author Jiahao Wang 
#' @export
read.bedgraph <- function(file, showProgress = FALSE, nrows = -1){ 
 
	Tx = read.faster(file = file, header = FALSE, sep = "\t", skip = 1, showProgress = showProgress, nrows = nrows) 
	return(Tx) 
} 
 
#' @title read.fasta 
#' @author Jiahao Wang 
#' @export
read.fasta <- function(filePath){ 
 
	N = as.numeric(subString(system(paste("cat", filePath, "| grep '>' -n"), intern = TRUE), 1, ":")) 
	lines = readLines(filePath) 
	seqs = c() 
	for(i in 1:length(N)){ 
 
		if(i != length(N)){ 
			xlines = readLines(filePath)[(N[i] + 1):(N[i + 1] - 1)] 
		} else{ 
			Nlines = getNLines(file) 
			xlines = readLines(filePath)[(N[i] + 1):Nlines] 
		} 
		seqs = c(seqs, paste(xlines, collapse = "")) 
	} 
 
	Pos = subString(lines[N], 2, ">")   
	names(seqs) = Pos 
	res = as.list(seqs) 
	return(res) 
} 
 
#' @title read.faster 
#' @author Jiahao Wang 
#' @export
read.faster <- function(file, header = FALSE, sep = "\t", showProgress = FALSE, skip = 0, nrows = -1){ 
 
	suppressPackageStartupMessages(library("data.table")) 
	Tx = data.frame(fread(file = file, header = header, sep = sep, showProgress = showProgress, skip = skip, nrows = nrows)) 
	return(Tx) 
} 
 
#' @title read.gmt 
#' @author Jiahao Wang 
#' @export
read.gmt <- function (file){ 
 
    if (!grepl("\\.gmt$", file)[1]) { 
        stop("Pathway information must be a .gmt file") 
    } 
    geneSetDB = readLines(file) 
    geneSetDB = strsplit(geneSetDB, "\t") 
    names(geneSetDB) = sapply(geneSetDB, "[", 1) 
    geneSetDB = lapply(geneSetDB, "[", -1:-2) 
    geneSetDB = lapply(geneSetDB, function(x) { 
        x[which(x != "")] 
    }) 
    return(geneSetDB) 
} 
 
#' @title read.gr 
#' @author Jiahao Wang 
#' @export
read.gr <- function(file, metaNames = NULL){ 
 
	format = tail(strsplit(file, "\\.")[[1]], 1) 
	if(format == "bed"){ 
		res = toGR(read.faster(file, sep = "\t", showProgress = FALSE)) 
	} else if(format == "RData"){ 
		res = get(load(file)) 
	} else if(format == "bedGraph"){ 
		res = bedGraphToGR(read.bedgraph(file)) 
	} 
	if(!is.null(metaNames)) 
		colnames(mcols(res)) = metaNames 
	return(res) 
} 
 
#' @title reload 
#' @author Jiahao Wang 
#' @export
reload <- function(){ 
	rm("MyEvaluation", envir = .GlobalEnv) 
	source("~/.yyds2") 
} 
 
#' @title removeRowsAllNa 
#' @author Jiahao Wang 
#' @export
removeRowsAllNa <- function(x) x[apply(x, 1, function(y) any(!is.na(y))),] 
 
#' @title rep2 
#' @author Jiahao Wang 
#' @export
rep2 <- function(x, each){ 
    if(length(each) == 1){ 
        return(rep(x, each = each)) 
    }else{ 
        res = c() 
        for(i in 1:length(each)) res = c(res, rep(x[i], each = each[i])) 
        return(res) 
    } 
} 
 
#' @title report 
#' @author Jiahao Wang 
#' @export
report <- function(info, type = "I", print = TRUE, file = NULL, append = TRUE) { 
    type_mapping = c("I" = "INFO", "W" = "WARN", "E" = "ERROR") 
    if(!type %in% names(type_mapping)){ 
        stop(paste0("Unsupported typs: ", type, ". Available types: I(Info)|W(Warning)|E(Error)")) 
    } 
    type = type_mapping[type] 
    if(print){ 
        cat(paste0(type, " [", as.character(Sys.time()), "] ", info, "\n")) 
    } 
 
    if(!is.null(file)){  
        cat(paste0(type, " [", as.character(Sys.time()), "] ", info, "\n"), file = file, append = append) 
    } 
} 
 
#' @title sortBy 
#' @author Jiahao Wang 
#' @export
sortBy <- function(Tx, by, decreasing = FALSE){ 
 
	Tx = Tx[order(Tx[, by[1]], decreasing = decreasing[1]), ] 
	if(length(by) > 1){ 
		UNI_Name = unique(Tx[, by[1]]) 
		for(i in 1:length(UNI_Name)){ 
			xTx = Tx[Tx[, by[1]] == UNI_Name[i], ] 
			xTx_sorted = xTx[order(xTx[, by[2]], decreasing = decreasing[2]), ] 
			if(i == 1){ 
				res = xTx_sorted 
			} else{ 
				res = rbind(res, xTx_sorted) 
			} 
		} 
	} 
	return(res) 
} 
 
#' @title rm2 
#' @author Jiahao Wang 
#' @export
rm2 <- function(){ 
	rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv) 
} 
 
#' @title subset2 
#' @author Jiahao Wang 
#' @export
subset2 <- function(x, condition) { 
  condition_call <- substitute(condition) 
  r <- eval(condition_call, x, parent.frame()) 
  x[r, ] 
} 
 
#' @title toTitle 
#' @author Jiahao Wang 
#' @export
toTitle <- function(string){ 
 
	Title = toupper(subString(string, 1)) 
	res = paste0(Title, subString(tolower(string), 2:nchar(string))) 
    return(res) 
} 
 
#' @title uid 
#' @author Jiahao Wang 
#' @export
uid <- function(){ 
    loadp(ids) 
    return(uuid()) 
} 
 
#' @title write 
#' @author Jiahao Wang 
#' @export
write <- function(type, data, file, sep = "\t"){ 
	if(!type %in% c("rc", "r", "c", "x")) 
		stop("\nOnly support below types: 
			\t- rc: save rownames and colnames 
			\t-  r: only save  
			\t-  c: only save colnames 
			\t-  x: save neither rownames nor colnames\n") 
	switch(type, 
		rc = write.table(data, file, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE), 
		r  = write.table(data, file, row.names = TRUE, col.names = FALSE, sep = "\t", quote = FALSE), 
		c  = write.table(data, file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE), 
		x  = write.table(data, file, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) 
	) 
} 
 
#' @title write.gmt 
#' @author Jiahao Wang 
#' @export
write.gmt <- function(gsList, gmt_file){ 
 
	sink(gmt_file) 
	for (i in 1:length(gsList)){ 
 
		cat(names(gsList)[i]) 
		cat('\tNA\t') 
		cat(paste(gsList[[i]], collapse = '\t')) 
		cat('\n') 
	} 
	sink() 
} 
 
#' @title write.line 
#' @author Jiahao Wang 
#' @export
write.line <- function(x, file){ 
	write.table(x, file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t") 
} 
 
