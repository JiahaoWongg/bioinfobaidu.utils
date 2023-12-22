#' @title countPattern 
#' @author Jiahao Wang 
#' @export
countPattern  <- function(string, pattern, ignore.case = FALSE){ 
 
	res = nrow(matchPattern2(string, pattern, ignore.case = ignore.case)) 
	return(res) 
} 
 
#' @title countProp 
#' @author Jiahao Wang 
#' @export
countProp  <- function(string, patterns, ignore.case = FALSE){ 
 
	if(ignore.case){ 
		string  = tolower(string) 
		pattern = tolower(pattern) 
	} 
	res = countPattern(string, patterns, ignore.case = FALSE) / nchar(string) 
	return(res) 
} 
 
#' @title gsub2 
#' @author Jiahao Wang 
#' @export
gsub2 <- function(from, to, obj){ 
	for(i in 1:length(from)){ 
		if(length(to) == 1){ 
			obj = gsub(from[i], to, obj) 
		}else{ 
			obj = gsub(from[i], to[i], obj) 
		} 
	} 
	obj 
} 
 
#' @title matchPattern2 
#' @author Jiahao Wang 
#' @export
matchPattern2 <- function(string, pattern, ignore.case = FALSE){ 
 
	string = as.character(string) 
	if(ignore.case){ 
		string = tolower(string) 
		pattern = tolower(pattern)		 
	} 
	seqs = subsplit(string, nchar(pattern)) 
	start = which(seqs == pattern) 
	end = start + nchar(pattern) - 1 
	res = data.frame(start, end) 
	return(res) 
} 
 
#' @title paster 
#' @author Jiahao Wang 
#' @export
paster <- function(strs1, strs2, sortBy = 2, start = "", end = "", sep = ""){ 
 
    res1 = unlist(lapply(strs1, function(x) paste0(x, sep, strs2))) 
    res2 = paste0(start, res1, end) 
    # if(sortBy == 2){ 
	   #  res2 = as.character(sapply(strs2, function(x) grep(x, res2, value = TRUE))) 
    # } 
    return(res2) 
} 
 
#' @title printSeq 
#' @author Jiahao Wang 
#' @export
printSeq <- function(Seq, binwidth = 50){ 
	N_line = nchar(Seq) / binwidth 
	if(N_line > as.integer(N_line)){ 
		N =  as.integer(N_line) + 1 
	} else{ 
		N =  as.integer(N_line) 
	} 
 
	bases = strsplit(Seq, "")[[1]] 
	for(i in 0:(N-1)){ 
		xseq = paste0(bases[(1 + binwidth*i):(binwidth + i*binwidth)], collapse = "")  
		if(i != (N-1)){ 
			xseq = paste0(bases[(1 + binwidth*i):(binwidth + i*binwidth)], collapse = "")  
		} else{ 
			xseq = paste0(bases[(1 + binwidth*i):length(bases)], collapse = "")  
		} 
		cat(xseq, "\n") 
	} 
} 
 
#' @title revCompString 
#' @author Jiahao Wang 
#' @export
revCompString <- function(strings){ 
 
	strings = toupper(strings)	 
	idx = data.frame(to = c("T", "C", "G", "A", "N"), row.names = c("A", "G", "C", "T", "N")) 
	res = as.character(sapply(strings, function(x) revString(paste(as.character(idx[strsplit(x, "")[[1]], "to"]), collapse = "")))) 
	# chars = strsplit(strings, "")[[1]] 
	# res = revString(paste(as.character(idx[chars, "to"]), collapse = "")) 
	return(res) 
} 
 
#' @title revString 
#' @author Jiahao Wang 
#' @export
revString <- function(strings){ 
 
	strings = as.character(strings) 
	res = as.character(lapply(strings, function(x) paste(rev(strsplit(x, "")[[1]]), collapse = ""))) 
	return(res) 
} 
 
#' @title strsplit2 
#' @author Jiahao Wang 
#' @export
strsplit2 <- function(x, sep, idx, psep = NULL){ 
	res = sapply(x, function(y) strsplit(y, sep)[[1]][idx], simplify = FALSE, USE.NAMES = FALSE) 
	if(!is.null(psep)){ 
		sapply(res, function(z) paste(z, collapse = psep), simplify = TRUE, USE.NAMES = FALSE)	 
	} else{ 
		as.character(res) 
	} 
} 
 
#' @title subString 
#' @author Jiahao Wang 
#' @export
subString <- function(strings, idx, sep = NULL, rev = FALSE, collapse = NULL){ 
 
	strings = as.character(strings) 
	if(is.null(sep)){ 
		if(rev){ 
			res = as.character(sapply(strings, function(x) paste(rev(rev(strsplit(x, "")[[1]])[idx]), collapse = ""))) 
		} else { 
			res = as.character(sapply(strings, function(x) paste(strsplit(x, "")[[1]][idx], collapse = ""))) 
		} 
	} else{ 
		if(rev){ 
			res = sapply(strsplit(strings, sep), function(x) paste(rev(rev(x)[idx]), collapse = collapse)) 
		} else { 
			res = sapply(strsplit(strings, sep), function(x) paste(x[idx], collapse = collapse))	 
		} 
	} 
	return(res) 
} 
 
#' @title subsplit 
#' @author Jiahao Wang 
#' @export
subsplit <- function(string, N = 2){ 
 
	string = as.character(string) 
	seqs = c() 
	for(i in 1:(nchar(string) - N + 1)){ 
		seqs = c(seqs, subString(string, i:(i + N - 1))) 
	} 
	return(seqs) 
} 
