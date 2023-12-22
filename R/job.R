 
#' @title buildSlurm 
#' @author Jiahao Wang 
#' @export
buildSlurm  <- function( 
					  cmd, # Command 
					  jobName,     # job namne  
					  shellFolder, # Folder to save slurm file 
					  logFolder,   # Folder to save log file 
					  N = 1,       # Number of node 
					  n = 1,       # Number of thread 
					  c = 1        # Number of core 
					  ){ 
 
	sink(paste0(shellFolder, "/", jobName, ".sh")) 
	cat("#!/bin/bash\n") 
	cat("#SBATCH -p amd_512\n") 
	cat("#SBATCH -N", N, "\n") 
	cat("#SBATCH -n", n, "\n") 
	cat("#SBATCH -c", c, "\n") 
	cat("#SBATCH -J", jobName, "\n") 
	cat("#SBATCH -o", paste0(logFolder, "/", jobName, ".out"), "\n") 
	cat("#SBATCH -e", paste0(logFolder, "/", jobName, ".err"), "\n") 
	cat("\n") 
 
	cat("echo 'Job running ...'\n") 
	cat(cmd, "||", paste0("{ echo 'Job error:'; tail -n 40 ", logFolder,"/", jobName, ".err; exit 0; }\n")) 
	cat("echo 'Job succeed!'\n") 
	cat("\n") 
	sink() 
} 
 
#' @title checkJobsNumber 
#' @author Jiahao Wang 
#' @export
checkJobsNumber  <- function(cycle, Nlimit = 40, sleep = 10, gap = 0.5){ 
 
	# deleteEqw() 
	if((cycle %% 10) == 0) 
	    Sys.sleep(gap) 
 
	check = TRUE 
	while(check){ 
		 
		USER  = system("echo $USER", intern = TRUE) 
		N_job = system(paste("qstat | grep", USER, "| wc -l"), intern = TRUE) 
		check = !as.numeric(N_job) < Nlimit 
		if(check){ 
			cat("\t", Nlimit, "jobs are running, wait", sleep, "secs\n") 
			deleteEqw() 
			Sys.sleep(sleep) 
		} 
	} 
} 
 
#' @title deleteEqw 
#' @author Jiahao Wang 
#' @export
deleteEqw  <- function(){ 
 
	jobNm = system("qstat | grep Eqw | awk '{print $3}'", intern = TRUE) 
	jobID = system("qstat | grep Eqw | awk '{print $1}'", intern = TRUE) 
	sapply(jobID, function(x) system(paste("qdel", x), ignore.stdout = TRUE)) 
	return(jobNm) 
} 
 
#' @title deleteJobs 
#' @author Jiahao Wang 
#' @export
deleteJobs  <- function(pattern = ""){ 
 
	jobNm = system(paste("qstat | grep", pattern,  "| awk '{print $3}'"), intern = TRUE) 
	jobID = system(paste("qstat | grep", pattern,  "| awk '{print $1}'"), intern = TRUE) 
	sapply(jobID, function(x) system(paste("qdel", x), ignore.stdout = TRUE)) 
 
} 
 
#' @title delJob 
#' @author Jiahao Wang 
#' @export
delJob  <- function(jobID){ 
	jobs = subString(subString(system("squeue", intern = TRUE)[-1], 2, "           "), 1, "   ") 
	if(jobID %in% jobs){ 
		system(paste("scancel", jobID)) 
		cat("\tJob", jobID, "deleted\n") 
	} else{ 
		cat("\tError: Job", jobID, "not exists\n") 
	} 
} 
 
#' @title JobEnd 
#' @author Jiahao Wang 
#' @export
JobEnd <- function(flag){ 
	check = suppressWarnings(system(paste("squeue --format='%.18i %.9P %.80j %.8u %.8T %.10M %.9l %.6D %R' | grep RUNNING | awk '{print $3}' | grep", flag), intern = TRUE)) 
	if(length(check) == 0){ 
		return(TRUE) 
	} else { 
		return(FALSE) 
	} 
} 
 
#' @title limitJobNumber 
#' @author Jiahao Wang 
#' @export
limitJobNumber <- function(Nlimit = 480, sleep = 60, flag = "amd_512", quiet = FALSE){ 
 
	check = TRUE 
	while(check){	 
		N_job = as.numeric(system(paste("squeue --format='%.18i %.9P %.80j %.8u %.8T %.10M %.9l %.6D %R' | grep", flag, "| grep -v COMPLETI | wc -l"), intern = TRUE)) - 1 
		check = !(as.numeric(N_job) < as.numeric(Nlimit)) 
		if(check){ 
			if(!quiet) cat("\t", Nlimit, "jobs are running, wait", sleep, "secs\n") 
			Sys.sleep(sleep) 
		} 
	} 
	Sys.sleep(1) 
} 
 
#' @title loopControlor 
#' @author Jiahao Wang 
#' @export
loopControlor <- function( 
						  Ns,           # Number of total cycle 
						  i,            # Current cycle index 
						  tag,          # Job name 
						  Nlimit = 480, # Maximun number of task 
						  sleep = 60,   # Waitting time, second 
						  flag = "amd_512", 
						  quota = 2, 
						  quiet = FALSE 
						  ){ 
	if(!quiet) cat("\tSubmitting job", tag, "Remain", Ns - i, "\n") 
	checkQuota(quota) 
	limitJobNumber(Nlimit, sleep, flag, quiet) 
} 
 
#' @title submitCMD 
#' @author Jiahao Wang 
#' @export
submitCMD <- function(cmd, jobName, c = 1, n = 1, N = 1, exclusive = NULL, verbose = FALSE, 
					  check = TRUE){ 
	 
	if(check){ 
		logFile = paste0("./log/", jobName, ".out") 
		if(file.exists(logFile)){ 
			cmdTmp = paste("grep 'Job succeed'", logFile) 
			succeedInfo = suppressWarnings(as.character(system(cmdTmp, intern = TRUE))) 
			if(length(succeedInfo) != 0){ 
				cat(paste0("\tWarnning: Job '", jobName, "' had succeed already, omit!\n")) 
				return("!") 
			} else { 
				cmdTmp = paste("squeue --format='%.8i %.9P %.60j %.8u %.8T %.10M %.9l %.6D %R' | awk '{print $3}'") 
				jobNames = suppressWarnings(as.character(system(cmdTmp, intern = TRUE))[-1]) 
				if(jobName %in% jobNames){ 
					cat("Job", jobName, "is running, please delete it or submit latter!\n") 
				} 
			} 
		} 
	} 
 
	code = "/public3/home/scg0045/SHARE/code/utils/submitCMD.sh" 
	if(is.null(exclusive)){ 
		CMD = paste0(code, " -c ", c, " -N ", N, " -n ", n, " -J ", jobName, " \"", cmd, "\"") 
	} else { 
		CMD = paste0(code, " -e -c ", c, " -N ", N, " -n ", n, " -J ", jobName, " \"", cmd, "\"") 
	} 
 
	info = system(CMD, intern = TRUE) 
	if(verbose){ 
		cat("\t", info, "\n") 
	} 
	invisible(info) 
} 
 
#' @title submitJob 
#' @author Jiahao Wang 
#' @export
submitJob <- function( 
					  cmd,         # Command 
					  jobName,     # job namne  
					  shellFolder, # Folder to save slurm file 
					  logFolder,   # Folder to save log file 
					  N = 1,       # Number of node 
					  n = 1,       # Number of thread 
					  c = 1        # Number of core 
					  ){ 
 
	shellFile = paste0(shellFolder, "/", jobName, ".sh") 
	sink(shellFile) 
	cat("#!/bin/bash\n") 
	cat("#SBATCH -p amd_512\n") 
	cat("#SBATCH -N", N, "\n") 
	cat("#SBATCH -n", n, "\n") 
	cat("#SBATCH -c", c, "\n") 
	cat("#SBATCH -J", jobName, "\n") 
	cat("#SBATCH -o", paste0(logFolder, "/", jobName, ".out"), "\n") 
	cat("#SBATCH -e", paste0(logFolder, "/", jobName, ".err"), "\n") 
	cat("\n") 
	cat(cmd) 
	cat("\n") 
	sink() 
 
	system(paste("chmod 770", shellFile)) 
	jobID = subString(system(paste("sbatch", shellFile), intern = TRUE), 4, " ") 
	cat("\tJob", jobName, paste0("(ID: ", jobID, ") submmitted\n")) 
	cat(as.character(Sys.time()), "\tJob", jobName, paste0("(ID: ", jobID, ") submmitted\n"), file = "/public3/home/scg5695/study/job", append = TRUE) 
} 
 
#' @title waitJobEnd 
#' @author Jiahao Wang 
#' @export
waitJobEnd <- function(flag, verbose = TRUE){ 
	if(verbose) cat("\tWaiting job end, flag:", flag, "\n") 
	Sys.sleep(10) 
	while(TRUE){ 
		check = suppressWarnings(system(paste("squeue --format='%.18i %.9P %.80j %.8u %.8T %.10M %.9l %.6D %R' | grep -v COMPLETI | awk '{print $3}' | grep", flag), intern = TRUE)) 
		if(!length(check)){ 
			# cat("Job end!\n") 
			break 
		}else{ 
			Sys.sleep(10) 
		} 
	} 
} 
