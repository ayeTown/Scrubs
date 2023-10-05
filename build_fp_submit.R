

library(optparse)


# assumes all directories have a / at end
# must have a valid account number
# must have a logs directory
# assume you can run for more than one gene set full ranks file


option_list <- list(
	make_option(
		c("-r","--runjob"), 
		action = "store_true",
        default = FALSE, 
        help = "", 
        metavar = "character"
    ),
	make_option(
		c("-s","--system"),
		type = "character",
		default = "andes", 
        help = "", 
        metavar = "character"
    ),
    make_option(
		c("-a","--account"), 
		type = "character",
		default = "syb111", 
        help = "", 
        metavar = "character"
    ),
    make_option(
		c("-j","--jobname"), 
		type = "character",
		default = "functional-partitioning", 
        help = "", 
        metavar = "character"
    ),
    make_option(
		c("-d","--codedir"), 
		type = "character",
		default = NULL, 
        help = "", 
        metavar = "character"
    ),
    make_option(
		c("-c","--codename"), 
		type = "character",
		default = NULL, 
        help = "", 
        metavar = "character"
    ),
    make_option(
		c("-l","--logdir"), 
		type = "character",
		default = NULL, 
        help = "", 
        metavar = "character"
    ),
    make_option(
		c("-t","--time"), 
		type = "character",
		default = "2:00:00", # assuming default is andes 
        help = "", 
        metavar = "character"
    ),
  #   make_option(
		# c("-m","--multiplex"), 
		# type = "character",
		# default = NULL,
  #       help = "", 
  #       metavar = "character"
  #   ),
    make_option(
		c("-g","--genesets"), 
		type = "character",
		default = NULL, # needs to be a string separated by spaces
        help = "", 
        metavar = "character"
    ),
    make_option(
		c("-n","--genesetnames"), 
		type = "character",
		default = NULL, # needs to be a string separated by spaces
        help = "", 
        metavar = "character"
    ),
    make_option(
		c("-z","--thresholds"), 
		type = "character",
		default = NULL, # needs to be a string separated by spaces
        help = "", 
        metavar = "character"
    ),
    make_option(
		c("-o","--dirout"), 
		type = "character",
		default = NULL,
        help = "", 
        metavar = "character"
    )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


check_dir <- function(path) {

	if(!dir.exists(path)) {
		dir.create(path,showWarnings = FALSE)
	}

}


submit_job <- function(path) {

    command <- paste0("sbatch ",path)
    system(command)

}


build_fp_submit_script <- function(
		
		run = FALSE, # run or don't run the job script (default: FALSE)
		system = "andes", # system to use (default: andes)
		account = "syb111", # account to charge (default: syb111)
		job_name = NULL, # job name (default: NULL)
		code_dir = NULL, # directory for where to save code file (default: NULL)
		code_name = NULL, # file name for code file to save (default: NULL)
		log_dir = NULL, # directory for code log files (default: NULL)
		time_per_run = NULL, # time to use each node for (default: NULL)
		# multiplex = NULL, # location of multiplex network (default: NULL)
		gene_sets = NULL, # gene set full rank folders/files (default: NULL)
		gene_set_names = NULL, # gene set names to give each run (default: NULL)
		thresholds = NULL, # thresholds to use (default: NULL)
		dir_out = NULL # base directory to save the different runs (default: NULL)

	) {

	print("running functional partitioning job script creation...")

	gene_sets <- do.call("c",strsplit(gene_sets," "))
	gene_set_names <- do.call("c",strsplit(gene_set_names," "))
	thresholds <- do.call("c",strsplit(thresholds," "))

	if(system == "andes"|system == "frontier") { # if system is andes or frontier

		top <- c(
			"#!/bin/bash",
	        paste0("#SBATCH -A ",account),
	        paste0("#SBATCH -J ",job_name),
	        paste0("#SBATCH -N ",length(gene_sets)*length(thresholds)),
	        "#SBATCH \\-\\-mem=0",
	        paste0("#SBATCH -t ",time_per_run),
	        paste0("#SBATCH -o ",log_dir,job_name,".%J.out"),
	        paste0("#SBATCH -e ",log_dir,job_name,".%J.err")
		)
		if(system == "frontier") {
			top <- c(top,"#SBATCH -p batch")
		}
		top <- c(
			top,
			"\n",
	        "source /ccs/home/atown/Scripts/loadCondaOnAndes.sh",
			"conda activate functional-partitioning",
			"\n"
		)
		if(length(gene_sets)*length(thresholds) > 1) {
			parallel <- "srun -n1 -N1 -c8 --cpu-bind=threads --exclusive "
		} else {
			parallel <- ""
		}
		code_file_name <- paste0(code_name,".slurm")

	} else { # else the system is summit

		top <- c(
			"#!/bin/bash",
	        paste0("#BSUB -P ",account),
	        paste0("#BSUB -J ",job_name),
	        paste0("#BSUB -nnodes ",length(gene_sets)*length(thresholds)),
	        "#BSUB \\-\\-mem=0",
	        paste0("#BSUB -W ",time_per_run),
	        paste0("#BSUB -o ",log_dir,job_name,".%J.out"),
	        paste0("#BSUB -e ",log_dir,job_name,".%J.err"),
	        "\n",
	        "source /ccs/home/atown/Scripts/loadCondaOnSummit.sh",
			"conda activate functional-partitioning",
			"\n"
		)
		if(length(gene_sets)*length(thresholds) > 1) {
			parallel <- "jsrun -n1 -a1 -c8 -g1 -dpacked "
		} else {
			parallel <- ""
		}
		code_file_name <- paste0(code_name,".lsf")

	}

	fp_commands <- do.call("c",lapply(1:length(gene_sets),function(x){
		dir_outs <- paste0(dir_out,"gene-set-",gene_set_names[x],"_thresh-",thresholds,"/")
		lapply(dir_outs,check_dir)
		return(
			paste0(
				parallel,
				"functional_partitioning ",
				"--rwr-fullranks ",gene_sets[x]," ",
				"--partition ",
				"--cut-method hard ",
				"--cut-threshold ",thresholds," ",
				"--outdir ",dir_outs," ",
				"--no-plot ",
				"--verbose & \n"
			)
		)
	}))

	if(length(fp_commands) == 1) fp_commands <- gsub(" & \n","",fp_commands) else fp_commands <- c(fp_commands,"wait")

	script <- c(top,fp_commands)

	write.table(
		script,
		file = paste0(code_dir,code_file_name),
		sep = "\n",
		row.names = FALSE,
		col.names = FALSE,
		quote = FALSE
	)

	if(run) {
		submit_job(paste0(code_dir,code_file_name))
	}

	print("beep bop bop done.")

}


build_fp_submit_script(
		
	run = opt$runjob, 
	system = opt$system,
	account = opt$account,
	job_name = opt$jobname,
	code_dir = opt$codedir,
	code_name = opt$codename,
	log_dir = opt$logdir,
	time_per_run = opt$time,
	# multiplex = opt$multiplex,
	gene_sets = opt$genesets,
	gene_set_names = opt$genesetnames,
	thresholds = opt$thresholds,
	dir_out = opt$dirout

)


# c("-r","--runjob") 
# c("-s","--system") 
# c("-a","--account") 
# c("-j","--jobname") 
# c("-d","--codedir") 
# c("-c","--codename") 
# c("-l","--logdir")
# c("-t","--time")
# c("-m","--multiplex") 
# c("-g","--genesets") 
# c("-n","--genesetnames") 
# c("-z","--thresholds") 
# c("-o","--dirout")

