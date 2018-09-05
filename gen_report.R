library(knitr)
library(rmarkdown)
args = commandArgs(TRUE)

params = list(source='upload',tv_file=args[1],start_day = as.numeric(args[2]),analysis_day=as.numeric(args[4]),
              sel_group=unlist(strsplit(args[5],',')),ref_group=args[6], posthoc_cmps=args[7], analysis_var=args[3],tgi_def=args[8])
Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc")
rmarkdown::render('tv_analysis.Rmd', output_file='InVivo_stat_analysis.html',
                  params = params,envir = new.env(parent = globalenv()))