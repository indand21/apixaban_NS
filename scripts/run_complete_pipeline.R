# Run from the project root. Set NS_WORKERS to adjust local CPU concurrency.
source('scripts/load_project.R')
dir.create('output/revised',recursive=TRUE,showWarnings=FALSE)
dir.create('output/_run_logs',recursive=TRUE,showWarnings=FALSE)
skip_analysis <- '--figures-only' %in% commandArgs(trailingOnly=TRUE)
if(!skip_analysis) {
  dll<-file.path('src',paste0('qsp_core',.Platform$dynlib.ext))
  if(!file.exists(dll)||file.mtime('src/qsp_core.c')>file.mtime(dll)) {
    status<-system2(file.path(R.home('bin'),'R'),c('CMD','SHLIB','src/qsp_core.c'))
    if(status!=0) stop('C solver compilation failed; install the R-compatible build tools.')
  }
  source('scripts/run_all_tests.R')
  source('scripts/run_revised_analysis.R')
}
source('scripts/finalize_results.R')
source('scripts/build_revised_figures.R')
cat('Pipeline complete. Results are in output/revised; figures are in output/revised/figures.\n')
cat('Manuscript generation is not part of this repository.\n')
