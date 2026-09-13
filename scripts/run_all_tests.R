dir.create('output/_run_logs',recursive=TRUE,showWarnings=FALSE)
files<-list.files('tests',pattern='^test_.*[.]R$',full.names=TRUE)
results<-list()
for(f in files) {
  log<-file.path('output/_run_logs',paste0(tools::file_path_sans_ext(basename(f)),'.log'))
  status<-system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',shQuote(f)),stdout=log,stderr=log)
  txt<-readLines(log,warn=FALSE)
  # Several original test files only printed FAIL; enforce failure centrally.
  failure<-any(grepl('(^FAIL |[.]?[.][.] FAIL|[1-9][0-9]* FAIL ===|[1-9][0-9]* failed$)',txt))
  results[[f]]<-data.frame(test_file=f,exit_status=status,pass=status==0&&!failure)
  cat(f,if(status==0&&!failure)'PASS' else 'FAIL','\n')
}
res<-do.call(rbind,results)
write.csv(res,'output/revised/test_suite_status.csv',row.names=FALSE)
if(any(!res$pass)) quit(status=1)
