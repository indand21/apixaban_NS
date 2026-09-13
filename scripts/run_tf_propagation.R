source('scripts/load_project.R')
stages<-c('Normal','NS_Mild','NS_Moderate','NS_Severe')
set.seed(20260910);tf_draws<-rlnorm(200,log(2),.95)
tf_tasks<-expand.grid(draw=seq_along(tf_draws),stage=stages,stringsAsFactors=FALSE)
tf_signature<-digest::digest(c(tools::md5sum(c(list.files('R',full.names=TRUE),list.files('data',full.names=TRUE),'src/qsp_core.c')),tf_draws))
tf_checkpoint<-file.path('output/revised/checkpoints',paste0('tf_',tf_signature))
dir.create(tf_checkpoint,recursive=TRUE,showWarnings=FALSE)
tf_workers<-as.integer(Sys.getenv('NS_WORKERS','6'))
tf_cluster<-parallel::makePSOCKcluster(tf_workers,outfile='output/revised/tf_workers.log')
parallel::clusterExport(tf_cluster,c('tf_tasks','tf_draws','tf_checkpoint'))
parallel::clusterEvalQ(tf_cluster,{source('scripts/load_project.R');load_fast_solver();NULL})
tf_results<-parallel::parLapplyLB(tf_cluster,seq_len(nrow(tf_tasks)),function(j){
  file<-file.path(tf_checkpoint,sprintf('case_%04d.csv',j))
  if(file.exists(file)) return(read.csv(file))
  i<-tf_tasks$draw[j];st<-tf_tasks$stage[j]
  pr<-prepare_patient(reference_patient(st,100),tf_draws[i]);o<-optimize_prepared(pr)
  ans<-data.frame(draw=i,stage=st,TF_pM=tf_draws[i],dose=o$optimal_dose,score=o$max_TI,
    boundary=o$boundary,weak_objective=o$max_TI<.001,score_etp=o$TI_etp)
  write.csv(ans,paste0(file,'.tmp'),row.names=FALSE);stopifnot(file.rename(paste0(file,'.tmp'),file))
  cat('Completed TF case',j,'\n');ans
},chunk.size=1)
parallel::stopCluster(tf_cluster)
write.csv(do.call(rbind,tf_results),'output/revised/tf_prior_propagation.csv',row.names=FALSE)
cat('TF PRIOR PROPAGATION COMPLETE\n')
