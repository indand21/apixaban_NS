source('scripts/load_project.R')
p<-prepare_patient(reference_patient());sp<-c(.coag_species_names,.plt_species_names)
for(a in c(1e-10,1e-12,1e-14)) {
  x<-simulate_fast(p$ic,p$plt,p$rates,p$coupling,10,atol=a,rtol=a*100)
  pp<-c(p$rates,p$plt,p$coupling,list(Cp_free_nM=10,Ki_apixaban=.25))
  pp$alpha_assembly<-0;pp$alpha_catalysis<-0
  xr<-as.data.frame(deSolve::ode(y=c(p$ic,get_platelet_initial_conditions(p$plt)),times=x$time,
    func=coupled_odes,parms=pp,atol=a,rtol=a*100,method='lsoda',maxsteps=100000))
  disc<-abs(as.matrix(x[,sp])-as.matrix(xr[,sp]))/(1+abs(as.matrix(xr[,sp])))
  ix<-which(disc==max(disc),arr.ind=TRUE)[1,]
  cat('atol',a,'ETP',compute_tga_metrics(x)$ETP_nM_min,'HC',compute_hemostatic_metrics(x)$HC,
    'RHC',compute_hemostatic_metrics(xr)$HC,'disc',max(disc),'time',x$time[ix[1]],'species',sp[ix[2]],'\n')
  v<-as.double(unlist(x[400,sp]));names(v)<-sp
  dc<-.C('derivs',as.integer(40),as.double(200),as.double(v),dy=double(40),double(0),as.integer(0),PACKAGE='qsp_core')$dy
  dr<-coupled_odes(200,v,pp)[[1]]
  cat('DERIVATIVE ERROR',max(abs(dc-dr)),'species',sp[which.max(abs(dc-dr))],'\n')
}
