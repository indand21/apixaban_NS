source('scripts/load_project.R')
dir.create('output/revised',recursive=TRUE,showWarnings=FALSE)
checks <- list()
check <- function(name,ok,value=NA_real_,tolerance=NA_real_) {
  checks[[length(checks)+1]] <<- data.frame(check=name,pass=isTRUE(ok),value=value,tolerance=tolerance)
  cat(if(isTRUE(ok)) 'PASS' else 'FAIL',name,'value=',value,'\n')
}
p <- prepare_patient(reference_patient())
sp <- c(.coag_species_names,.plt_species_names)
for(mode in c(FALSE,TRUE)) {
  x <- simulate_fast(p$ic,p$plt,p$rates,p$coupling,10,coupled=mode)
  cp <- p$coupling
  if(!mode) {cp$alpha_assembly<-0;cp$alpha_catalysis<-0}
  xr <- simulate_coupled(p$ic,get_platelet_initial_conditions(p$plt),p$rates,p$plt,cp,10)$full_result
  discrepancy <- max(abs(as.matrix(x[,sp])-as.matrix(xr[,sp]))/(1+abs(as.matrix(xr[,sp]))))
  check(paste('R/C trajectory agreement',mode),discrepancy<1e-4,discrepancy,1e-4)
  check(paste('Nonnegative finite complete trajectory',mode),all(is.finite(as.matrix(x)))&&min(as.matrix(x[,sp]))>-1e-7&&tail(x$time,1)==1200,min(as.matrix(x[,sp])),-1e-7)
  totals <- list(TF=c('TF','TF_VII','TF_VIIa','TF_VIIa_X','TF_VIIa_Xa','TF_VIIa_IX','TFPI_Xa_TF_VIIa'),
    FII=c('II','IIa','mIIa','Va_Xa_II','ATIII_IIa'),
    FV=c('V','Va','Va_Xa','Va_Xa_II','Va_i'),
    FVII=c('VII','VIIa','TF_VII','TF_VIIa','TF_VIIa_X','TF_VIIa_Xa','TF_VIIa_IX','TFPI_Xa_TF_VIIa'),
    FVIII=c('VIII','VIIIa','IXa_VIIIa','IXa_VIIIa_X','VIIIa_i'),
    FIX=c('IX','IXa','TF_VIIa_IX','IXa_VIIIa','IXa_VIIIa_X','ATIII_IXa'),
    FX=c('X','Xa','TF_VIIa_X','TF_VIIa_Xa','IXa_VIIIa_X','Va_Xa','Va_Xa_II','TFPI_Xa','TFPI_Xa_TF_VIIa','ATIII_Xa'),
    AT=c('ATIII','ATIII_IIa','ATIII_Xa','ATIII_IXa'),PC=c('PC','APC'),TFPI=c('TFPI','TFPI_Xa','TFPI_Xa_TF_VIIa'))
  for(nm in names(totals)) {
    v<-rowSums(x[,totals[[nm]]]);err<-max(abs(v/v[1]-1))
    check(paste('Moiety conservation',nm,mode),err<1e-6,err,1e-6)
  }
  tight <- simulate_fast(p$ic,p$plt,p$rates,p$coupling,10,coupled=mode,dt=.25,rtol=1e-13,atol=1e-16)
  m1<-c(compute_tga_metrics(x)$ETP,compute_hemostatic_metrics(x)$HC)
  m2<-c(compute_tga_metrics(tight)$ETP,compute_hemostatic_metrics(tight)$HC)
  err<-max(abs(m1/m2-1));check(paste('Endpoint convergence',mode),err<1e-4,err,1e-4)
  platelet_total<-rowSums(x[,c('PLT_rest','PLT_act','PLT_agg')])
  sink<-p$plt$k_deact*compute_auc(x$time,x$PLT_act)
  err<-abs(tail(platelet_total,1)+sink-platelet_total[1])/platelet_total[1]
  check(paste('Platelet balance includes desensitized sink',mode),err<1e-5,err,1e-5)
}
pk2<-compute_pk_metrics(simulate_pbpk(2,n_doses=28,params=p$params))
check('Linear dose scaling',abs(pk2$AUCtau_total/(2*p$pk_metrics_unit$AUCtau_total)-1)<1e-5)
pk56<-compute_pk_metrics(simulate_pbpk(1,n_doses=56,params=p$params))
err<-abs(pk56$AUCtau_total/p$pk_metrics_unit$AUCtau_total-1)
check('28 versus 56 dose steady-state convergence',err<1e-5,err,1e-5)
v<-c(A_gut=1,A_liver=.3,A_kidney=.1,A_plasma=.8,A_peripheral=.5)
q<-p$params;dd<-pbpk_odes(1,v,q)[[1]]
loss<-q$CLint_hepatic*q$fu_plasma*v['A_liver']/q$V_liver/q$Kp_liver+
  (q$CLint_renal*q$fu_plasma+q$CL_proteinuria*(1-q$fu_plasma))*v['A_kidney']/q$V_kidney/q$Kp_kidney
check('PBPK instantaneous mass balance',abs(sum(dd)+loss)<1e-12,abs(sum(dd)+loss),1e-12)
pr<-get_pbpk_params(120)
check('Systemic reference clearance is not multiplied by F',abs(pr$CL_total-3.3)<1e-12&&abs(pr$CL_renal-.9)<1e-12)
check('Primary assay excludes implicit thrombomodulin',p$rates$k39==0)
check('Inactive rate constants removed',!any(c('k17','k42')%in%names(p$rates)))
check('Fixed TF across NS scenarios',all(vapply(c('Normal','NS_Mild','NS_Moderate','NS_Severe'),function(s)apply_ckd_modifiers(p$ic,s)['TF']==.005,logical(1))))
z<-evaluate_prepared(p,0)
check('Matched no-drug reference gives score zero and retention one',z$therapeutic_index$TI==0&&z$therapeutic_index$safety_score==1)
check('ETP seconds-to-minutes conversion',abs(z$tga_metrics$ETP/60-z$tga_metrics$ETP_nM_min)<1e-12)
pt<-reference_patient();pt$pk_eta$V<-1.2;pt$pk_eta$ka<-.8;pt$plt_eta<-c(k_thr=1.3,k_agg=.7)
p2<-prepare_patient(pt)
check('Sampled volume effect applied',abs(p2$params$V_plasma/p$params$V_plasma-1.2)<1e-12)
check('Sampled absorption effect applied',abs(p2$params$ka/p$params$ka-.8)<1e-12)
check('Sampled platelet effects applied',abs(p2$plt$k_thr/p$plt$k_thr-1.3)<1e-12&&abs(p2$plt$k_agg/p$plt$k_agg-.7)<1e-12)
check('Invalid ETP-to-stroke mapping disabled',inherits(try(etp_to_stroke_hazard(1000),silent=TRUE),'try-error'))
h1<-suppressWarnings(hc_to_bleeding_hazard(100,hc_healthy_nodrug=100))
h2<-suppressWarnings(hc_to_bleeding_hazard(50,hc_healthy_nodrug=100))
check('Toy hazard sign regression (not clinical validation)',abs(h2/h1-2)<1e-12)
set.seed(80);draws<-replicate(50000,lognormal_mean_one(.34))
check('Lognormal mean-one parameterization',abs(mean(draws)-1)<.005,mean(draws),.005)
check('Lognormal target CV',abs(sd(draws)/mean(draws)-.34)<.005,sd(draws)/mean(draws),.005)
flat<-p;flat$rates$k16<-0;flat$rates$k23<-0;flat$cache<-new.env(parent=emptyenv())
fo<-optimize_prepared(flat)
check('Exactly flat score chooses no drug without an invented optimum',fo$max_TI==0&&fo$optimal_dose==0&&fo$boundary)
check('Dose search explicitly includes both endpoints',any(fo$grid$dose==0)&&any(fo$grid$dose==20))
# --- Tier-1 revision checks (objective, PBPK structure, NS severity scale) ---
z5<-evaluate_prepared(p,5)
ti5<-z5$therapeutic_index
check('Reported efficacy endpoint is peak thrombin suppression',
  identical(ti5$objective,'peak')&&isTRUE(all.equal(ti5$efficacy_score,ti5$efficacy_score_peak)))
check('Legacy ETP-based score retained for audit',
  is.finite(ti5$efficacy_score_etp)&&is.finite(ti5$TI_etp))
# The ETP objective was flat (about 1e-4 at therapeutic exposure), which is what
# drove the dose search onto its search boundary. The peak endpoint must be
# materially more responsive or the fix has not worked.
check('Peak-thrombin objective is not degenerate at therapeutic exposure',
  ti5$efficacy_score_peak>0.05&&ti5$efficacy_score_peak>50*ti5$efficacy_score_etp,
  ti5$efficacy_score_peak,0.05)
# The absorption lag delays input; it must not create or destroy drug, so
# steady-state unbound exposure stays F x Dose / CLint regardless of the lag.
pf<-get_pbpk_params(120,70,'Normal',40)
auc_free<-function(tl){pp<-pf;pp$Tlag<-tl
  compute_pk_metrics(simulate_pbpk(5,n_doses=28,params=recompute_derived_params(pp)))$AUCtau_free}
d2<-abs(auc_free(pf$Tlag)/auc_free(0)-1)
check('Absorption lag does not change steady-state unbound exposure',d2<1e-4,d2,1e-4)
check('Absorption lag is calibrated to a nonzero delay',pf$Tlag>0.1,pf$Tlag)
raw<-read.csv('data/ns_factor_levels.csv',stringsAsFactors=FALSE)
m1<-get_ckd_modifiers('NS_Severe',severity_scale=1)
check('Severity scale of one reproduces the raw factor table',
  max(abs(m1-setNames(raw$NS_Severe,raw$species)))<1e-12)
m0<-get_ckd_modifiers('NS_Severe',severity_scale=0)
check('Severity scale of zero removes the NS coagulation gradient',max(abs(m0-1))<1e-12)
check('Severity scale is calibrated, not left at the raw table',
  abs(get_ns_severity_scale()-1)>1e-6,get_ns_severity_scale())
# Calibrated PK must reproduce the Frost 2013 aggregate moments it was fitted to.
er<-read.csv('data/external_reference.csv',stringsAsFactors=FALSE)
fr<-setNames(er$value[er$study=='Frost2013'],er$metric[er$study=='Frost2013'])
mf<-compute_pk_metrics(simulate_pbpk(5,n_doses=28,
  params=get_pbpk_params(120,70,'Normal',40)))
eauc<-abs(mf$AUCtau_total/fr['AUCtau_total']-1)
check('Calibrated AUCtau matches the Frost observation',eauc<0.005,eauc,0.005)
ecmax<-abs(mf$Cmax_total/fr['Cmax_total']-1)
check('Calibrated Cmax within 12 percent of observation',ecmax<0.12,ecmax,0.12)
ectr<-abs(mf$Ctrough_total/fr['Ctrough_total']-1)
check('Calibrated Ctrough within 12 percent of observation',ectr<0.12,ectr,0.12)
# Tmax is a median of n=6 with a reported range of 2-4 h. The inherited model
# gave 1.4 h, below that range; the calibrated model must at least reach it.
check('Calibrated Tmax reaches the observed range',mf$Tmax>=2&&mf$Tmax<=4,mf$Tmax)
# Calibrated severity scale must reproduce the observed NS:healthy ETP ratio.
kv<-er[er$study=='Kelddal2025'&er$metric=='ETP_nM_min',]
tr<-kv$value[kv$cohort=='NS']/kv$value[kv$cohort=='Healthy']
etp_at<-function(st,eg) evaluate_prepared(prepare_patient(reference_patient(st,eg,70,51)),5,
  compute_ti=FALSE)$tga_metrics$ETP_nM_min
er_ratio<-abs(etp_at('NS_Severe',81)/etp_at('Normal',107)/tr-1)
check('Calibrated NS severity reproduces the observed ETP ratio',er_ratio<0.01,er_ratio,0.01)
csv_files<-list.files('data',pattern='[.]csv$',full.names=TRUE)
csv_ok<-vapply(csv_files,function(f){n<-count.fields(f,sep=',',quote='\"',comment.char='');n<-n[n>0];length(unique(n))==1},logical(1))
check('All parameter and reference CSV rows conform to their headers',all(csv_ok))
res<-do.call(rbind,checks);write.csv(res,'output/revised/revision_tests.csv',row.names=FALSE)
cat(sprintf('REVISION TESTS: %d passed, %d failed\n',sum(res$pass),sum(!res$pass)))
if(any(!res$pass)) quit(status=1)
