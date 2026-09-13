# Synthetic design distributions; not a calibrated NS virtual population.
lognormal_mean_one <- function(cv) {
  sdlog <- sqrt(log1p(cv^2))
  exp(rnorm(1, -sdlog^2/2, sdlog))
}
sample_virtual_patient <- function(ckd_stage="Normal", egfr=100, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  pk_cv <- c(CL=.34,V=.27,ka=.50,fu=.12)
  coag_cv <- c(VII=.15,VIII=.25,V=.15,II=.10,ATIII=.12,TFPI=.15)
  plt_cv <- c(k_thr=.35,k_agg=.35,k_txa_syn=.30,k_rel_adp=.30)
  list(body_weight=max(45,min(120,rnorm(1,75,15))),
       age=max(18,min(85,rnorm(1,50,15))),ckd_stage=ckd_stage,egfr=egfr,
       pk_eta=as.list(sapply(pk_cv,lognormal_mean_one)),
       coag_eta=sapply(coag_cv,lognormal_mean_one),
       plt_eta=sapply(plt_cv,lognormal_mean_one))
}
reference_patient <- function(stage="Normal", egfr=100, body_weight=70, age=40) {
  list(body_weight=body_weight,age=age,ckd_stage=stage,egfr=egfr,
       pk_eta=list(CL=1,V=1,ka=1,fu=1),coag_eta=numeric(),plt_eta=numeric())
}
prepare_patient <- function(patient, TF_pM=5, pk_override=NULL, rates=NULL, coupling=NULL,
                            t_end=1200) {
  p <- get_pbpk_params(patient$egfr, patient$body_weight, patient$ckd_stage, patient$age)
  if (!is.null(pk_override)) p <- modifyList(p,pk_override)
  p$CLint_renal <- p$CLint_renal*patient$pk_eta$CL
  p$CLint_hepatic <- p$CLint_hepatic*patient$pk_eta$CL
  p$fu_plasma <- min(.99,p$fu_plasma*patient$pk_eta$fu)
  p$ka <- p$ka*patient$pk_eta$ka
  for (nm in c("V_liver","V_kidney","V_plasma","V_peripheral")) p[[nm]] <- p[[nm]]*patient$pk_eta$V
  p <- recompute_derived_params(p)
  pk <- simulate_pbpk(1,tau=12,n_doses=28,params=p)
  ic <- apply_ckd_modifiers(get_coag_initial_conditions(TF_pM),patient$ckd_stage)
  ic <- apply_age_coag_modifiers(ic,patient$age)
  for (nm in intersect(names(patient$coag_eta),names(ic))) ic[nm] <- ic[nm]*patient$coag_eta[nm]
  pp <- apply_ckd_platelet_modifiers(get_platelet_params(),patient$ckd_stage)
  pp <- apply_age_platelet_modifiers(pp,patient$age)
  for(nm in intersect(names(patient$plt_eta),names(pp))) pp[[nm]] <- pp[[nm]]*patient$plt_eta[nm]
  list(patient=patient,params=p,pk_unit=pk,pk_metrics_unit=compute_pk_metrics(pk),
       cp_unit=extract_trough_nM(pk),ic=ic,plt=pp,
       rates=if(is.null(rates)) get_coag_rate_constants() else rates,
       coupling=if(is.null(coupling)) get_coupling_params() else coupling,
       t_end=t_end,TF_pM=TF_pM,cache=new.env(parent=emptyenv()))
}
evaluate_prepared <- function(prep,dose_mg,coupled=FALSE,compute_ti=TRUE,return_series=FALSE,
                              exposure="trough") {
  cp <- switch(exposure,trough=prep$cp_unit,
      peak=prep$pk_metrics_unit$Cmax_free*1000/prep$params$MW,
      mean=prep$pk_metrics_unit$AUCtau_free/12*1000/prep$params$MW)
  out <- simulate_fast(prep$ic,prep$plt,prep$rates,prep$coupling,cp*dose_mg,
                       coupled=coupled,t_end=prep$t_end)
  tg <- compute_tga_metrics(out)
  hc <- compute_hemostatic_metrics(out)
  ti <- NULL
  if(compute_ti) {
    key <- paste(coupled,prep$t_end)
    if(!exists(key,prep$cache,inherits=FALSE))
      assign(key,evaluate_prepared(prep,0,coupled,FALSE),prep$cache)
    base <- get(key,prep$cache,inherits=FALSE)
    ti <- compute_therapeutic_index(tg,base$tga_metrics,hc,base$hemostatic_metrics)
  }
  pk <- prep$pk_metrics_unit
  for(nm in setdiff(names(pk),"Tmax")) pk[[nm]] <- pk[[nm]]*dose_mg
  ans <- list(pk_metrics=pk,tga_metrics=tg,hemostatic_metrics=hc,therapeutic_index=ti,
              Cp_free_nM=cp*dose_mg,patient=prep$patient,coupled=coupled)
  if(return_series) ans$series <- out
  ans
}
# Grid brackets every local maximum, includes both endpoints and the no-drug comparator.
optimize_prepared <- function(prep,coupled=FALSE,interval=c(0,20),exposure="trough") {
  grid <- sort(unique(c(interval,0,exp(seq(log(.01),log(interval[2]),length.out=35)))))
  grid <- grid[grid>=interval[1] & grid<=interval[2]]
  objective <- function(d) evaluate_prepared(prep,d,coupled,exposure=exposure)$therapeutic_index$TI
  vals <- vapply(grid,objective,numeric(1))
  candidates <- data.frame(dose=grid,score=vals)
  for(i in 2:(length(grid)-1)) {
    if(vals[i]>=vals[i-1] && vals[i]>=vals[i+1] &&
       (vals[i]>vals[i-1] || vals[i]>vals[i+1])) {
      o <- optimize(objective,c(grid[i-1],grid[i+1]),maximum=TRUE,tol=.001)
      candidates <- rbind(candidates,data.frame(dose=o$maximum,score=o$objective))
    }
  }
  best <- candidates[which.max(candidates$score),]
  v <- evaluate_prepared(prep,best$dose,coupled,exposure=exposure)
  list(optimal_dose=best$dose,max_TI=best$score,
       efficacy_score=v$therapeutic_index$efficacy_score,
       safety_score=v$therapeutic_index$safety_score,
       peak_suppression=v$therapeutic_index$efficacy_score_peak,
       ETP_suppression=v$therapeutic_index$efficacy_score_etp,
       TI_etp=v$therapeutic_index$TI_etp,
       ETP_nM_min=v$tga_metrics$ETP_nM_min,HC=v$hemostatic_metrics$HC,
       boundary=abs(best$dose-interval[1])<.002 || abs(best$dose-interval[2])<.002,
       grid=data.frame(dose=grid,score=vals))
}
run_patient_simulation <- function(patient,dose_mg,TF_pM=5,compute_ti=TRUE) {
  evaluate_prepared(prepare_patient(patient,TF_pM),dose_mg,FALSE,compute_ti)
}
run_coupled_patient_simulation <- function(patient,dose_mg,TF_pM=5,compute_ti=TRUE) {
  evaluate_prepared(prepare_patient(patient,TF_pM),dose_mg,TRUE,compute_ti)
}
find_patient_optimal_dose <- function(patient,TF_pM=5,dose_interval=c(0,20)) {
  optimize_prepared(prepare_patient(patient,TF_pM),FALSE,dose_interval)
}
