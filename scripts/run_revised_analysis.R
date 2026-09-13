# Authoritative reproducible NS analysis. No clinical dose recommendations.
source('scripts/load_project.R')
dir.create('output/revised',recursive=TRUE,showWarnings=FALSE)
outdir <- 'output/revised'
save_csv <- function(x,name) write.csv(x,file.path(outdir,paste0(name,'.csv')),row.names=FALSE)
stages <- c('Normal','NS_Mild','NS_Moderate','NS_Severe')
set.seed(20260909)
cat('Reference scenarios\n')
scenarios <- list(); profiles <- list(); qprofiles <- list(); opts <- list(); curves <- list()
for(st in stages) {
  prep <- prepare_patient(reference_patient(st,100))
  ss <- extract_steady_state(prep$pk_unit)
  profiles[[st]] <- data.frame(stage=st,time_h=ss$time_in_interval,
    total_ng_mL=ss$C_total_ng_mL*5,free_ng_mL=ss$C_free_ng_mL*5)
  for(mode in c(FALSE,TRUE)) {
    label <- if(mode) 'coupled' else 'feedforward'
    for(d in c(0,2.5,5,10)) {
      x <- evaluate_prepared(prep,d,mode,return_series=TRUE)
      key <- paste(st,label,d)
      scenarios[[key]] <- data.frame(stage=st,mode=label,dose_mg_BID=d,
        as.data.frame(x$pk_metrics),ETP_nM_min=x$tga_metrics$ETP_nM_min,
        peak_thrombin_nM=x$tga_metrics$peak_thrombin,lag_s=x$tga_metrics$lag_time,
        peak_time_s=x$tga_metrics$time_to_peak,HC_model_units_s=x$hemostatic_metrics$HC,
        ETP_suppression=x$therapeutic_index$efficacy_score_etp,
        peak_suppression=x$therapeutic_index$efficacy_score_peak,
        HC_retention=x$therapeutic_index$safety_score,score=x$therapeutic_index$TI,
        score_etp=x$therapeutic_index$TI_etp)
      if(d %in% c(0,5)) qprofiles[[key]] <- data.frame(stage=st,mode=label,dose=d,
        time_s=x$series$time,IIa=x$series$IIa,PLT_agg=x$series$PLT_agg)
    }
    o <- optimize_prepared(prep,mode)
    opts[[paste(st,label)]] <- data.frame(stage=st,mode=label,
      as.data.frame(o[setdiff(names(o),'grid')]),weak_objective=o$max_TI<.001)
    curves[[paste(st,label)]] <- data.frame(stage=st,mode=label,o$grid)
  }
}
scenarios <- do.call(rbind,scenarios); optima <- do.call(rbind,opts)
save_csv(scenarios,'reference_scenarios');save_csv(optima,'surrogate_optima')
save_csv(do.call(rbind,profiles),'pk_profiles');save_csv(do.call(rbind,qprofiles),'qsp_profiles')
save_csv(do.call(rbind,curves),'score_curves')
cat('External comparisons\n')
ref <- read.csv('data/external_reference.csv')
ref$prediction <- NA_real_
for(i in seq_len(nrow(ref))) {
  st <- if(ref$cohort[i]=='NS') 'NS_Severe' else 'Normal'
  egfr <- if(ref$study[i]=='Frost2013') 120 else if(st=='Normal') 107 else 81
  age <- if(ref$study[i]=='Frost2013') 40 else 51
  x <- evaluate_prepared(prepare_patient(reference_patient(st,egfr,70,age)),5)
  vals <- c(x$pk_metrics,list(ETP_nM_min=x$tga_metrics$ETP_nM_min))
  ref$prediction[i] <- vals[[ref$metric[i]]]
}
ref$predicted_observed_ratio <- ref$prediction/ref$value
ref$relative_error_pct <- 100*(ref$prediction/ref$value-1)
save_csv(ref,'external_comparison')
cat('PK mechanism decomposition\n')
pkdec <- list()
for(st in stages) {
  pt <- reference_patient(st,100); p <- get_pbpk_params(100,70,st,40)
  for(h in c(1,.82,.70)) {
    for(loss in c(0,p$CL_proteinuria)) {
      pr <- prepare_patient(pt,pk_override=list(CLint_hepatic=p$CLint_hepatic*h,CL_proteinuria=loss))
      x <- evaluate_prepared(pr,5,FALSE,FALSE)
      pkdec[[length(pkdec)+1]] <- data.frame(stage=st,nonrenal_CLint_multiplier=h,
        bound_loss_L_h=loss,as.data.frame(x$pk_metrics))
    }
  }
}
save_csv(unique(do.call(rbind,pkdec)),'pk_mechanism_scenarios')
cat('Structural and endpoint sensitivity\n')
structure_rows <- list()
for(st in stages) for(tf in c(1,5,25)) for(pc in c(0,6.2e-6)) for(duration in c(600,1200,2400)) {
  rates <- get_coag_rate_constants();rates$k39 <- pc
  pr <- prepare_patient(reference_patient(st,100),tf,rates=rates,t_end=duration)
  for(ex in c('trough','peak','mean')) {
    x <- evaluate_prepared(pr,5,exposure=ex)
    structure_rows[[length(structure_rows)+1]] <- data.frame(stage=st,TF_pM=tf,
      k39_nM_s=pc,duration_s=duration,exposure=ex,ETP=x$tga_metrics$ETP_nM_min,
      peak=x$tga_metrics$peak_thrombin,
      ETP_suppression=x$therapeutic_index$efficacy_score_etp,
      peak_suppression=x$therapeutic_index$efficacy_score_peak,
      HC_retention=x$therapeutic_index$safety_score,score=x$therapeutic_index$TI,
      score_etp=x$therapeutic_index$TI_etp)
  }
}
save_csv(do.call(rbind,structure_rows),'structural_sensitivity')
cat('Local parameter sensitivity\n')
sa <- list()
for(st in c('Normal','NS_Severe')) {
  pr <- prepare_patient(reference_patient(st,100)); base <- evaluate_prepared(pr,5)
  for(nm in c('fu_plasma','CLint_hepatic','CLint_renal','ka','k23','k28','k30','k32','k33','k_thr','k_agg')) for(fac in c(.5,1.5)) {
    alt <- pr
    if(nm %in% names(pr$params)) {
      override <- setNames(list(pr$params[[nm]]*fac),nm)
      alt <- prepare_patient(reference_patient(st,100),pk_override=override)
    } else if(nm %in% names(pr$rates)) alt$rates[[nm]] <- alt$rates[[nm]]*fac
    else alt$plt[[nm]] <- alt$plt[[nm]]*fac
    alt$cache <- new.env(parent=emptyenv())
    x <- evaluate_prepared(alt,5)
    for(endpoint in c('ETP','HC')) {
      a <- if(endpoint=='ETP') x$tga_metrics$ETP else x$hemostatic_metrics$HC
      b <- if(endpoint=='ETP') base$tga_metrics$ETP else base$hemostatic_metrics$HC
      sa[[length(sa)+1]] <- data.frame(stage=st,parameter=nm,multiplier=fac,endpoint=endpoint,
        change_pct=100*(a/b-1),normalized_sensitivity=(a/b-1)/(fac-1))
    }
  }
}
save_csv(do.call(rbind,sa),'local_sensitivity')
cat('Matched synthetic population and grid\n')
source('scripts/run_population_analysis.R')
cat('TF prior propagation\n')
source('scripts/run_tf_propagation.R')
source('scripts/finalize_results.R')
cat('REVISED ANALYSIS COMPLETE\n')
