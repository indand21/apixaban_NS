# Compiled constant-exposure solver; equation identity is tested against R.
.fast_parameter_names <- c("k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11","k12","k13","k14","k15","k16","k18","k19","k20","k21","k22","k23","k24","k25","k26","k27","k28","k29","k30","k31","k32","k33","k34","k35","k36","k37","k38","k39","k40","k41","k_thr","Km_thr","n_hill","k_adp","Km_adp","k_txa","Km_txa","k_agg","Km_fbg","FBG","k_rel_adp","k_txa_syn","k_ps","k_adp_deg","k_txa_deg","k_deact","PLT_rest_0","Km_PS","n_ps","alpha_assembly","alpha_catalysis","Cp_free_nM","Ki_apixaban")
load_fast_solver <- function() {
  if ("qsp_core" %in% names(getLoadedDLLs())) return(invisible(NULL))
  ext <- .Platform$dynlib.ext
  path <- file.path("src", paste0("qsp_core", ext))
  if (!file.exists(path)) stop("Compile src/qsp_core.c with R CMD SHLIB before analysis.")
  if (!"qsp_core" %in% names(getLoadedDLLs())) dyn.load(normalizePath(path))
}
simulate_fast <- function(ic_coag, plt_params, rate_constants = get_coag_rate_constants(),
                          coupling_params = get_coupling_params(), Cp_free_nM = 0,
                          coupled = FALSE, t_end = 1200, dt = 0.5,
                          Ki_apixaban = 0.25, rtol = 1e-12, atol = 1e-14) {
  load_fast_solver()
  if (!coupled) { coupling_params$alpha_assembly <- 0; coupling_params$alpha_catalysis <- 0 }
  p <- c(rate_constants, plt_params, coupling_params,
         list(Cp_free_nM=Cp_free_nM, Ki_apixaban=Ki_apixaban))
  vals <- unlist(p[.fast_parameter_names], use.names=FALSE)
  stopifnot(length(vals)==length(.fast_parameter_names), all(is.finite(vals)), Cp_free_nM>=0)
  ic <- c(ic_coag, get_platelet_initial_conditions(plt_params))
  out <- deSolve::ode(y=ic, times=unique(c(seq(0,t_end,by=dt),t_end)),
    func="derivs", parms=vals, dllname="qsp_core", initfunc="initmod",
    method="lsoda", rtol=rtol, atol=atol, maxsteps=100000)
  if (tail(out[,1],1) < t_end || any(!is.finite(out))) stop("Incomplete/nonfinite ODE solution")
  as.data.frame(out)
}
