# Shared component implementation prevents coupled/sequential equation drift.
get_coupling_params <- function(data_dir = "data/") {
  df <- read.csv(file.path(data_dir, "coupling_parameters.csv"))
  as.list(setNames(df$value, df$parameter))
}
compute_ps_enhancement <- function(PS_exp, alpha, Km_PS, n_ps) {
  ps <- max(0, PS_exp)^n_ps
  1 + alpha * ps / (Km_PS^n_ps + ps)
}
coupled_odes <- function(t, state, params) {
  asm <- compute_ps_enhancement(state["PS_exp"], params$alpha_assembly, params$Km_PS, params$n_ps)
  cat <- compute_ps_enhancement(state["PS_exp"], params$alpha_catalysis, params$Km_PS, params$n_ps)
  cp <- params
  cp$k19 <- params$k19 * asm
  cp$k23 <- params$k23 * cat
  coag <- coag_odes(t, state[1:34], cp)
  pp <- params
  pp$IIa_forcing <- max(0, state["IIa"])
  plt <- platelet_odes(t, state[35:40], pp)
  c(list(c(coag[[1]], plt[[1]])), coag[-1],
    list(ps_enhance_asm = unname(asm), ps_enhance_cat = unname(cat),
         IIa_input = unname(state["IIa"])))
}
