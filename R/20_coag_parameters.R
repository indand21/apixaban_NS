# CSV is the authoritative reduced-cascade parameter source.
get_coag_rate_constants <- function(data_dir = "data/") {
  df <- read.csv(file.path(data_dir, "hockin_mann_parameters.csv"))
  stopifnot(!anyDuplicated(df$param_id), all(is.finite(df$value)), all(df$value >= 0))
  values <- ifelse(df$units == "1/(M*s)", df$value / 1e9, df$value)
  as.list(setNames(values, df$param_id))
}
get_coag_initial_conditions <- function(TF_pM = 5, data_dir = "data/") {
  df <- read.csv(file.path(data_dir, "hockin_mann_species.csv"))
  df <- df[order(df$species_id), ]
  ic <- setNames(df$initial_nM, df$species_name)
  ic["TF"] <- TF_pM / 1000
  stopifnot(length(ic) == 34, all(ic >= 0))
  ic
}
