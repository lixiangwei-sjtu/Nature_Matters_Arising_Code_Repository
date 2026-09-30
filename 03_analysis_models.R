# Survey-model specifications for the Nature Matters Arising analysis.
# No participant-level data are distributed with this repository.

suppressPackageStartupMessages({
  library(survey)
  library(mitools)
})

zscore <- function(x) as.numeric(scale(x))

# -----------------------------------------------------------------------------
# NHANES survey-weighted model
# -----------------------------------------------------------------------------
# Expected analytic columns (rename as needed in your authorized/public dataset):
# outcome, exposure_z, age, sex, race_ethnicity, education, smoking, alcohol,
# physical_activity, marital_status, socioeconomic_position, total_energy,
# dnam_weight, strata, psu

fit_nhanes_model <- function(dat, outcome, exposure) {
  d <- subset(dat, !is.na(dat[[outcome]]) & !is.na(dat[[exposure]]) & dnam_weight > 0)
  des <- svydesign(
    ids = ~psu,
    strata = ~strata,
    weights = ~dnam_weight,
    nest = TRUE,
    data = d
  )
  f <- as.formula(paste0(
    outcome, " ~ ", exposure,
    " + age + sex + race_ethnicity + education + smoking + alcohol",
    " + physical_activity + marital_status + socioeconomic_position + total_energy"
  ))
  svyglm(f, design = des)
}

# -----------------------------------------------------------------------------
# HRS survey-weighted model (single completed dataset)
# -----------------------------------------------------------------------------
# Replace hrs_weight/hrs_strata/hrs_psu with the exact design variables used in
# the locked analysis.

fit_hrs_model <- function(dat, outcome, exposure) {
  d <- subset(dat, !is.na(dat[[outcome]]) & !is.na(dat[[exposure]]))
  des <- svydesign(
    ids = ~hrs_psu,
    strata = ~hrs_strata,
    weights = ~hrs_weight,
    nest = TRUE,
    data = d
  )
  f <- as.formula(paste0(
    outcome, " ~ ", exposure,
    " + age + sex + race_ethnicity + education + smoking + alcohol",
    " + physical_activity + marital_status + socioeconomic_position + total_energy"
  ))
  svyglm(f, design = des)
}

# -----------------------------------------------------------------------------
# Rubin pooling across completed HRS datasets
# -----------------------------------------------------------------------------
# `completed_data` should be a list of completed participant-level datasets.
# This generic function fits the same model to each completed dataset and pools
# coefficients/covariances with MIcombine.

pool_hrs_mi <- function(completed_data, outcome, exposure) {
  fits <- lapply(completed_data, function(d) fit_hrs_model(d, outcome, exposure))
  MIcombine(fits)
}
