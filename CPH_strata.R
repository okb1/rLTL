library(survival)

# File config.csv must be in the current directory
config_variables <- read.csv("config_CPH_strata.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

dir_in <- config_variables["dir_in", "Value"]
dir_out <- config_variables["dir_out", "Value"]
file_name <- config_variables["file_name", "Value"]
age_shift <- as.numeric(config_variables["age_shift", "Value"])
str_cause <- config_variables["str_cause", "Value"]
name_rLTL <- config_variables["name_rLTL", "Value"]
name_event <- config_variables["name_event", "Value"]
name_sex <- config_variables["name_sex", "Value"]
name_age0 <- config_variables["name_age0", "Value"]
name_time <- config_variables["name_time", "Value"]
exclude_sex <- as.numeric(config_variables["exclude_sex", "Value"])
age0gr_name <- config_variables["age0gr_name", "Value"]

str_strata <- sprintf("+ strata(%s)", age0gr_name)
ages_S <- seq(50, 90, 10)
num_ages_S <- length(ages_S)
value_rLTL0 <- 0
value_rLTL1 <- -1
model_abbr <- "age0m_rLTL0"
str_covars <- "age0m + rLTL0 "
study_name <- basename(file_name)

outputcox <- function(myCPH, study_name, modelName, out_fname, out_fname_zph = "") {
  sum_myCPH <- summary(myCPH)
  zph_myCPH <- cox.zph(myCPH)
  d2 <- data.frame(sum_myCPH$coefficients, rbind(sum_myCPH$conf.int[, 3:4]), rbind(zph_myCPH$table[1:max(1, nrow(zph_myCPH$table) - 1), ]))
  d2$Cause <- str_cause
  d2$parameter <- rownames(d2)
  names(d2)[which(names(d2) == "coef")] <- paste("ParameterEstimate")
  names(d2)[which(names(d2) == "se.coef.")] <- paste("StandardErrors")
  names(d2)[which(names(d2) == "Pr...z..")] <- paste("p-value")
  names(d2)[which(names(d2) == "exp.coef.")] <- paste("HazardRatio")
  names(d2)[which(names(d2) == "lower..95")] <- paste("HR_low95")
  names(d2)[which(names(d2) == "upper..95")] <- paste("HR_up95")
  names(d2)[which(names(d2) == "p")] <- paste("p-valuePropHazards")
  d2$Study <- study_name
  d2$model <- modelName
  d2$R2 <- sum_myCPH$rsq[1]
  d2$R2_nagelkerke <- sum_myCPH$rsq[1] / sum_myCPH$rsq[2]
  d2$HazardRatioMinus <- 1 / d2$HazardRatio
  d2$HR_low95_minus <- 1 / d2$HR_up95
  d2$HR_up95_minus <- 1 / d2$HR_low95

  dout <- d2[c(
    "Cause", "Study", "model", "parameter", "ParameterEstimate", "StandardErrors", "p-value",
    "HazardRatio", "HR_low95", "HR_up95", "p-valuePropHazards", "R2", "R2_nagelkerke",
    "HazardRatioMinus", "HR_low95_minus", "HR_up95_minus"
  )]

  write.csv(dout, file = out_fname, row.names = FALSE, na = "")

  if (out_fname_zph != "") {
    test_ph <- as.data.frame(zph_myCPH$table, row.names = row.names(zph_myCPH$table))
    write.csv(test_ph, file = out_fname_zph, row.names = TRUE, na = "")
  }
}

data_study <- read.csv(sprintf("%s%s.csv", dir_in, file_name), header = TRUE, stringsAsFactors = FALSE)
data_study$age0m <- data_study[, name_age0] - age_shift
data_study$rLTL0 <- data_study[, name_rLTL]
data_study$sexM <- ifelse(data_study[, name_sex] == 2, 0, data_study[, name_sex])
max_strata <- max(data_study[age0gr_name])
min_strata <- min(data_study[age0gr_name])
num_strata <- max_strata - min_strata + 1
strata_levels <- seq(min_strata, max_strata)

if (exclude_sex == 1) {
  str_sex <- ""
} else {
  str_sex <- "+ sexM"
}

dir_out2 <- dirname(sprintf("%s%s/", dir_out, file_name))
dir.create(dir_out2, recursive = TRUE, showWarnings = FALSE)

model_i <- c(sprintf("Surv(%s,%s) ~ %s %s %s ", name_time, name_event, str_covars, str_sex, str_strata))
myCPH_i <- coxph(as.formula(model_i), ties = "efron", data = data_study)

out_fname <- sprintf("%s/%s_cox_time_s.csv", dir_out2, basename(file_name))
out_fname_zph <- sprintf("%s/%s_cox_time_zph_s.csv", dir_out2, basename(file_name))

outputcox(myCPH_i, study_name, model_abbr, out_fname, out_fname_zph)

survfit_S_data <- data.frame()
if (exclude_sex == 1) {
  for (i_rLTL0 in value_rLTL0:value_rLTL1) {
    for (i_age in 1:num_ages_S) {
      newdata_cph_S_i <- data.frame(age0m = rep(ages_S[i_age] - age_shift, num_strata), rLTL0 = rep(i_rLTL0, num_strata))
      newdata_cph_S_i[age0gr_name] <- seq(min_strata, max_strata)

      survfit_S_i <- survfit(myCPH_i, newdata = newdata_cph_S_i, censor = FALSE)
      survfit_S_data_i <- data.frame(
        age0m = ages_S[i_age] - age_shift, rLTL0 = i_rLTL0, sexM = 0, time = survfit_S_i$time,
        n_risk = survfit_S_i$n.risk, n_event = survfit_S_i$n.event, n_censor = survfit_S_i$n.censor,
        surv = survfit_S_i$surv, cumhaz = survfit_S_i$cumhaz, std_err = survfit_S_i$std.err,
        lower = survfit_S_i$lower, upper = survfit_S_i$upper
      )

      end_row_strata <- cumsum(survfit_S_i$strata)
      start_row_strata <- end_row_strata - c(survfit_S_i$strata[1], survfit_S_i$strata[-1]) + 1

      strata_S_i <- rep(NA, length(survfit_S_data_i$time))
      for (i_strata in 1:num_strata) {
        strata_S_i[start_row_strata[i_strata]:end_row_strata[i_strata]] <- strata_levels[i_strata]
      }

      survfit_S_data_i[age0gr_name] <- strata_S_i

      survfit_S_data <- rbind(survfit_S_data, survfit_S_data_i)
    }
  }
} else {
  for (i_sex in 0:1) {
    for (i_rLTL0 in value_rLTL0:value_rLTL1) {
      for (i_age in 1:num_ages_S) {
        newdata_cph_S_i <- data.frame(
          age0m = rep(ages_S[i_age] - age_shift, num_strata), rLTL0 = rep(i_rLTL0, num_strata),
          sexM = rep(i_sex, num_strata)
        )
        newdata_cph_S_i[age0gr_name] <- seq(min_strata, max_strata)

        survfit_S_i <- survfit(myCPH_i, newdata = newdata_cph_S_i, censor = FALSE)
        survfit_S_data_i <- data.frame(
          age0m = ages_S[i_age] - age_shift, rLTL0 = i_rLTL0, sexM = i_sex, time = survfit_S_i$time,
          n_risk = survfit_S_i$n.risk, n_event = survfit_S_i$n.event, n_censor = survfit_S_i$n.censor,
          surv = survfit_S_i$surv, cumhaz = survfit_S_i$cumhaz, std_err = survfit_S_i$std.err,
          lower = survfit_S_i$lower, upper = survfit_S_i$upper
        )

        end_row_strata <- cumsum(survfit_S_i$strata)
        start_row_strata <- end_row_strata - c(survfit_S_i$strata[1], survfit_S_i$strata[-1]) + 1

        strata_S_i <- rep(NA, length(survfit_S_data_i$time))
        for (i_strata in 1:num_strata) {
          strata_S_i[start_row_strata[i_strata]:end_row_strata[i_strata]] <- strata_levels[i_strata]
        }

        survfit_S_data_i[age0gr_name] <- strata_S_i

        survfit_S_data <- rbind(survfit_S_data, survfit_S_data_i)
      }
    }
  }
}

write.csv(survfit_S_data, sprintf("%s/%s_S_s.csv", dir_out2, basename(file_name)), row.names = FALSE, na = "")

var_betas <- as.data.frame(myCPH_i$var)
names(var_betas) <- names(myCPH_i$coefficients)
row.names(var_betas) <- names(var_betas)

out_fname_var <- sprintf("%s/%s_cox_time_var_s.csv", dir_out2, basename(file_name))
write.csv(var_betas, file = out_fname_var, row.names = TRUE, na = "")
