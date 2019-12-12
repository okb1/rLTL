library(splines)
library(survival)

# File config.csv must be in the current directory
config_variables <- read.csv("config_CPH_splines.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

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
name_aget <- config_variables["name_aget", "Value"]
n_knots <- as.numeric(config_variables["n_knots", "Value"])

if (n_knots == 3) {
  str_splines <- "BA1 + BA2 + BA3 "
} else if (n_knots == 5) {
  str_splines <- "BA1 + BA2 + BA3 + BA4 + BA5 "
} else {
  stop("n_knots must equal 3 or 5")
}


ages_spline_out_all <- seq(50, 90, 0.01)
ages_S <- seq(50, 90, 10)
num_ages_S <- length(ages_S)
value_rLTL0 <- 0
value_rLTL1 <- -1
study_name <- basename(file_name)
model_abbr <- "age0m_rLTL0"
str_covars <- sprintf(sprintf("%s+ rLTL0 ", str_splines))

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

data_study_died <- data_study[data_study[, name_event] == 1, ]

min_aget <- min(data_study[, name_aget])
min_aget

max_aget <- max(data_study[, name_aget])
max_aget

quant_aget_5 <- round(quantile(data_study[, name_aget], probs = c(0.2, 0.4, 0.6, 0.8)))
quant_aget_5

quant_aget_died_5 <- round(quantile(data_study_died[, name_aget], probs = c(0.2, 0.4, 0.6, 0.8)))
quant_aget_died_5

quant_aget_3 <- round(quantile(data_study[, name_aget], probs = c(1 / 3, 2 / 3)))
quant_aget_3

quant_aget_died_3 <- round(quantile(data_study_died[, name_aget], probs = c(1 / 3, 2 / 3)))
quant_aget_died_3

my_spline_basis_5 <- splines::ns(data_study[, name_age0], knots = quant_aget_died_5, Boundary.knots = c(min_aget, max_aget))
my_spline_basis_3 <- splines::ns(data_study[, name_age0], knots = quant_aget_died_3, Boundary.knots = c(min_aget, max_aget))

if (n_knots == 3) {
  my_spline_basis <- my_spline_basis_3
  quant_aget <- quant_aget_3
  data_study$BA1 <- my_spline_basis[, 1]
  data_study$BA2 <- my_spline_basis[, 2]
  data_study$BA3 <- my_spline_basis[, 3]
}
if (n_knots == 5) {
  my_spline_basis <- my_spline_basis_5
  quant_aget <- quant_aget_5
  data_study$BA1 <- my_spline_basis[, 1]
  data_study$BA2 <- my_spline_basis[, 2]
  data_study$BA3 <- my_spline_basis[, 3]
  data_study$BA4 <- my_spline_basis[, 4]
  data_study$BA5 <- my_spline_basis[, 5]
}

my_spline_basis_S <- splines::ns(ages_S, knots = quant_aget, Boundary.knots = c(min_aget, max_aget))
row.names(my_spline_basis_S) <- ages_S

if (exclude_sex == 1) {
  str_sex <- ""
} else {
  str_sex <- "+ sexM"
}

dir_out2 <- dirname(sprintf("%s%s/", dir_out, file_name))
dir.create(dir_out2, recursive = TRUE, showWarnings = FALSE)

model_i <- c(sprintf("Surv(%s, %s) ~ %s %s ", name_time, name_event, str_covars, str_sex))
myCPH_i <- coxph(as.formula(model_i), ties = "efron", data = data_study)
myCPH_i

out_fname <- sprintf("%s/%s_cox_time_spl%dba.csv", dir_out2, basename(file_name), n_knots)
out_fname_zph <- sprintf("%s/%s_cox_time_zph_spl%dba.csv", dir_out2, basename(file_name), n_knots)

outputcox(myCPH_i, study_name, model_abbr, out_fname, out_fname_zph)

survfit_S_data <- data.frame()
if (exclude_sex == 1) {
  for (i_rLTL0 in value_rLTL0:value_rLTL1) {
    for (i_age in 1:num_ages_S) {
      if (n_knots == 3) {
        newdata_cph_S_i <- data.frame(
          BA1 = my_spline_basis_S[i_age, 1], BA2 = my_spline_basis_S[i_age, 2],
          BA3 = my_spline_basis_S[i_age, 3], rLTL0 = i_rLTL0
        )
      }
      if (n_knots == 5) {
        newdata_cph_S_i <- data.frame(
          BA1 = my_spline_basis_S[i_age, 1], BA2 = my_spline_basis_S[i_age, 2],
          BA3 = my_spline_basis_S[i_age, 3], BA4 = my_spline_basis_S[i_age, 4],
          BA5 = my_spline_basis_S[i_age, 5], rLTL0 = i_rLTL0
        )
      }

      survfit_S_i <- survfit(myCPH_i, newdata = newdata_cph_S_i, censor = FALSE)
      survfit_S_data_i <- data.frame(
        age0 = ages_S[i_age], rLTL0 = i_rLTL0, sexM = 0, time = survfit_S_i$time,
        n_risk = survfit_S_i$n.risk, n_event = survfit_S_i$n.event,
        n_censor = survfit_S_i$n.censor, surv = survfit_S_i$surv,
        cumhaz = survfit_S_i$cumhaz, std_err = survfit_S_i$std.err,
        lower = survfit_S_i$lower, upper = survfit_S_i$upper
      )

      survfit_S_data <- rbind(survfit_S_data, survfit_S_data_i)
    }
  }
} else {
  for (i_sex in 0:1) {
    for (i_rLTL0 in value_rLTL0:value_rLTL1) {
      for (i_age in 1:num_ages_S) {
        if (n_knots == 3) {
          newdata_cph_S_i <- data.frame(
            BA1 = my_spline_basis_S[i_age, 1], BA2 = my_spline_basis_S[i_age, 2],
            BA3 = my_spline_basis_S[i_age, 3], rLTL0 = i_rLTL0, sexM = i_sex
          )
        }
        if (n_knots == 5) {
          newdata_cph_S_i <- data_frame(
            BA1 = my_spline_basis_S[i_age, 1], BA2 = my_spline_basis_S[i_age, 2],
            BA3 = my_spline_basis_S[i_age, 3], BA4 = my_spline_basis_S[i_age, 4],
            BA5 = my_spline_basis_S[i_age, 5], rLTL0 = i_rLTL0, sexM = i_sex
          )
        }

        survfit_S_i <- survfit(myCPH_i, newdata = newdata_cph_S_i, censor = FALSE)
        survfit_S_data_i <- data.frame(
          age0 = ages_S[i_age], rLTL0 = i_rLTL0, sexM = i_sex, time = survfit_S_i$time,
          n_risk = survfit_S_i$n.risk, n_event = survfit_S_i$n.event,
          n_censor = survfit_S_i$n.censor, surv = survfit_S_i$surv,
          cumhaz = survfit_S_i$cumhaz, std_err = survfit_S_i$std.err,
          lower = survfit_S_i$lower, upper = survfit_S_i$upper
        )

        survfit_S_data <- rbind(survfit_S_data, survfit_S_data_i)
      }
    }
  }
}

write.csv(survfit_S_data, sprintf("%s/%s_S_spl%dba_all.csv", dir_out2, basename(file_name), n_knots), row.names = FALSE, na = "")

if (exclude_sex == 1) {
  if (n_knots == 3) {
    newdata_cph <- data.frame(BA1 = 0, BA2 = 0, BA3 = 0, rLTL0 = 0)
  }
  if (n_knots == 5) {
    newdata_cph <- data.frame(BA1 = 0, BA2 = 0, BA3 = 0, BA4 = 0, BA5 = 0, rLTL0 = 0)
  }
} else {
  if (n_knots == 3) {
    newdata_cph <- data.frame(BA1 = 0, BA2 = 0, BA3 = 0, rLTL0 = 0, sexM = 0)
  }
  if (n_knots == 5) {
    newdata_cph <- data.frame(BA1 = 0, BA2 = 0, BA3 = 0, BA4 = 0, BA5 = 0, rLTL0 = 0, sexM = 0)
  }
}

survfit_S0 <- survfit(myCPH_i, newdata = newdata_cph, censor = FALSE)
survfit_S0_data <- data.frame(
  time = survfit_S0$time, n_risk = survfit_S0$n.risk, n_event = survfit_S0$n.event,
  n_censor = survfit_S0$n.censor, surv = survfit_S0$surv, cumhaz = survfit_S0$cumhaz,
  std_err = survfit_S0$std.err, lower = survfit_S0$lower, upper = survfit_S0$upper
)

write.csv(survfit_S0_data, sprintf("%s/%s_S0_spl%dba_all.csv", dir_out2, basename(file_name), n_knots), row.names = FALSE, na = "")

var_betas <- as.data.frame(myCPH_i$var)
names(var_betas) <- names(myCPH_i$coefficients)
row.names(var_betas) <- names(var_betas)

out_fname_var <- sprintf("%s/%s_var_spl%dba_all.csv", dir_out2, basename(file_name), n_knots)

write.csv(var_betas, file = out_fname_var, row.names = TRUE, na = "")

my_spline_basis_out_all <- splines::ns(ages_spline_out_all, knots = quant_aget, Boundary.knots = c(min_aget, max_aget))
row.names(my_spline_basis_out_all) <- sprintf("%.2f", ages_spline_out_all)

write.csv(my_spline_basis_out_all, sprintf("%s/%s_age0_spl%dba_all.csv", dir_out2, basename(file_name), n_knots), row.names = FALSE, na = "")

