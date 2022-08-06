#' Evaluate MENDEL Output for Simulated Data
#' 
#' Retrieves the estimated distribution means and variances separately for 
#' carriers and non-carriers and compares them to the true values.
#' 
#' @param distribution cancer distribution function, one of 
#' `c("exponential","binomial","normal")`.
#' @param age numeric age value for which to retrieve the cumulative risk. Should
#' be `NA` if distribution is `"binomial"`.
#' @param c.true.mean numeric mean value of the carrier cancer distribution.
#' @param c.true.var numeric variance value of the carrier cancer distribution.
#' Optional for distributions where the variance is just a transformation of the 
#' mean.
#' @param nc.true.mean numeric mean value of the non-carrier cancer distribution.
#' @param nc.true.var numeric variance value of the non-carrier cancer distribution.
#' Optional for distributions where the variance is just a transformation of the 
#' mean.
#' @param intercept numeric value of the GRAND PARAMETER from the MENDEL output.
#' @param age.coef numeric value of the AGE PARAMETER from the MENDEL output.
#' @param carrier.coef numeric value of the 1/2&2/2 PARAMETER from the MENDEL output.
#' Note, specifying the non-carrier coefficient is unnecessary because it is 
#' always the negative of `carrier.coef`.
#' @returns a list with elements:
#'  - Summary: a two row data frame where one row is for the carrier distribution 
#'  and the other row is for the non-carrier distribution with columns:
#'   - DistType: same as `distribution` argument
#'   - TrueMean: same as `c.true.mean` and `nc.true.mean` arguments.
#'   - EstimatedMean: estimated distribution means.
#'   - DiffMeans: the true mean minus the estimated mean.
#'   - PercDiffMeans: percent difference of the estimated mean from the true mean.
#'   - Age: same as `age` argument.
#'   - TrueRisk: cumulative cancer risk based on the true distribution parameters.
#'   - EstimatedRisk: estimated cumulative cancer risk based on the MENDEL output.
#'   - DiffRisk: the true risk minus the estimated risk.
#'   - PercDiffRis: the percent difference of the estimated risk from the true risk.
#'   - TrueVar: the distribution's true variance.
#'   - EstimatedVar: the estimated distribution variance based on the MENDEL output.
#'   - DiffVars: the true variance minus the estimated variance.
#'   - PercDiffVars: the percent difference of the estimated variance from the 
#'   true variance.
#'  - Data: a data frame with columns:
#'   - Age: numeric, ranges from 1 to 95.
#'   - Carrier_Status: string, "Carrier" or "Non-carrier".
#'   - Type: string, "True" or "Estimated"
#'   - Mean: numeric value of the cancer distribution mean.
#'   - Variance: numeric value of the cancer distribution variance.
#'   - Penetrance: age-conditional cancer penetrances (the PDF).
#'   - Risk: cumulative cancer risk (the CDF).
#'  - Plot: a two facet ggplot line plot of the cancer distributions where age 
#'  is the x-axis and cumulative probability is the y-axis. Each facet has two
#'  lines:
#'   - facet 1: true and estimated carrier CDF
#'   - facet 2: true and estimated non-carrier CDF
test.MENDEL.output <- function(distribution, age = NA, 
                               c.true.mean, c.true.var = NA, 
                               nc.true.mean, nc.true.var = NA,
                               intercept, age.coef, carrier.coef){
  
  max.age <- 95
  age.range <- 1:max.age
  
  # repeat analysis for carriers and non-carriers
  for(c.status in c("c","nc")){
    
    # get parameters
    if(c.status == "c"){
      true.mean <- c.true.mean
      true.var <- c.true.var
      gene.coef <- carrier.coef
    } else if(c.status == "nc"){
      true.mean <- nc.true.mean
      true.var <- nc.true.var
      gene.coef <- -1 * carrier.coef
    }
    
    # true cancer risk
    if(distribution == "exponential"){
      # CDF
      true.risk <- 1 - exp(-age / true.mean)
    } else if(distribution == "normal"){
      # CDF
      true.risk <- pnorm(q = age, mean = true.mean, sd = sqrt(true.var))
    } else if(distribution == "binomial"){
      # identity
      true.risk <- true.mean
    }
    
    # estimated linear function value
    if(distribution %in% c("binomial","negative binomial")){
      linear.expected <- intercept + gene.coef
    } else {
      linear.expected <- intercept + gene.coef + age.coef * age
    }
    
    # estimated cancer risk using inverse link function
    if(distribution %in% c("exponential", "poisson", "gamma")){ 
      # inverse log
      estimated.risk <- exp(linear.expected)
    } else if(distribution %in% c("inverse gaussian", "logistic", "lognormal", "normal")){ 
      # identity
      estimated.risk <- linear.expected
    } else if(distribution %in% c("binomial", "negative binomial")){ 
      # inverse logit
      estimated.risk <- exp(linear.expected) / (1 + exp(linear.expected))
    }
    
    # estimated mean
    if(distribution == "exponential"){
      # solve CDF for mu using estimated risk
      estimated.mean <- -age / log(1 - estimated.risk)
    } else if(distribution == "normal"){
      # use points along the CDF curve to get the PDF curve, then estimate
      cdf.points <- intercept + gene.coef + age.coef * 0:10000
      pdf.points <- c(cdf.points[1], diff(cdf.points))
      estimated.mean <- mean(pdf.points)
    } else if(distribution == "binomial"){
      # identity
      estimated.mean <- estimated.risk
    }
    
    # if true variances not provided for distributions where the variances is
    # a transformation of the mean, calculate them
    if(is.na(true.var)){
      if(distribution == "exponential"){
        true.var <- true.mean^2
      } else if(distribution == "binomial"){
        true.var <- true.mean * (1-true.mean)
      }
    }
    
    # estimated variance
    if(distribution == "exponential"){
      estimated.var <- estimated.mean^2
    } else if(distribution == "normal"){
      estimated.var <- var(pdf.points)
    } else if(distribution == "binomial"){
      estimated.var <- estimated.mean * (1-estimated.mean)
    }
    
    # compare means and risk
    diff.means <- true.mean - estimated.mean
    perc.diff.means <- abs(diff.means) / true.mean * 100
    diff.risk <- true.risk - estimated.risk
    perc.diff.risk <- abs(diff.risk) / true.risk * 100
    diff.vars <- true.var - estimated.var
    perc.diff.vars <- abs(diff.vars) / true.var * 100
    
    # store results
    carrier.status <- ifelse(c.status == "c", "Carrier", "Non-carrier")
    tmp.df <- data.frame(DistType = distribution, 
                         CarrierStatus = carrier.status,
                         TrueMean = true.mean, 
                         EstimatedMean = round(estimated.mean, 3), 
                         DiffMeans = round(diff.means, 3),
                         PercDiffMeans = round(perc.diff.means, 1),
                         Age = age,
                         TrueRisk = round(true.risk, 3),
                         EstimatedRisk = round(estimated.risk, 3),
                         DiffRisk = round(diff.risk, 3),
                         PercDiffRisk = round(perc.diff.risk, 1),
                         TrueVar = round(true.var, 3),
                         EstimatedVar = round(estimated.var, 3),
                         DiffVars = round(diff.vars, 3),
                         PercDiffVars = round(perc.diff.vars, 1))
    tmp.plot.df <- data.frame(Age = rep(age.range, times = 2),
                              Carrier_Status = rep(carrier.status, times = 2*max.age),
                              Type = rep(c("True","Estimated"), each = max.age),
                              Mean = rep(c(true.mean, estimated.mean), each = max.age),
                              Variance = rep(c(true.var, estimated.var), each = max.age),
                              Penetrance = rep(0, times = 2*max.age),
                              Risk = rep(0, times = 2*max.age))
    tmp.plot.df <- 
      tmp.plot.df %>%
      mutate(Penetrance = dexp(x = Age, rate = 1/Mean)) %>%
      mutate(Risk = pexp(q = Age, rate = 1/Mean))
    if(c.status == "c"){
      df <- tmp.df
      plot.df <- tmp.plot.df
    } else if(c.status == "nc"){
      df <- rbind(df, tmp.df)
      plot.df <- rbind(plot.df, tmp.plot.df)
    }
  }
  
  ## plot the CDFs
  plot <- 
    ggplot(plot.df, aes(x = Age, y = Risk, color = Type)) +
    geom_line() +
    facet_grid(. ~ Carrier_Status) +
    scale_x_continuous(limits = c(0,100)) +
    scale_y_continuous(limits = c(0,1)) +
    labs(title = "Cumulative Cancer Risk by Age and Carrier Status",
         x = "Age",
         y = "Cumulative Risk") +
    theme_bw()
  
  list(Summary = df, Data = plot.df, Plot = plot)
}

