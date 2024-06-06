# Function to align the alleles so the effects are for the same allele.
align_effect = function(dat, ind, flip = F) {
  if (length(ind) == 0) { return(dat) }
  res = dat
  if (!flip) {
    tmp = res$effect_allele[ind]
    res$effect_allele[ind] = res$other_allele[ind]
    res$other_allele[ind] = tmp
  }
  res$MAF[ind] = 1 - res$MAF[ind]
  res$effect[ind] = -res$effect[ind]
  return(res)
}

# Function to run Mendelian randomisation (X = exposure, Y = outcome)
# dat1 = exposure, dat2 = outcome
run_mr = function(dat1, dat2) {
  input = mr_input(bx = dat1$effect, bxse = dat1$SE, by = dat2$effect, byse = dat2$SE)
  # Run various MR methods:
  ivw = mr_ivw(input)
  med = mr_median(input, weighting = 'weighted', iterations = 10000, distribution = 'normal')
  egger = mr_egger(input)
  # Pull out relevant stats:
  ivw_res = data.frame(method = 'IVW',
                       estimate = ivw@Estimate,
                       se = ivw@StdError,
                       lower = ivw@CILower,
                       upper = ivw@CIUpper,
                       p = ivw@Pvalue)
  med_res = data.frame(method = 'Weighted Median',
                       estimate = med@Estimate,
                       se = med@StdError,
                       lower = med@CILower,
                       upper = med@CIUpper,
                       p = med@Pvalue)
  egger_res = data.frame(method = c('MR-Egger', 'MR-Egger Intercept'),
                         estimate = c(egger@Estimate, egger@Intercept),
                         se = c(egger@StdError.Est, egger@StdError.Int),
                         lower = c(egger@CILower.Est, egger@CILower.Int),
                         upper = c(egger@CIUpper.Est, egger@CIUpper.Int),
                         p = c(egger@Pvalue.Est, egger@Pvalue.Int))
  res = rbind(ivw_res, med_res, egger_res)
  return(res)
}

