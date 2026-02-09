################################################################################
# APPENDIX S2.1
#
# Episodic Birth-Death model: CoMET implemented in R 
# CoMET with rjMCMC (Reversible Jump MCMC algorithm) to estimate simultaneously 
# changes in diversification rates and mass extinction events
#
################################################################################

# Load packages ape and coda and desolve as dependencies
library(TESS)

treefile <- "species"

# Use this line for the first two analyses
analysis_name <- sprintf("CoMET_%s",treefile) 

# Use this line for the third analyses: Incomplete Taxon Sampling Strategy
analysis_name <- sprintf("CoMET_%s_ME-Div",treefile) 

# Conditioning on survival of all lineages in empirical phylogeny
CDT <- "survival"

# This assigns the prior for changes in speciation and extinction (or net diversification and turnover) as Poisson process priors with 0.5 probability assigned to 0 MEE or 0 rate shifts events
EXPECTED_NUM_EVENTS <- 2

MCMC_ITERATIONS <- 1000000

# Read pruned tree containing only 1 tip per species in phylogeny without outgroups (the other 3 subgenera) from 338 phylogeny, leaving only 322 tips and 321 internal nodes
tree <- read.tree(file="Esula_tree_no_outgroups.tre")
tree
plot(tree)
axisPhylo()

# Assign sampling value at present
# Since there is incomplete taxon sampling, estimate this as fraction of total. THERE ARE 322 SPECIES IN PHYLOGENY (322/489 = 0.66)
total <- 489
rho <- (tree$Nnode+1)/total # rho <- 0.66
priorForms <- c("lognormal","normal","gamma")

# RUN ANALYSIS #
# Run two different analysis, one with MEE sampling events (integrating out rate shifts) and one without MEEs, estimating posterior of rate shifts
# First, prepare the names for the files, using this "if" condition: if ALLOW_MASS_EXTINCTION is FALSE (!TRUE), we use "analysis_name" como CoMET_species; if ALLOW_MASS_EXTINCTION is TRUE (== TRUE), we use  "analysis_name" as "CoMET_species_MEE"
ALLOW_MASS_EXTINCTION <- !TRUE
if ( ALLOW_MASS_EXTINCTION == TRUE ) {
  analysis_name <- sprintf("CoMET_%s_ME",treefile)
} else {
  analysis_name <- sprintf("CoMET_%s",treefile)
}

# To run an analysis with no MEEs, use the code [1] below. This will produce an output called "CoMET_species" containing the files with parameter values. It will also produce at the end a PDF figure named "CoMET_species.pdf" containing vignettes for 4 figures. This function will create a directory with the name of "analysis_name" given above in the current working directory getwd()
# IMPORTANT! Before running the code [2], change the name of the analysis' folder [1], to allow overwritting files. Eg: CoMET_species_noMEE for code [1]
#[1]
tess.analysis(tree=tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS, numExpectedMassExtinctions=EXPECTED_NUM_EVENTS, initialSpeciationRate=2.0, initialExtinctionRate=1.0, empiricalHyperPriorInflation = 10.0, empiricalHyperPriorForm = priorForms, samplingProbability=rho, estimateMassExtinctionTimes = ALLOW_MASS_EXTINCTION, estimateNumberMassExtinctions = ALLOW_MASS_EXTINCTION, MAX_ITERATIONS = MCMC_ITERATIONS, THINNING = 100,  MAX_TIME = Inf, MIN_ESS = 1000, CONDITION=CDT, dir = analysis_name)

# PROCESS ANALYSIS OUTPUT #
out <- tess.process.output(analysis_name, tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS)

# CHECK MCMC DIAGNOSTICS WITH CODA #
#If only one chain and with EES diagnostic
#Remember to give only one parameter ("speciation rates") or divide into vignettes with "par" function. Otherwise, plots will be overlaid. Similarly, choose between "EES" and "geweke" to avoid overlaying of values. There are other arguments, which we do not consider here: (col=NULL, xaxt="n", yaxt="s")

# Run this line only when I had the three analyses with folders rename different.  Run this line with all lines above and end.
# Export plots as .pdf with A4 size
par(mfrow= c(4, 2))
TessPlot <- tess.plot.singlechain.diagnostics(out, parameters=c("speciation rates", "extinction rates", "speciation shift times", "extinction shift times", "net-diversification rates", "relative-extinction rates", "mass extinction times"), diagnostics=c("ESS"), ess.crit=c(100,200), correction="bonferroni", xlab="million years ago", pch=19)
dev.off()

# PLOT RESULTS IN VIGNETTES #
# Results will be plotted as vignettes within a PDF: nrows = 2, ncolumns= FIG_TYPES/2
# Figures will be plotted as vignettes: number of figures depends on Fig_Types. Rate shift times and mass extinction times are plotted together with Bayes Factors (2lnBF) significance levels at 2,6,10 on the right y axis (this can be assigned a name with argument "yaxt")
NUM_FIGS <- 6
FIG_TYPES <- c("speciation rates","extinction rates","net-diversification rates","relative-extinction rates", "speciation shift times","extinction shift times")

pdf(sprintf("%s.pdf",analysis_name))
layout.mat <- matrix(1:NUM_FIGS, nrow=2, ncol=NUM_FIGS / 2)
layout(layout.mat)
tess.plot.output(out, fig.types=FIG_TYPES, las=1)
dev.off()

# OPTIONS #
# To run an analysis with MEEs, use the code [2] below. This will produce an output folder called "CoMET_ME_species" containing the files with parameter posterior MCMC values. It will also produce at the end a PDF figure named "CoMET_ME_species.pdf" containing vignettes for six figures.
# IMPORTANT! Before running the code "CHANGE SAMPLING STRATEGY FROM UNIFORM TO DIVERSIFIED, change the name of the analysis' folder [2], to allow overwritting files. Eg: CoMET_species_MEE for code [2]
#[2]
ALLOW_MASS_EXTINCTION <- TRUE

tess.analysis(tree=tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS, numExpectedMassExtinctions=EXPECTED_NUM_EVENTS, initialSpeciationRate=2.0, initialExtinctionRate=1.0, empiricalHyperPriorInflation = 10.0, empiricalHyperPriorForm = priorForms, samplingProbability=rho, estimateMassExtinctionTimes = ALLOW_MASS_EXTINCTION, estimateNumberMassExtinctions = ALLOW_MASS_EXTINCTION, MAX_ITERATIONS = MCMC_ITERATIONS, THINNING = 100,  MAX_TIME = Inf, MIN_ESS = 1000, CONDITION=CDT, dir = analysis_name)

# PROCESS ANALYSIS OUTPUT #
outME <- tess.process.output(analysis_name, tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS)
NUM_FIGS <- 6
FIG_TYPES <- c("speciation rates", "extinction rates", "net-diversification rates", "relative-extinction rates", "mass extinction Bayes factors", "mass extinction times")
pdf(sprintf("%s.pdf",analysis_name))
layout.mat <- matrix(1:NUM_FIGS, nrow=2, ncol=NUM_FIGS / 2)
layout(layout.mat)
tess.plot.output(outME, fig.types=FIG_TYPES, las=1)
dev.off()

# CHANGE SAMPLING STRATEGY FROM UNIFORM TO DIVERSIFIED #
# INCOMPLETE TAXON SAMPLING Strategy. Default in tess.analysis function is "samplingStrategy=uniform" (species are sampled uniformly at random), which means that each species has the same probability (rho) of being sampled in the phylogeny. 
# To apply alternative: "samplingStrategy=diversified" (species are sampled to maximize the diversity sampled in the phylogeny: i.e., only the oldest 25% of divergence events are included in the reconstructed phylogeny with sampling probability rho, and all later divergence events are excluded), we can use the function "tess.analysis.diversified.R" created by modifying the argument "samplingStrategy=diversified" in the "tess.likelihood" function part of the code. First, source it:
# IMPORTANT! 1) Before running the next lines, NO run the line 10: analysis_name <- sprintf("CoMET_%s",treefile), don't run it and copy, paste in its place, and run the line: analysis_name <- sprintf("CoMET_%s_ME-Div",treefile), and don't run it again on the next part of the analysis (below).
# IMPORTANT! 2) Only run the first 60 lines, until the end of the loop of: ALLOW_MASS_EXTINCTION <- !TRUE ....   analysis_name <- sprintf("CoMET_%s",treefile)}})
source("tess.analysis.diversified.R") # See the end of APPENDIX S2 (tess.analysis.diversified.R) to find the function

# Next, run the analysis. Specify the analysis_name with the command below.
analysis_name <- sprintf("CoMET_%s_ME-Div",treefile)

tess.analysis.diversified(tree=tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS, numExpectedMassExtinctions=EXPECTED_NUM_EVENTS, initialSpeciationRate=2.0, initialExtinctionRate=1.0, empiricalHyperPriorInflation = 10.0, empiricalHyperPriorForm = priorForms, samplingProbability=rho, estimateMassExtinctionTimes = ALLOW_MASS_EXTINCTION, estimateNumberMassExtinctions = ALLOW_MASS_EXTINCTION, MAX_ITERATIONS = MCMC_ITERATIONS, THINNING = 100,  MAX_TIME = Inf, MIN_ESS = 1000, CONDITION=CDT, dir = analysis_name)   

# PROCESS ANALYSIS OUTPUT #
outDiv <- tess.process.output(analysis_name, tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS)
NUM_FIGS <- 6
FIG_TYPES <- c("net-diversification rates", "relative-extinction rates", "speciation shift times", "extinction shift times", "mass extinction Bayes factors", "mass extinction times")

pdf(sprintf("%s.pdf",analysis_name))
layout.mat <- matrix(1:NUM_FIGS, nrow=2, ncol=NUM_FIGS / 2)
layout(layout.mat)
tess.plot.output(outDiv, fig.types=FIG_TYPES, las=1)
dev.off()



###################################
### tess.analysis.diversified.R ###
###################################

# Script with the necessary function to run the CoMET analysis in R (APPENDIX S2.1):
tess.analysis.diversified <- function (tree, initialSpeciationRate, initialExtinctionRate, 
                                       empiricalHyperPriors = TRUE, empiricalHyperPriorInflation = 10, 
                                       empiricalHyperPriorForm = c("lognormal", "normal", "gamma"), 
                                       speciationRatePriorMean = 0, speciationRatePriorStDev = 1, 
                                       extinctionRatePriorMean = 0, extinctionRatePriorStDev = 1, 
                                       initialSpeciationRateChangeTime = c(), initialExtinctionRateChangeTime = c(), 
                                       estimateNumberRateChanges = TRUE, numExpectedRateChanges = 2, 
                                       samplingProbability = 1, missingSpecies = c(), timesMissingSpecies = c(), 
                                       tInitialMassExtinction = c(), pInitialMassExtinction = c(), 
                                       pMassExtinctionPriorShape1 = 5, pMassExtinctionPriorShape2 = 95, 
                                       estimateMassExtinctionTimes = TRUE, numExpectedMassExtinctions = 2, 
                                       estimateNumberMassExtinctions = TRUE, MRCA = TRUE, CONDITION = "survival", 
                                       BURNIN = 10000, MAX_ITERATIONS = 2e+05, THINNING = 100, OPTIMIZATION_FREQUENCY = 500, 
                                       CONVERGENCE_FREQUENCY = 1000, MAX_TIME = Inf, MIN_ESS = 500, 
                                       ADAPTIVE = TRUE, dir = "", priorOnly = FALSE, verbose = TRUE) 
{
  orgDir <- getwd()
  if (dir != "") {
    if (file.exists(dir)) {
      setwd(file.path(dir))
    }
    else {
      dir.create(file.path(dir))
      setwd(file.path(dir))
    }
  }
  if (BURNIN > 0 & BURNIN < 1) {
    burn <- MAX_ITERATIONS * BURNIN
  }
  else {
    burn <- BURNIN
  }
  MAX_TRIES <- 1000
  times <- branching.times(tree)
  if (!is.ultrametric(tree)) {
    stop("The likelihood function is only defined for ultrametric trees!")
  }
  write.nexus(tree, file = "input_tree.tre")
  if (empiricalHyperPriors == TRUE) {
    if (verbose) {
      cat("Estimating empirical hyper-parameters.\n\n")
    }
    empirical.prior.likelihood <- function(params) {
      b <- params[1]/(1 - params[2])
      d <- params[2] * b
      lnl <- tess.likelihood(times, b, d, missingSpecies = missingSpecies, 
                             timesMissingSpecies = timesMissingSpecies, samplingStrategy = "diversified", 
                             samplingProbability = samplingProbability, MRCA = MRCA, 
                             CONDITION = CONDITION, log = TRUE)
      return(lnl)
    }
    prior.diversification <- function(x) {
      dexp(x, rate = max(times)/log(length(times)), log = TRUE)
    }
    prior.turnover <- function(x) {
      dbeta(x, shape1 = 1.5, shape2 = 3, log = TRUE)
    }
    priors <- c(diversification = prior.diversification, 
                turnover = prior.turnover)
    tries <- 1
    repeat {
      starting.values <- runif(2, 0, 1)
      starting.lnl <- empirical.prior.likelihood(starting.values)
      if (is.finite(starting.lnl) || tries >= MAX_TRIES) {
        break
      }
      tries <- tries + 1
    }
    if (tries >= MAX_TRIES) {
      stop("Could not find valid starting values after ", 
           tries, " attemps.\n")
    }
    samples <- tess.mcmc(empirical.prior.likelihood, priors, 
                         starting.values, logTransforms = c(TRUE, TRUE), delta = c(0.1, 
                                                                                   0.1), iterations = 20000, burnin = 2000, verbose = verbose)
    samples.lambda <- samples[, 1]/(1 - samples[, 2])
    m.lambda <- mean(samples.lambda)
    v.lambda <- empiricalHyperPriorInflation * var(samples.lambda)
    samples.mu <- samples.lambda * samples[, 2]
    m.mu <- mean(samples.mu)
    v.mu <- empiricalHyperPriorInflation * var(samples.mu)
    lambda.hyperprior.fits <- c(lognormal = Inf, normal = Inf, 
                                gamma = Inf)
    mu.hyperprior.fits <- c(lognormal = Inf, normal = Inf, 
                            gamma = Inf)
    if ("lognormal" %in% empiricalHyperPriorForm) {
      lambda.hyperprior.fits["lognormal"] <- suppressWarnings(ks.test(samples.lambda, 
                                                                      "plnorm", meanlog = log((m.lambda^2)/sqrt(v.lambda + 
                                                                                                                  m.lambda^2)), sdlog = sqrt(log(1 + v.lambda/(m.lambda^2)))))$statistic
      mu.hyperprior.fits["lognormal"] <- suppressWarnings(ks.test(samples.mu, 
                                                                  "plnorm", meanlog = log((m.mu^2)/sqrt(v.mu + 
                                                                                                          m.mu^2)), sdlog = sqrt(log(1 + v.mu/(m.mu^2)))))$statistic
    }
    if ("normal" %in% empiricalHyperPriorForm) {
      lambda.hyperprior.fits["normal"] <- suppressWarnings(ks.test(samples.lambda, 
                                                                   "pnorm", mean = m.lambda, sd = sqrt(v.lambda)))$statistic
      mu.hyperprior.fits["normal"] <- suppressWarnings(ks.test(samples.mu, 
                                                               "pnorm", mean = m.mu, sd = sqrt(v.mu)))$statistic
    }
    if ("gamma" %in% empiricalHyperPriorForm) {
      lambda.hyperprior.fits["gamma"] <- suppressWarnings(ks.test(samples.lambda, 
                                                                  "pgamma", shape = (m.lambda^2)/v.lambda, rate = m.lambda/v.lambda))$statistic
      mu.hyperprior.fits["gamma"] <- suppressWarnings(ks.test(samples.mu, 
                                                              "pgamma", shape = (m.mu^2)/v.mu, rate = m.mu/v.mu))$statistic
    }
    lambda.hyper.prior.form <- names(which.min(lambda.hyperprior.fits[empiricalHyperPriorForm]))
    mu.hyper.prior.form <- names(which.min(mu.hyperprior.fits[empiricalHyperPriorForm]))
    if (lambda.hyper.prior.form == "lognormal") {
      speciationRatePriorMean <- log((m.lambda^2)/sqrt(v.lambda + 
                                                         m.lambda^2))
      speciationRatePriorStDev <- sqrt(log(1 + v.lambda/(m.lambda^2)))
    }
    else if (lambda.hyper.prior.form == "normal") {
      speciationRatePriorMean <- m.lambda
      speciationRatePriorStDev <- sqrt(v.lambda)
    }
    else if (lambda.hyper.prior.form == "gamma") {
      speciationRatePriorMean <- (m.lambda^2)/v.lambda
      speciationRatePriorStDev <- m.lambda/v.lambda
    }
    if (mu.hyper.prior.form == "lognormal") {
      extinctionRatePriorMean <- log((m.mu^2)/sqrt(v.mu + 
                                                     m.mu^2))
      extinctionRatePriorStDev <- sqrt(log(1 + v.mu/(m.mu^2)))
    }
    else if (mu.hyper.prior.form == "normal") {
      extinctionRatePriorMean <- m.mu
      extinctionRatePriorStDev <- sqrt(v.mu)
    }
    else if (mu.hyper.prior.form == "gamma") {
      extinctionRatePriorMean <- (m.mu^2)/v.mu
      extinctionRatePriorStDev <- m.mu/v.mu
    }
    initialSpeciationRate <- m.lambda
    initialExtinctionRate <- m.mu
    empirical.hyperprior.samples <- data.frame(1:nrow(samples), 
                                               samples.lambda, samples.mu)
    write.table(empirical.hyperprior.samples[-1, ], file = "samplesHyperprior.txt", 
                row.names = FALSE, col.names = c("Iteration", "speciation", 
                                                 "extinction"), sep = "\t")
    summary.file <- "empiricalHyperpriorSummary.txt"
    cat("Summary of the empirical hyperprior analysis", file = summary.file)
    cat("\n============================================\n\n", 
        file = summary.file, append = TRUE)
    cat("Specation rate hyperprior\n", file = summary.file, 
        append = TRUE)
    cat("-------------------------\n", file = summary.file, 
        append = TRUE)
    cat("Form:\t", lambda.hyper.prior.form, "\n", file = summary.file, 
        append = TRUE)
    if (lambda.hyper.prior.form == "lognormal") {
      cat("meanlog:\t", speciationRatePriorMean, "\n", 
          file = summary.file, append = TRUE)
      cat("sdlog:\t\t", speciationRatePriorStDev, "\n", 
          file = summary.file, append = TRUE)
    }
    else if (lambda.hyper.prior.form == "normal") {
      cat("mean:\t", speciationRatePriorMean, "\n", file = summary.file, 
          append = TRUE)
      cat("sd:\t\t", speciationRatePriorStDev, "\n", file = summary.file, 
          append = TRUE)
    }
    else if (lambda.hyper.prior.form == "gamma") {
      cat("shape:\t", speciationRatePriorMean, "\n", file = summary.file, 
          append = TRUE)
      cat("rate:\t\t", speciationRatePriorStDev, "\n", 
          file = summary.file, append = TRUE)
    }
    cat("\nExtinction rate hyperprior\n", file = summary.file, 
        append = TRUE)
    cat("-------------------------\n", file = summary.file, 
        append = TRUE)
    cat("Form:\t", mu.hyper.prior.form, "\n", file = summary.file, 
        append = TRUE)
    if (mu.hyper.prior.form == "lognormal") {
      cat("meanlog:\t", extinctionRatePriorMean, "\n", 
          file = summary.file, append = TRUE)
      cat("sdlog:\t\t", extinctionRatePriorStDev, "\n", 
          file = summary.file, append = TRUE)
    }
    else if (mu.hyper.prior.form == "normal") {
      cat("mean:\t", extinctionRatePriorMean, "\n", file = summary.file, 
          append = TRUE)
      cat("sd:\t\t", extinctionRatePriorStDev, "\n", file = summary.file, 
          append = TRUE)
    }
    else if (mu.hyper.prior.form == "gamma") {
      cat("shape:\t", extinctionRatePriorMean, "\n", file = summary.file, 
          append = TRUE)
      cat("rate:\t\t", extinctionRatePriorStDev, "\n", 
          file = summary.file, append = TRUE)
    }
  }
  else {
    lambda.hyper.prior.form <- mu.hyper.prior.form <- empiricalHyperPriorForm[1]
  }
  likelihood <- function(lambda, lambdaTimes, mu, muTimes, 
                         tMassExtinction, pSurvival) {
    if (priorOnly == FALSE) {
      if (any(lambda < 0) | any(mu < 0)) {
        lnl <- -Inf
      }
      else {
        lnl <- tess.likelihood.rateshift(times, lambda, 
                                         mu, rateChangeTimesLambda = lambdaTimes, rateChangeTimesMu = muTimes, 
                                         massExtinctionTimes = tMassExtinction, massExtinctionSurvivalProbabilities = pSurvival, 
                                         missingSpecies = missingSpecies, timesMissingSpecies = timesMissingSpecies, 
                                         samplingStrategy = "uniform", samplingProbability = samplingProbability, 
                                         MRCA = MRCA, CONDITION = CONDITION, log = TRUE)
      }
    }
    else {
      lnl <- 0
    }
    return(lnl)
  }
  prior <- function(timesLambda, valuesLambda, timesMu, valuesMu, 
                    timesMassExtinction, valuesMassExtinction) {
    if (any(valuesLambda < 0) | any(valuesMu < 0)) {
      return(-Inf)
    }
    k <- length(timesLambda)
    lnp <- dpois(k, numExpectedRateChanges, log = TRUE)
    if (k > 0) {
      lnp <- lnp + sum(dunif(timesLambda, 0, AGE, log = TRUE))
    }
    if (lambda.hyper.prior.form == "lognormal") {
      lnp <- lnp + sum(dlnorm(valuesLambda, speciationRatePriorMean, 
                              speciationRatePriorStDev, log = TRUE))
    }
    else if (lambda.hyper.prior.form == "normal") {
      lnp <- lnp + sum(dnorm(valuesLambda, speciationRatePriorMean, 
                             speciationRatePriorStDev, log = TRUE))
    }
    else if (lambda.hyper.prior.form == "gamma") {
      lnp <- lnp + sum(dgamma(valuesLambda, speciationRatePriorMean, 
                              speciationRatePriorStDev, log = TRUE))
    }
    k <- length(timesMu)
    lnp <- lnp + dpois(k, numExpectedRateChanges, log = TRUE)
    if (k > 0) {
      lnp <- lnp + sum(dunif(timesMu, 0, AGE, log = TRUE))
    }
    if (mu.hyper.prior.form == "lognormal") {
      lnp <- lnp + sum(dlnorm(valuesMu, extinctionRatePriorMean, 
                              extinctionRatePriorStDev, log = TRUE))
    }
    else if (mu.hyper.prior.form == "normal") {
      lnp <- lnp + sum(dnorm(valuesMu, extinctionRatePriorMean, 
                             extinctionRatePriorStDev, log = TRUE))
    }
    else if (mu.hyper.prior.form == "gamma") {
      lnp <- lnp + sum(dgamma(valuesMu, extinctionRatePriorMean, 
                              extinctionRatePriorStDev, log = TRUE))
    }
    k <- length(timesMassExtinction)
    lnp <- lnp + dpois(k, numExpectedMassExtinctions, log = TRUE)
    if (k > 0) {
      lnp <- lnp + sum(dunif(timesMassExtinction, 0, AGE, 
                             log = TRUE))
      lnp <- lnp + sum(dbeta(valuesMassExtinction, pMassExtinctionPriorShape1, 
                             pMassExtinctionPriorShape2, log = TRUE))
    }
    return(lnp)
  }
  AGE <- max(as.numeric(branching.times(tree)))
  posterior <- c()
  kLambda <- c()
  kMu <- c()
  kMassExtinction <- c()
  speciationRateValues <- list()
  speciationRateChangeTimes <- list()
  extinctionRateValues <- list()
  extinctionRateChangeTimes <- list()
  survivalProbability <- list()
  massExtinctionTime <- list()
  cat("SpeciationRates\n", file = "SpeciationRates.txt")
  cat("SpeciationRateChanges\n", file = "SpeciationRateChanges.txt")
  cat("ExtinctionRates\n", file = "ExtinctionRates.txt")
  cat("ExtinctionRateChangeTimes\n", file = "ExtinctionRateChanges.txt")
  cat("SurvivalProbabilities\n", file = "SurvivalProbabilities.txt")
  cat("MassExtinctionTimes\n", file = "MassExtinctionTimes.txt")
  cat(paste("Iteration", "posterior", "NumSpeciation", "numExtinction", 
            "numMassExtinctions", sep = "\t"), "\n", file = "samples_numCategories.txt")
  lambda <- initialSpeciationRate
  mu <- initialExtinctionRate
  lambdaChangeTimes <- initialSpeciationRateChangeTime
  muChangeTimes <- initialExtinctionRateChangeTime
  pMassExtinction <- pInitialMassExtinction
  tMassExtinction <- tInitialMassExtinction
  initialPP <- prior(lambdaChangeTimes, lambda, muChangeTimes, 
                     mu, tMassExtinction, pMassExtinction) + likelihood(lambda, 
                                                                        lambdaChangeTimes, mu, muChangeTimes, tMassExtinction, 
                                                                        pMassExtinction)
  currentPP <- initialPP
  MAX_TRIES <- 1000
  for (i in 1:MAX_TRIES) {
    tmpPP <- prior(lambdaChangeTimes, lambda, muChangeTimes, 
                   mu, tMassExtinction, pMassExtinction) + likelihood(lambda, 
                                                                      lambdaChangeTimes, mu, muChangeTimes, tMassExtinction, 
                                                                      pMassExtinction)
    if (is.finite(tmpPP)) 
      break
    lambdaChangeTimes <- c()
    lambda <- rlnorm(1, speciationRatePriorMean, speciationRatePriorStDev)
    muChangeTimes <- c()
    mu <- rlnorm(1, extinctionRatePriorMean, extinctionRatePriorStDev)
    tMassExtinction <- c()
    pMassExtinction <- c()
  }
  deltaLambdaTime <- 0.1
  lambdaTimeAccepted <- 0
  lambdaTimeTried <- 0
  deltaLambdaValue <- 0.1
  lambdaValueAccepted <- 0
  lambdaValueTried <- 0
  deltaMuTime <- 0.1
  muTimeAccepted <- 0
  muTimeTried <- 0
  deltaMuValue <- 0.1
  muValueAccepted <- 0
  muValueTried <- 0
  deltaMassExtinctionTime <- 0.1
  massExtinctionTimeAccepted <- 0
  massExtinctionTimeTried <- 0
  deltaMassExtinctionValue <- 0.1
  massExtinctionValueAccepted <- 0
  massExtinctionValueTried <- 0
  finished <- FALSE
  SAMPLE <- FALSE
  min.ess <- 0
  startTime <- Sys.time()
  i <- 1
  if (verbose) {
    cat("\nPerforming CoMET analysis.\n\n")
    if (burn > 0) {
      cat("Burning-in the chain ...\n")
    }
    else {
      cat("Running the chain ... \n")
    }
    cat("0--------25--------50--------75--------100\n")
    bar <- txtProgressBar(style = 1, width = 42)
  }
  while (!finished) {
    if (length(lambdaChangeTimes) > 0) {
      idx <- floor(runif(1, 0, length(lambdaChangeTimes))) + 
        1
      t <- lambdaChangeTimes[idx]
      oldPP <- currentPP
      tPrime <- rnorm(1, t, deltaLambdaTime)
      lambdaChangeTimes[idx] <- tPrime
      if (tPrime > 0 && tPrime < AGE) {
        newPP <- prior(lambdaChangeTimes, lambda, muChangeTimes, 
                       mu, tMassExtinction, pMassExtinction) + likelihood(lambda, 
                                                                          lambdaChangeTimes, mu, muChangeTimes, tMassExtinction, 
                                                                          pMassExtinction)
      }
      else {
        newPP <- -Inf
      }
      if (is.finite(newPP - oldPP) == FALSE || log(runif(1, 
                                                         0, 1)) > (newPP - oldPP)) {
        lambdaChangeTimes[idx] <- t
      }
      else {
        currentPP <- newPP
        lambdaTimeAccepted <- lambdaTimeAccepted + 1
      }
      lambdaTimeTried <- lambdaTimeTried + 1
    }
    idx <- floor(runif(1, 0, length(lambda))) + 1
    v <- lambda[idx]
    oldPP <- currentPP
    u <- rnorm(1, 0, deltaLambdaValue)
    scalingFactor <- exp(u)
    vPrime <- v * scalingFactor
    lambda[idx] <- vPrime
    lnHastingsratio <- log(scalingFactor)
    newPP <- prior(lambdaChangeTimes, lambda, muChangeTimes, 
                   mu, tMassExtinction, pMassExtinction) + likelihood(lambda, 
                                                                      lambdaChangeTimes, mu, muChangeTimes, tMassExtinction, 
                                                                      pMassExtinction)
    if (is.finite(newPP - oldPP + lnHastingsratio) == FALSE || 
        log(runif(1, 0, 1)) > (newPP - oldPP + lnHastingsratio)) {
      lambda[idx] <- v
    }
    else {
      currentPP <- newPP
      lambdaValueAccepted <- lambdaValueAccepted + 1
    }
    lambdaValueTried <- lambdaValueTried + 1
    if (estimateNumberRateChanges == TRUE) {
      u <- runif(1, 0, 1)
      if (u > 0.5) {
        oldPP <- currentPP
        t <- runif(1, 0, AGE)
        if (lambda.hyper.prior.form == "lognormal") {
          v <- rlnorm(1, speciationRatePriorMean, speciationRatePriorStDev)
        }
        else if (lambda.hyper.prior.form == "normal") {
          v <- rnorm(1, speciationRatePriorMean, speciationRatePriorStDev)
        }
        else if (lambda.hyper.prior.form == "gamma") {
          v <- rgamma(1, speciationRatePriorMean, speciationRatePriorStDev)
        }
        lambdaChangeTimes[length(lambdaChangeTimes) + 
                            1] <- t
        lambda[length(lambda) + 1] <- v
        newPP <- prior(lambdaChangeTimes, lambda, muChangeTimes, 
                       mu, tMassExtinction, pMassExtinction) + likelihood(lambda, 
                                                                          lambdaChangeTimes, mu, muChangeTimes, tMassExtinction, 
                                                                          pMassExtinction)
        kPrime <- length(lambdaChangeTimes)
        if (lambda.hyper.prior.form == "lognormal") {
          pr <- 1/(dunif(t, 0, AGE) * dlnorm(v, speciationRatePriorMean, 
                                             speciationRatePriorStDev))
        }
        else if (lambda.hyper.prior.form == "normal") {
          pr <- 1/(dunif(t, 0, AGE) * dnorm(v, speciationRatePriorMean, 
                                            speciationRatePriorStDev))
        }
        else if (lambda.hyper.prior.form == "gamma") {
          pr <- 1/(dunif(t, 0, AGE) * dgamma(v, speciationRatePriorMean, 
                                             speciationRatePriorStDev))
        }
        if (is.finite(newPP - oldPP + log(pr)) == FALSE || 
            log(runif(1, 0, 1)) > (newPP - oldPP + log(pr))) {
          lambdaChangeTimes <- lambdaChangeTimes[-length(lambdaChangeTimes)]
          lambda <- lambda[-length(lambda)]
        }
        else {
          currentPP <- newPP
        }
      }
      else {
        if (length(lambdaChangeTimes) > 0) {
          oldPP <- currentPP
          idx <- floor(runif(1, 0, length(lambdaChangeTimes))) + 
            1
          t <- lambdaChangeTimes
          v <- lambda
          lambdaChangeTimes <- lambdaChangeTimes[-idx]
          lambda <- lambda[-(idx + 1)]
          newPP <- prior(lambdaChangeTimes, lambda, muChangeTimes, 
                         mu, tMassExtinction, pMassExtinction) + likelihood(lambda, 
                                                                            lambdaChangeTimes, mu, muChangeTimes, tMassExtinction, 
                                                                            pMassExtinction)
          kPrime <- length(lambdaChangeTimes) + 1
          if (lambda.hyper.prior.form == "lognormal") {
            pr2 <- dunif(t[idx], 0, AGE) * dlnorm(v[idx + 
                                                      1], speciationRatePriorMean, speciationRatePriorStDev)
          }
          else if (lambda.hyper.prior.form == "normal") {
            pr2 <- dunif(t[idx], 0, AGE) * dnorm(v[idx + 
                                                     1], speciationRatePriorMean, speciationRatePriorStDev)
          }
          else if (lambda.hyper.prior.form == "gamma") {
            pr2 <- dunif(t[idx], 0, AGE) * dgamma(v[idx + 
                                                      1], speciationRatePriorMean, speciationRatePriorStDev)
          }
          if (is.finite(newPP - oldPP + log(pr2)) == 
              FALSE || log(runif(1, 0, 1)) > (newPP - oldPP + 
                                              log(pr2))) {
            lambdaChangeTimes <- t
            lambda <- v
          }
          else {
            currentPP <- newPP
          }
        }
      }
    }
    if (length(muChangeTimes) > 0) {
      idx <- floor(runif(1, 0, length(muChangeTimes))) + 
        1
      t <- muChangeTimes[idx]
      oldPP <- currentPP
      tPrime <- rnorm(1, t, deltaMuTime)
      muChangeTimes[idx] <- tPrime
      if (tPrime > 0 && tPrime < AGE) {
        newPP <- prior(lambdaChangeTimes, lambda, muChangeTimes, 
                       mu, tMassExtinction, pMassExtinction) + likelihood(lambda, 
                                                                          lambdaChangeTimes, mu, muChangeTimes, tMassExtinction, 
                                                                          pMassExtinction)
      }
      else {
        newPP <- -Inf
      }
      if (is.finite(newPP - oldPP) == FALSE || log(runif(1, 
                                                         0, 1)) > (newPP - oldPP)) {
        muChangeTimes[idx] <- t
      }
      else {
        currentPP <- newPP
        muTimeAccepted <- muTimeAccepted + 1
      }
      muTimeTried <- muTimeTried + 1
    }
    idx <- floor(runif(1, 0, length(mu))) + 1
    v <- mu[idx]
    oldPP <- currentPP
    u <- rnorm(1, 0, deltaMuValue)
    scalingFactor <- exp(u)
    vPrime <- v * scalingFactor
    mu[idx] <- vPrime
    lnHastingsratio <- log(scalingFactor)
    newPP <- prior(lambdaChangeTimes, lambda, muChangeTimes, 
                   mu, tMassExtinction, pMassExtinction) + likelihood(lambda, 
                                                                      lambdaChangeTimes, mu, muChangeTimes, tMassExtinction, 
                                                                      pMassExtinction)
    if (is.finite(newPP - oldPP + lnHastingsratio) == FALSE || 
        log(runif(1, 0, 1)) > (newPP - oldPP + lnHastingsratio)) {
      mu[idx] <- v
    }
    else {
      currentPP <- newPP
      muValueAccepted <- muValueAccepted + 1
    }
    muValueTried <- muValueTried + 1
    if (estimateNumberRateChanges == TRUE) {
      u <- runif(1, 0, 1)
      if (u > 0.5) {
        oldPP <- currentPP
        t <- runif(1, 0, AGE)
        if (mu.hyper.prior.form == "lognormal") {
          v <- rlnorm(1, extinctionRatePriorMean, extinctionRatePriorStDev)
        }
        else if (mu.hyper.prior.form == "normal") {
          v <- rnorm(1, extinctionRatePriorMean, extinctionRatePriorStDev)
        }
        else if (mu.hyper.prior.form == "gamma") {
          v <- rgamma(1, extinctionRatePriorMean, extinctionRatePriorStDev)
        }
        muChangeTimes[length(muChangeTimes) + 1] <- t
        mu[length(mu) + 1] <- v
        newPP <- prior(lambdaChangeTimes, lambda, muChangeTimes, 
                       mu, tMassExtinction, pMassExtinction) + likelihood(lambda, 
                                                                          lambdaChangeTimes, mu, muChangeTimes, tMassExtinction, 
                                                                          pMassExtinction)
        kPrime <- length(muChangeTimes)
        if (mu.hyper.prior.form == "lognormal") {
          pr <- 1/(dunif(t, 0, AGE) * dlnorm(v, extinctionRatePriorMean, 
                                             extinctionRatePriorStDev))
        }
        else if (mu.hyper.prior.form == "normal") {
          pr <- 1/(dunif(t, 0, AGE) * dnorm(v, extinctionRatePriorMean, 
                                            extinctionRatePriorStDev))
        }
        else if (mu.hyper.prior.form == "gamma") {
          pr <- 1/(dunif(t, 0, AGE) * dgamma(v, extinctionRatePriorMean, 
                                             extinctionRatePriorStDev))
        }
        if (is.finite(newPP - oldPP + log(pr)) == FALSE || 
            log(runif(1, 0, 1)) > (newPP - oldPP + log(pr))) {
          muChangeTimes <- muChangeTimes[-length(muChangeTimes)]
          mu <- mu[-length(mu)]
        }
        else {
          currentPP <- newPP
        }
      }
      else {
        if (length(muChangeTimes) > 0) {
          oldPP <- currentPP
          idx <- floor(runif(1, 0, length(muChangeTimes))) + 
            1
          t <- muChangeTimes
          v <- mu
          muChangeTimes <- muChangeTimes[-idx]
          mu <- mu[-(idx + 1)]
          newPP <- prior(lambdaChangeTimes, lambda, muChangeTimes, 
                         mu, tMassExtinction, pMassExtinction) + likelihood(lambda, 
                                                                            lambdaChangeTimes, mu, muChangeTimes, tMassExtinction, 
                                                                            pMassExtinction)
          kPrime <- length(muChangeTimes) + 1
          if (mu.hyper.prior.form == "lognormal") {
            pr2 <- dunif(t[idx], 0, AGE) * dlnorm(v[idx + 
                                                      1], extinctionRatePriorMean, extinctionRatePriorStDev)
          }
          else if (mu.hyper.prior.form == "normal") {
            pr2 <- dunif(t[idx], 0, AGE) * dnorm(v[idx + 
                                                     1], extinctionRatePriorMean, extinctionRatePriorStDev)
          }
          else if (mu.hyper.prior.form == "gamma") {
            pr2 <- dunif(t[idx], 0, AGE) * dgamma(v[idx + 
                                                      1], extinctionRatePriorMean, extinctionRatePriorStDev)
          }
          if (is.finite(newPP - oldPP + log(pr2)) == 
              FALSE || log(runif(1, 0, 1)) > (newPP - oldPP + 
                                              log(pr2))) {
            muChangeTimes <- t
            mu <- v
          }
          else {
            currentPP <- newPP
          }
        }
      }
    }
    if (estimateMassExtinctionTimes == TRUE && length(tMassExtinction) > 
        0) {
      idx <- floor(runif(1, 0, length(tMassExtinction))) + 
        1
      t <- tMassExtinction[idx]
      oldPP <- currentPP
      tPrime <- rnorm(1, t, deltaMassExtinctionTime)
      tMassExtinction[idx] <- tPrime
      if (tPrime > 0 && tPrime < AGE) {
        newPP <- prior(lambdaChangeTimes, lambda, muChangeTimes, 
                       mu, tMassExtinction, pMassExtinction) + likelihood(lambda, 
                                                                          lambdaChangeTimes, mu, muChangeTimes, tMassExtinction, 
                                                                          pMassExtinction)
      }
      else {
        newPP <- -Inf
      }
      if (is.finite(newPP - oldPP) == FALSE || log(runif(1, 
                                                         0, 1)) > (newPP - oldPP)) {
        tMassExtinction[idx] <- t
      }
      else {
        currentPP <- newPP
        massExtinctionTimeAccepted <- massExtinctionTimeAccepted + 
          1
      }
      massExtinctionTimeTried <- massExtinctionTimeTried + 
        1
    }
    if (length(pMassExtinction) > 0) {
      idx <- floor(runif(1, 0, length(pMassExtinction))) + 
        1
      v <- pMassExtinction[idx]
      oldPP <- currentPP
      vPrime <- rnorm(1, v, deltaMassExtinctionValue)
      pMassExtinction[idx] <- vPrime
      if (vPrime > 0 && vPrime < 1) {
        newPP <- prior(lambdaChangeTimes, lambda, muChangeTimes, 
                       mu, tMassExtinction, pMassExtinction) + likelihood(lambda, 
                                                                          lambdaChangeTimes, mu, muChangeTimes, tMassExtinction, 
                                                                          pMassExtinction)
      }
      else {
        newPP <- -Inf
      }
      if (is.finite(newPP - oldPP) == FALSE || log(runif(1, 
                                                         0, 1)) > (newPP - oldPP)) {
        pMassExtinction[idx] <- v
      }
      else {
        currentPP <- newPP
        massExtinctionValueAccepted <- massExtinctionValueAccepted + 
          1
      }
      massExtinctionValueTried <- massExtinctionValueTried + 
        1
    }
    if (estimateNumberMassExtinctions == TRUE) {
      u <- runif(1, 0, 1)
      if (u > 0.5) {
        oldPP <- currentPP
        t <- runif(1, 0, AGE)
        v <- rbeta(1, pMassExtinctionPriorShape1, pMassExtinctionPriorShape2)
        tMassExtinction[length(tMassExtinction) + 1] <- t
        pMassExtinction[length(pMassExtinction) + 1] <- v
        newPP <- prior(lambdaChangeTimes, lambda, muChangeTimes, 
                       mu, tMassExtinction, pMassExtinction) + likelihood(lambda, 
                                                                          lambdaChangeTimes, mu, muChangeTimes, tMassExtinction, 
                                                                          pMassExtinction)
        kPrime <- length(pMassExtinction)
        pr <- 1/(dunif(t, 0, AGE) * dbeta(v, pMassExtinctionPriorShape1, 
                                          pMassExtinctionPriorShape2))
        if (is.finite(newPP - oldPP + log(pr)) == FALSE || 
            log(runif(1, 0, 1)) > (newPP - oldPP + log(pr))) {
          tMassExtinction <- tMassExtinction[-length(tMassExtinction)]
          pMassExtinction <- pMassExtinction[-length(pMassExtinction)]
        }
        else {
          currentPP <- newPP
        }
      }
      else {
        if (length(tMassExtinction) > 0) {
          oldPP <- currentPP
          idx <- floor(runif(1, 0, length(tMassExtinction))) + 
            1
          t <- tMassExtinction
          v <- pMassExtinction
          tMassExtinction <- tMassExtinction[-idx]
          pMassExtinction <- pMassExtinction[-idx]
          newPP <- prior(lambdaChangeTimes, lambda, muChangeTimes, 
                         mu, tMassExtinction, pMassExtinction) + likelihood(lambda, 
                                                                            lambdaChangeTimes, mu, muChangeTimes, tMassExtinction, 
                                                                            pMassExtinction)
          kPrime <- length(tMassExtinction) + 1
          pr2 <- dunif(t[idx], 0, AGE) * dbeta(v[idx], 
                                               pMassExtinctionPriorShape1, pMassExtinctionPriorShape2)
          if (is.finite(newPP - oldPP + log(pr2)) == 
              FALSE || log(runif(1, 0, 1)) > (newPP - oldPP + 
                                              log(pr2))) {
            tMassExtinction <- t
            pMassExtinction <- v
          }
          else {
            currentPP <- newPP
          }
        }
      }
    }
    if (ADAPTIVE == TRUE & SAMPLE == FALSE) {
      if (i%%OPTIMIZATION_FREQUENCY == 0) {
        if (lambdaTimeTried > 0) {
          rate <- lambdaTimeAccepted/lambdaTimeTried
          if (rate > 0.44) {
            deltaLambdaTime <- deltaLambdaTime * (1 + 
                                                    ((rate - 0.44)/0.56))
          }
          else {
            deltaLambdaTime <- deltaLambdaTime/(2 - rate/0.44)
          }
          lambdaTimeTried <- 0
          lambdaTimeAccepted <- 0
        }
        rate <- lambdaValueAccepted/lambdaValueTried
        if (rate > 0.44) {
          deltaLambdaValue <- deltaLambdaValue * (1 + 
                                                    ((rate - 0.44)/0.56))
        }
        else {
          deltaLambdaValue <- deltaLambdaValue/(2 - rate/0.44)
        }
        lambdaValueTried <- 0
        lambdaValueAccepted <- 0
        if (muTimeTried > 0) {
          rate <- muTimeAccepted/muTimeTried
          if (rate > 0.44) {
            deltaMuTime <- deltaMuTime * (1 + ((rate - 
                                                  0.44)/0.56))
          }
          else {
            deltaMuTime <- deltaMuTime/(2 - rate/0.44)
          }
          muTimeTried <- 0
          muTimeAccepted <- 0
        }
        rate <- muValueAccepted/muValueTried
        if (rate > 0.44) {
          deltaMuValue <- deltaMuValue * (1 + ((rate - 
                                                  0.44)/0.56))
        }
        else {
          deltaMuValue <- deltaMuValue/(2 - rate/0.44)
        }
        muValueTried <- 0
        muValueAccepted <- 0
        if (massExtinctionTimeTried > 0) {
          rate <- massExtinctionTimeAccepted/massExtinctionTimeTried
          if (rate > 0.44) {
            deltaMassExtinctionTime <- deltaMassExtinctionTime * 
              (1 + ((rate - 0.44)/0.56))
          }
          else {
            deltaMassExtinctionTime <- deltaMassExtinctionTime/(2 - 
                                                                  rate/0.44)
          }
          massExtinctionTimeTried <- 0
          massExtinctionTimeAccepted <- 0
        }
        if (massExtinctionValueTried > 0) {
          rate <- massExtinctionValueAccepted/massExtinctionValueTried
          if (rate > 0.44) {
            deltaMassExtinctionValue <- deltaMassExtinctionValue * 
              (1 + ((rate - 0.44)/0.56))
          }
          else {
            deltaMassExtinctionValue <- deltaMassExtinctionValue/(2 - 
                                                                    rate/0.44)
          }
          massExtinctionValueTried <- 0
          massExtinctionValueAccepted <- 0
        }
      }
    }
    if (i%%THINNING == 0 & SAMPLE) {
      sampleIndex <- length(kLambda) + 1
      kLambda[sampleIndex] <- length(lambda)
      speciationRateValues[[sampleIndex]] <- lambda
      speciationRateChangeTimes[[sampleIndex]] <- lambdaChangeTimes
      kMu[sampleIndex] <- length(mu)
      extinctionRateValues[[sampleIndex]] <- mu
      extinctionRateChangeTimes[[sampleIndex]] <- muChangeTimes
      kMassExtinction[sampleIndex] <- length(pMassExtinction)
      survivalProbability[[sampleIndex]] <- pMassExtinction
      massExtinctionTime[[sampleIndex]] <- tMassExtinction
      posterior[sampleIndex] <- currentPP
      cat(lambda, sep = "\t", file = "SpeciationRates.txt", 
          append = TRUE)
      cat("\n", file = "SpeciationRates.txt", append = TRUE)
      cat(lambdaChangeTimes, sep = "\t", file = "SpeciationRateChanges.txt", 
          append = TRUE)
      cat("\n", file = "SpeciationRateChanges.txt", append = TRUE)
      cat(mu, sep = "\t", file = "ExtinctionRates.txt", 
          append = TRUE)
      cat("\n", file = "ExtinctionRates.txt", append = TRUE)
      cat(muChangeTimes, sep = "\t", file = "ExtinctionRateChanges.txt", 
          append = TRUE)
      cat("\n", file = "ExtinctionRateChanges.txt", append = TRUE)
      cat(pMassExtinction, sep = "\t", file = "SurvivalProbabilities.txt", 
          append = TRUE)
      cat("\n", file = "SurvivalProbabilities.txt", append = TRUE)
      cat(tMassExtinction, sep = "\t", file = "MassExtinctionTimes.txt", 
          append = TRUE)
      cat("\n", file = "MassExtinctionTimes.txt", append = TRUE)
      nSamples <- sampleIndex
      write.table(list(Iteration = nSamples * THINNING, 
                       posterior = posterior[nSamples], NumSpeciation = kLambda[nSamples], 
                       numExtinction = kMu[nSamples], numMassExtinctions = kMassExtinction[nSamples]), 
                  "samples_numCategories.txt", sep = "\t", row.names = FALSE, 
                  append = TRUE, quote = FALSE, col.names = FALSE)
    }
    if (i%%CONVERGENCE_FREQUENCY == 0 & SAMPLE & length(posterior) > 
        2) {
      samples <- length(posterior)
      if (estimateNumberRateChanges == TRUE) {
        spectralDensityLambda <- spectrum0.ar(kLambda)$spec
        spectralDensityMu <- spectrum0.ar(kMu)$spec
        if (spectralDensityLambda > 0 && spectralDensityMu > 
            0) {
          ess_kLambda <- (var(kLambda) * samples/spectralDensityLambda)
          ess_kMu <- (var(kMu) * samples/spectralDensityMu)
        }
        else {
          ess_kLambda <- ess_kMu <- Inf
        }
      }
      else {
        ess_kLambda <- ess_kMu <- Inf
      }
      spec <- c()
      for (j in 1:length(speciationRateValues)) {
        spec[j] <- speciationRateValues[[j]][1]
      }
      spectralDensity <- spectrum0.ar(spec)$spec
      if (spectralDensity > 0) {
        ess_spec <- (var(spec) * samples/spectralDensity)
      }
      else {
        ess_spec <- Inf
      }
      ext <- c()
      for (j in 1:length(extinctionRateValues)) {
        ext[j] <- extinctionRateValues[[j]][1]
      }
      spectralDensity <- spectrum0.ar(ext)$spec
      if (spectralDensity > 0) {
        ess_ext <- (var(ext) * samples/spectralDensity)
      }
      else {
        ess_ext <- Inf
      }
      if (estimateNumberMassExtinctions == TRUE) {
        ess_kMassExtinctions <- (var(kMassExtinction) * 
                                   samples/spectrum0.ar(kMassExtinction)$spec)
      }
      else {
        ess_kMassExtinctions <- Inf
      }
      min.ess <- min(c(ess_kLambda, ess_kMu, ess_spec, 
                       ess_ext, ess_kMassExtinctions))
    }
    if (SAMPLE) {
      iteration.progress <- i/MAX_ITERATIONS
      time.progress <- ((Sys.time() - startTime)/MAX_TIME)
      ess.progress <- min.ess/MIN_ESS
      max.progress <- max(c(iteration.progress, time.progress, 
                            ess.progress))
      if (verbose) {
        setTxtProgressBar(bar, pmin(max.progress, 1))
      }
      if (max.progress >= 1) {
        finished <- TRUE
      }
    }
    else {
      if (verbose) {
        setTxtProgressBar(bar, pmin(1, i/burn))
      }
    }
    i <- i + 1
    if (i == burn & !SAMPLE) {
      startTime <- Sys.time()
      i <- 1
      SAMPLE <- TRUE
      if (verbose) {
        cat("\n\nRunning the chain ... \n")
        cat("0--------25--------50--------75--------100\n")
        bar <- txtProgressBar(style = 1, width = 42)
      }
    }
  }
  cat("\n")
  setwd(file.path(orgDir))
}
