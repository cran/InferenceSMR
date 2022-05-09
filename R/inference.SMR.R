inference.SMR<-function (obs.death, normal = "log-smr", alpha = 0.05, contribution, 
    incid, cox, fuzz = 0.01, Poisson = FALSE, covnames) {
    expected = est.expDeath(contribution, incid, cox, fuzz, covnames)
    variance = var.expDeath(contribution, incid, cox, fuzz, Poisson, 
        covnames)
    results = list(expected = expected, obs.death = obs.death, 
        variance = variance, smr = smr)
    smr = obs.death/expected
    cat("\n\n\n********  INFERENCE ABOUT THE SMR  ********* \n\nObserved = ", 
        obs.death, " Expected = ", expected, "\nObs.var. = ", 
        obs.death, " Exp.var. = ", variance, "\nSMR = ", smr, 
        "\n\n", (1 - alpha) * 100, "% Confidence intervals with normality assumption at : \n\n")
    if ("smr" %in% normal) {
        smr.var = 1/expected^2 * obs.death + obs.death^2/expected^4 * 
            variance
        smr.ci = c(smr - qnorm(1 - alpha/2) * sqrt(smr.var), 
            smr + qnorm(1 - alpha/2) * sqrt(smr.var))
        cat("The SMR level : (", smr.ci, ")\n\n")
        results = c(results, list(smr.var = smr.var, smr.ci = smr.ci))
    }
    if ("log-smr" %in% normal) {
        logSMR.var = 1/obs.death + 1/expected^2 * variance
        logSMR.ci = c(smr * exp(-qnorm(1 - alpha/2) * sqrt(logSMR.var)), 
            smr * exp(qnorm(1 - alpha/2) * sqrt(logSMR.var)))
        cat("The log-SMR level : (", logSMR.ci, ")\n\n")
        results = c(results, list(logSMR.var = logSMR.var, logSMR.ci = logSMR.ci))
    }
    if ("root-smr" %in% normal) {
        rootSMR.var = 1/(4 * expected) + obs.death/(4 * expected^3) * 
            variance
        rootSMR.ci = c((sqrt(smr) - qnorm(1 - alpha/2) * sqrt(rootSMR.var))^2, 
            (sqrt(smr) + qnorm(1 - alpha/2) * sqrt(rootSMR.var))^2)
        cat("The root-SMR level : (", rootSMR.ci, ")\n\n")
        results = c(results, list(rootSMR.var = rootSMR.var, 
            rootSMR.ci = rootSMR.ci))
    }
    return(results)
}
