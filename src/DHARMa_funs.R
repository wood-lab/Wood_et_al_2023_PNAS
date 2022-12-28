getPt <- function(simulated, observed, alternative, plot = F){
  
  if(alternative == "greater") p = mean(simulated >= observed)
  if(alternative == "less") p = mean(simulated <= observed)
  if(alternative == "two.sided") p = min(min(mean(simulated <= observed), mean(simulated >= observed) ) * 2,1)
  
  if(plot == T){
    hist(simulated, xlim = range(simulated, observed), col = "lightgrey")
    abline(v = mean(simulated), col = 1, lwd = 2)
    abline(v = observed, col = "red", lwd = 2)
  }
  
  return(p)
}


plotQQunift <- function(simulatedOutput, ...){
  
  gap::qqunif(u = simulatedOutput$scaledResiduals, type = "unif", logscale = F, pch = 2,
              bty = "n",
              col = "black",
              xlim = c(0,1),
              ylim = c(0,1))
  
  # Ks test on scaled Residuals
  out <- list()
  scaledResiduals <- simulationOutput$scaledResiduals
  nSim <- simulationOutput$nSim
  
  out$ks.test <- suppressWarnings(ks.test(scaledResiduals, 
                                  "punif"))
  # DHARMa outlier test, modified
  
  outliers = sum(scaledResiduals <
                   (1 / (nSim + 1))) + sum(scaledResiduals >
                                             (1 - 1 / (nSim + 1)))
  
  outFreqH0 = 1/(nSim + 1) * 2
  out$outlier = binom.test(outliers, simulationOutput$nObs, p = outFreqH0, alternative = "two.sided")
  
  # DHARMa - based Dispersion test
  expectedVar = sd(simulationOutput$simulatedResponse)^2
  spread <- function(x) var(x - simulationOutput$fittedPredictedResponse)/expectedVar
  observed = as.numeric(spread(simulationOutput$observedResponse))
  simulated = apply(simulationOutput$simulatedResponse, 2, 
                    spread)
  
  
  p = getPt(simulated = simulated, observed = observed, alternative = "two.sided")
  out$dispersion$statistic = c(ratioObsSim = observed/mean(simulated))
  out$dispersion$p.value = p
  return(out)
}

plotResidualst <- function(simulationOutput, form = NULL, quantreg = NULL, rank = T, asFactor = NULL, smoothScatter = NULL, quantiles = c(0.25, 0.5, 0.75), ...){
  
  
  ##### Setup #####
  
  a <- list()
  a$xlab ="Model predictions (rank transformed)"
  a$ylab = "DHAMa residual"
  
 
  res = simulationOutput$scaledResiduals
   pred <- simulationOutput$fittedPredictedResponse
  ##### Rank transform and factor conversion#####
    
    if (rank == T){
      pred = rank(pred, ties.method = "average")
      pred = pred / max(pred)
      a$xlim =  c(0,1)
    }
    
    nuniq = length(unique(pred))
    ndata = length(pred)
    #if(is.null(asFactor)) asFactor = (nuniq == 1) | (nuniq < 10 & ndata / nuniq > 10)
    #if (asFactor) pred = factor(pred)
  #}
  
  ##### Residual scatter plots #####
  
  if(is.null(quantreg)) if (length(res) > 2000) quantreg = FALSE else quantreg = TRUE
  
  switchScatter = 10000
  if(is.null(smoothScatter)) if (length(res) > switchScatter) smoothScatter = TRUE else smoothScatter = FALSE
  
  blackcol = rgb(0,0,0, alpha = max(0.1, 1 - 3 * length(res) / switchScatter))
  
  
  # smooth scatter
  if (smoothScatter == TRUE) {
    defaultCol = ifelse(res == 0 | res == 1, 2,blackcol)
    graphics::smoothScatter(x = pred,
                            y = res,
                            ylim = c(0,1),
                            xlim = c(0,1),
                            ylab = a$ylab,
                            xlab = a$xlab,
                            axes = F,
                            colramp = colorRampPalette(c("white", "darkgrey"))
    )
    do.call(graphics::smoothScatter, append(list(x = pred, y = res , ylim = c(0,1), axes = FALSE, colramp = colorRampPalette(c("white", "darkgrey"))),
                                            a))
    points(pred[defaultCol == 2], res[defaultCol == 2], col = "red", cex = 0.5)
    
    axis(1)
    axis(2, at=c(0, quantiles, 1))
  }
  # normal plot
  else {
    defaultCol = ifelse(res == 0 | res == 1, 2,blackcol)
    defaultPch = ifelse(res == 0 | res == 1, 8,1)
    a$col = defaultCol #checkDots("col", defaultCol, ...)
    a$pch = defaultPch#checkDots("pch", defaultPch, ...)
    do.call("plot", append(list(res ~ pred, ylim = c(0,1), axes = FALSE), a))
    
    axis(1)
    axis(2, at=c(0, quantiles, 1))
  }
  
  ##### Quantile regressions #####
  
  main = ""# checkDots("main", ifelse(is.null(form), "Residual vs. predicted", "Residual vs. predictor"), ...)
  out = NULL
  
  if(is.numeric(pred)){
    if(quantreg == F){
      title(main = main, cex.main = 1)
      abline(h = quantiles, col = "black", lwd = 0.5, lty = 2)
      try({
        lines(smooth.spline(pred, res, df = 10), lty = 2, lwd = 2, col = "red")
        abline(h = 0.5, col = "red", lwd = 2)
      }, silent = T)
    }else{
      
      out = testQuantilest(scaledResiduals, fittedPredictedResponse, pred, quantiles = quantiles, plot = F)
      
      
      if(any(out$pvals < 0.05, na.rm = TRUE)){
        main = paste(main, "Quantile deviations detected (red curves)", sep ="\n")
        if(out$p.value <= 0.05){
          main = paste(main, "Combined adjusted quantile test significant", sep ="\n")
        } else {
          main = paste(main, "Combined adjusted quantile test n.s.", sep ="\n")
        }
        maincol = "red"
      } else {
        main = paste(main, "No significant problems detected", sep ="\n")
        maincol = "black"
      }
      title(main = main, cex.main = 0.8,
            col.main = maincol)
      
      for(i in 1:length(quantiles)){
        
        lineCol = ifelse(out$pvals[i] <= 0.05 & !(is.na(out$pvals[i])), "red", "black")
        filCol = ifelse(out$pvals[i] <= 0.05 & !(is.na(out$pvals[i])), "#FF000040", "#00000020")
        
        abline(h = quantiles[i], col = lineCol, lwd = 0.5, lty = 2)
        polygon(c(out$predictions$pred, rev(out$predictions$pred)),
                c(out$predictions[,2*i] - out$predictions[,2*i+1], rev(out$predictions[,2*i] + out$predictions[,2*i+1])),
                col = "#00000020", border = F)
        lines(out$predictions$pred, out$predictions[,2*i], col = lineCol, lwd = 2)
      }
      
      # legend("bottomright", c(paste("Quantile test: p=", round(out$p.value, digits = 5)), paste("Deviation ", ifelse(out$p.value < 0.05, "significant", "n.s."))), text.col = ifelse(out$p.value < 0.05, "red", "black" ), bty="n")
      
    }
  }
  invisible(out)
}

#x = 0.01
#x <= 0.05 & !(is.na(x))






testQuantilest <- function (simulationOutput,
                           predictor = NULL, 
                           quantiles = c(0.25, 0.5, 0.75), plot = T) 
{
  if (plot == F) {
    out = list()
#    out$data.name = deparse(substitute(simulationOutput))
#    simulationOutput = ensureDHARMa(simulationOutput, convert = T)
    res = simulationOutput$scaledResiduals
    pred = simulationOutput$fittedPredictedResponse
    dat = data.frame(res = scaledResiduals, 
                     pred = pred)
    quantileFits <- list()
    pval = rep(NA, length(quantiles))
    predictions = data.frame(pred = sort(dat$pred))
    predictions = cbind(predictions, matrix(ncol = 2 * length(quantiles), 
                                            nrow = nrow(dat)))
    for (i in 1:length(quantiles)) {
      datTemp = dat
      datTemp$res = datTemp$res - quantiles[i]
      dimSmooth = min(length(unique(datTemp$pred)), 10)
      quantResult = try(capture.output(quantileFits[[i]] <- qgam::qgam(res ~ 
                                                                         s(pred, k = dimSmooth), data = datTemp, qu = quantiles[i])), 
                        silent = T)
      if (inherits(quantResult, "try-error")) {
        message("Unable to calculate quantile regression for quantile ", 
                quantiles[i], ". Possibly to few (unique) data points / predictions. Will be ommited in plots and significance calculations.")
      }
      else {
        x = summary(quantileFits[[i]])
        pval[i] = min(p.adjust(c(x$p.table[1, 4], x$s.table[1, 
                                                            4]), method = "BH"))
        quantPre = predict(quantileFits[[i]], newdata = predictions, 
                           se = T)
        predictions[, 2 * i] = quantPre$fit + quantiles[i]
        predictions[, 2 * i + 1] = quantPre$se.fit
      }
    }
    out$method = "Test for location of quantiles via qgam"
    out$alternative = "both"
    out$pvals = pval
    out$p.value = min(p.adjust(pval, method = "BH"))
    out$predictions = predictions
    out$qgamFits = quantileFits
    class(out) = "htest"
  }
  else if (plot == T) {
    out <- plotResiduals(simulationOutput = simulationOutput, 
                         form = predictor, quantiles = quantiles, quantreg = TRUE)
  }
  return(out)
}


getQuantilet <- function (simulations, observed, integerResponse, method = c("PIT", 
                                                                                "traditional"), rotation = NULL) 
{
  method = match.arg(method)
  n = length(observed)
  if (nrow(simulations) != n) 
    stop("DHARMa::getquantile: wrong dimension of simulations")
  nSim = ncol(simulations)
  if (method == "traditional") {
    if (!is.null(rotation)) 
      stop("rotation can only be used with PIT residuals")
    if (integerResponse == F) {
      if (any(duplicated(observed))) 
        message("Model family was recognized or set as continuous, but duplicate values were detected in the response. Consider if you are fitting an appropriate model.")
      values = as.vector(simulations)[duplicated(as.vector(simulations))]
      if (length(values) > 0) {
        if (all(values%%1 == 0)) {
          integerResponse = T
          message("Model family was recognized or set as continuous, but duplicate values were detected in the simulation - changing to integer residuals (see ?simulateResiduals for details)")
        }
        else {
          message("Duplicate non-integer values found in the simulation. If this is because you are fitting a non-inter valued discrete response model, note that DHARMa does not perform appropriate randomization for such cases.")
        }
      }
    }
    scaledResiduals = rep(NA, n)
    for (i in 1:n) {
      if (integerResponse == T) {
        scaledResiduals[i] <- DHARMa.ecdf(simulations[i,] + 
                                            runif(nSim, -0.5, 0.5))(observed[i] + runif(1, 
                                                                                        -0.5, 0.5))
      }
      else {
        scaledResiduals[i] <- DHARMa.ecdf(simulations[i, 
        ])(observed[i])
      }
    }
  }
  else {
    if (!is.null(rotation)) {
      if (is.character(rotation) && rotation == "estimated") {
        covar = Matrix::nearPD(cov(t(simulations)))$mat
        L = t(as.matrix(Matrix::chol(covar)))
      }
      else if (is.matrix(rotation)) 
        L <- t(chol(rotation))
      else stop("DHARMa::getQuantile - wrong argument to rotation parameter")
      observed <- solve(L, observed)
      simulations = apply(simulations, 2, function(a) solve(L, 
                                                            a))
    }
    scaledResiduals = rep(NA, n)
    for (i in 1:n) {
      minSim <- mean(simulations[i, ] < observed[i])
      maxSim <- mean(simulations[i, ] <= observed[i])
      if (minSim == maxSim) 
        scaledResiduals[i] = minSim
      else scaledResiduals[i] = runif(1, minSim, maxSim)
    }
  }
  return(scaledResiduals)
}

simulateResidualst <- function(y, n = 250, refit = F, integerResponse = T, plot = F, seed = 123, method = c("PIT", "traditional"), rotation = NULL,simulations, fitted, ...){
  
  ######## general assertions and startup calculations ##########
  
  randomState <-getRandomState(seed)
  on.exit({randomState$restoreCurrent()})
  ptm <- proc.time()
  
  ####### extract model info ############
  
  out = list()
  
  family = "nbinom2"
  out$fittedModel = NULL
  out$modelClass = "glmmTMB"
  
  out$additionalParameters = list(...)
  
  out$nObs = length(y)
  out$nSim = n
  out$refit = refit
  out$observedResponse = y
  
  # this is to check for problem #325
  #  if (out$nObs < length(out$observedResponse)) stop("DHARMA::simulateResiduals: nobs(model) < nrow(model.frame). A possible reason is that you have observation with zero prior weights (ir binomial with n=0) in your data. Calculating residuals in this case wouldn't be sensible. Please remove zero-weight observations from your data and re-fit your model! If you believe that this is not the reason, please file an issue under https://github.com/florianhartig/DHARMa/issues")
  
  # if(is.null(integerResponse)){
  #   if (family$family %in% c("binomial", "poisson", "quasibinomial", "quasipoisson", "Negative Binom", "nbinom2", "nbinom1", "genpois", "compois", "truncated_poisson", "truncated_nbinom2", "truncated_nbinom1", "betabinomial", "Poisson", "Tpoisson", "COMPoisson", "negbin", "Tnegbin") | grepl("Negative Binomial",family$family) ) integerResponse = TRUE
  #   else integerResponse = FALSE
  # }
  out$integerResponse = integerResponse
  
  out$problems = list()
  
  out$fittedPredictedResponse = fitted
  #out$fittedFixedEffects = getFixedEffects(fittedModel)
  #out$fittedResiduals = getResiduals(fittedModel)
  
  ######## refit = F ##################
  
  if (refit == FALSE){
    
    out$method = method
    
    #out$simulatedResponse = getSimulations(fittedModel, nsim = n, type = "normal")
    out$simulatedResponse <- simulations
    #checkSimulations(out$simulatedResponse, out$nObs, out$nSim)
    
    out$scaledResiduals = getQuantile(simulations = out$simulatedResponse , observed = out$observedResponse , integerResponse = integerResponse, method = method, rotation = rotation)
  }
  ########### Wrapup ############
  
  out$time = proc.time() - ptm
  out$randomState = randomState
  
  class(out) = "DHARMa"
  
  if(plot == TRUE) plot(out)
  
  return(out)
}


plot.DHARMat <- function(simulationOutput, title = "DHARMa residual", ...){
  oldpar <- par(mfrow = c(1,2), las = 1,xpd = F)
  on.exit(par(oldpar))
  
  out1 <- plotQQunift(simulationOutput)
  out2 <-plotResidualst(simulationOutput, ...)
  return(list(qqout = out1, residout = out2))
  
  #mtext(title, outer = T)
}

