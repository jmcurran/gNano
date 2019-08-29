.dsn = function (x, xi = 0, omega = 1, alpha = 0, tau = 0, dp = NULL, log = FALSE){
  if (!is.null(dp)) {
    if (!missing(alpha))
      stop("You cannot set both 'dp' and component parameters")
    xi <- dp[1]
    omega <- dp[2]
    alpha <- dp[3]
    tau <- if (length(dp) > 3)
      dp[4]
    else 0
  }
  z <- (x - xi)/omega
  logN <- (-log(sqrt(2 * pi)) - logb(omega) - z^2/2)
  logS = rep(0, length(alpha))

  if (all(abs(alpha) < Inf))
    logS <- pnorm(tau * sqrt(1 + alpha^2) + alpha * z, log.p = TRUE)
  else{
    i = abs(alpha) < Inf
    logS = i
    logS[i] <- pnorm(tau[i] * sqrt(1 + alpha[i]^2) + alpha[i] * z[i], log.p = TRUE)
    logS[!i] <- log(as.numeric(sign(alpha[!i]) * z[!i] + tau[!i] > 0))
  }

  logPDF <- as.numeric(logN + logS - pnorm(tau, log.p = TRUE))
  logPDF <- replace(logPDF, abs(x) == Inf, -Inf)
  logPDF <- replace(logPDF, omega <= 0, NaN)

  if (log)
    logPDF
  else exp(logPDF)
}
