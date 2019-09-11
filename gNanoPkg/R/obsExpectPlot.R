#' Plot the Expected versus Observed or Pred-Res
#'
#' @param results.df A data frame from the loadResults function
#' @param responseDist one of \code{"g"},\code{"ln"}, or \code{"sn"}
#' @param predRes if \code{TRUE} then a pred-res (residuals vs. fitted) plot is displayed
#'
#' @return a list with the ggplot object and the data used to create the plot
#' @export
obsExpectPlot = function(results.df, responseDist = c("g", "ln", "sn"),
                         predRes = FALSE){
  gNano.df = readData()

  responseDist = match.arg(responseDist)

  plot.df = data.frame(
    mean = results.df %>%
      select(starts_with("pred")) %>%
      summarise_all(mean) %>%
      unlist(),

    median = results.df %>%
      select(starts_with("pred")) %>%
      summarise_all(median) %>%
      unlist(),

    lwr = results.df %>%
      select(starts_with("pred")) %>%
      summarise_all(function(x)quantile(x, probs = 0.025)) %>%
      unlist(),

    upr = results.df %>%
      select(starts_with("pred")) %>%
      summarise_all(function(x)quantile(x, probs = 0.975)) %>%
      unlist()
  )

  o = order(gNano.df$obs)
  plot.df = plot.df[o,]

  plot.df = plot.df %>%
    mutate(observed = log10(gNano.df$obs[o]))


  plot.df = if(responseDist == "g"){
    plot.df %>%
      mutate_at(vars(-observed), function(col)log10(col))
  }else{
    plot.df %>%
      mutate_at(vars(-observed), function(col)log10(exp(col)))
  }

  plot.df = plot.df %>%
    mutate(residuals = observed - mean)

  p = if(predRes){
        plot.df %>%
          ggplot(aes(x = mean, y = residuals)) +
          geom_point() +
          geom_hline(linetype = "dashed", colour = "grey", alpha = 0.2,
                     yintercept = 0) +
          geom_smooth() +
          xlab("Fitted (posterior mean)")
      }else{
        plot.df %>%
          ggplot(aes(y = mean, x = observed)) +
            geom_point() +
            geom_abline(intercept = 0, slope = 1) +
            geom_smooth(col = "blue") +
            geom_smooth(aes(y = median), col = "red") +
            geom_smooth(aes(y = lwr), col = "green") +
            geom_smooth(aes(y = upr), col = "green") +
            xlim(1, 5) +
            ylim(1, 5) +
            xlab(expression(log[10]~(observed~height))) +
            ylab(expression(log[10]~(expected~height))) +
            theme(text = element_text(size = 20),
                aspect.ratio = 1,
                plot.background = element_rect(color="black"))
      }

  print(p)
  invisible(list(p = p, plot.df = plot.df))
}
