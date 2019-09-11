#' Plot of the error bars for the effects
#'
#' @param results.df
#' @param effect
#'
#' @return
#' @export
#'
#' @examples
effectsPlot = function(results.df, effect = c("locus", "profile", "dye"),
                       xlab = NULL, ylab = NULL){
  effect = match.arg(effect)

  effectPrefix = if(effect == "locus"){
    "alpha"
  }else if(effect == "profile"){
    "beta"
  }else{
    "gamma"
  }

  if(is.null(xlab)){
    xlab = stringr::str_to_title(effect)
  }

  if(is.null(ylab)){
    ylab = "Effect size"
  }
  effect = paste0(effectPrefix, ".", effect)

  stats.df = data.frame(
    lb = results.df %>%
    select(starts_with(effect)) %>%
    summarise_all(quantile, prob = c(0.025)) %>%
    unlist(),

    med = results.df %>%
      select(starts_with(effect)) %>%
      summarise_all(quantile, prob = c(0.5)) %>%
      unlist(),

    ub = results.df %>%
      select(starts_with(effect)) %>%
      summarise_all(quantile, prob = c(0.975)) %>%
      unlist()
  )

  stats.df = stats.df %>%
    mutate(x = 1:nrow(stats.df))


  p = stats.df %>%
    ggplot(aes(x = x, y = med)) +
    geom_point(aes(col = "red"), show.legend = FALSE) +
    geom_errorbar(aes(ymin = lb, ymax = ub)) +
    geom_hline(aes(yintercept = 0)) +
    ylab(ylab) +
    xlab(xlab) +
    theme(text = element_text(size = 20), aspect.ratio = 1,
          plot.background = element_rect(color="black"))

  p
}
