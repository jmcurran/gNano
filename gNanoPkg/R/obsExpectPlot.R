obsExpectPlot = function(results.df, responseDist = c("g", "ln", "sn")){
  gNano.df = readData()

  responseDist = match.arg(responseDist)

  expected = results.df %>%
    select(starts_with("pred")) %>%
    summarise_all(mean) %>%
    unlist()

  plot.df = if(responseDist == "g"){
    data.frame(observed = gNano.df$obs, expected = expected)
  }else{
    data.frame(observed = log10(gNano.df$obs), expected = log10(exp(expected)))
  }

  p = plot.df %>%
    ggplot(aes(x = expected, y = observed)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1)
  p
}
