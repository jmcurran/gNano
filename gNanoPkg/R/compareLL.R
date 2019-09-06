#' Compare log-likelihoods
#'
#' @param \dots two or more data frames containing log likelihoods which come from calcLogLik
#' @param lables an optional vector of label for each model
#'
#' @return
#' @export
#'
#' @examples
compareLL = function(..., labels = NULL){

  # if(!missing(model1LL.df)){
  #   model1lab = ifelse(is.null(model1lab), deparse(substitute(model1LL.df)), model1lab)
  # }
  #
  # if(!missing(model2LL.df)){
  #   model2lab = ifelse(is.null(model2lab), deparse(substitute(model2LL.df)), model2lab)
  # }
  #
  #
  # plot.df = left_join(model1LL.df, model2LL.df, by = c("rep" = "rep")) %>%
  #   pivot_longer(cols = c(ll.x, ll.y), names_to = "model", values_to = "logLik") %>%
  #   mutate(model = case_when(
  #     grepl("^ll\\.x$", model) ~ model1lab,
  #     grepl("^ll\\.y$", model) ~ model2lab
  #   ))

  plot.df = concatLL(..., labels = labels)

  p = plot.df %>%
    ggplot(aes(x = logLik, color = model)) + geom_density()
  p
}
