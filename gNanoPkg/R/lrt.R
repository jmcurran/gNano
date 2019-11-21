#' Title
#'
#' @param model1LL
#' @param model2LL
#' @param df
#'
#' @return
#' @export
#'
#' @examples
lrt = function(model1LL, model2LL, df){
  l1 = model1LL %>% pull(ll) %>% mean()
  l2 = model2LL %>% pull(ll) %>% mean()

  X2 = -2*(l1 - l2)
  p.val = 1 - pchisq(X2, df)

  cat(paste0("Test stat: ", round(X2, 4), "\n"))
  cat(paste0("df:", df, "\n"))
  cat(paste0("Critical Value: ", round(qchisq(0.95, df), 3), "\n"))
  cat(paste0("P: ", round(p.val, 4), "\n"))

  invisible(p.val)

}
