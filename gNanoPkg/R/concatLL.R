concatLL = function(..., labels = NULL){
  model.ll.list = list(...)

  numModels = length(model.ll.list)

  if(numModels <= 1){
    stop("There must be more than one model to concatenate")
  }

  if(is.null(labels)){
    labels = paste0("Model ", 1:numModels)
  }

  #names(model.ll.list) = labels

  all = model.ll.list %>%
    reduce(left_join, by = c("rep" = "rep")) %>%
    select(-rep)

  names(all) = labels

  all = all %>%
    pivot_longer(everything(), names_to = "model", values_to = "logLik")


  return(all)
}
