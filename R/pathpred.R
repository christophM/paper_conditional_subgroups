# Return the paths of a ctree for each training data point
pathpred <- function(object, ...) {
  ## coerce to "party" object if necessary
  if (!inherits(object, "party")) object <- partykit::as.party(object)

  ## get rules for each node
  rls <- list.rules.party(object)

  ## get predicted node and select corresponding rule
  rules <- rls[as.character(predict(object, type = "node", ...))]
  rules <- gsub("&", "&\n", rules)

  return(rules)
}


