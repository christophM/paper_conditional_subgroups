#' Compute grouped feature importance
#'
#' Partitions the data, computes the feature importance per group and aggregates the results again
#'
#' @param pred iml::Predictor
#' @param group integer vector with group indicators; length should be nrow(pred$data$X)
#' @param loss loss function string, see ?iml::FeatureImp
#' @param fname character, feature name for which to compute the feature importance
#' @return data.frame with feature importances

grouped_pfi_feat = function(pred, group, loss, fname) {
  if(all(is.na(group))) group = rep("", times = length(group)) 
  ugroups = unique(group)
  dat = pred$data$X
  target = pred$data$y[[1]]
  res = lapply(ugroups, function(g) {
    gi = which(group == g)
    d = data.frame(dat[gi,])
    y = target[gi]
    pred_sub = iml::Predictor$new(predict.fun = pred$predict, data = d, y = y)
    loss_orig = loss(y, pred_sub$predict(d)[[1]])
    d_perm = d
    # The sample function works differently when only one number is given
    if (nrow(d) > 1) d_perm[fname] = myshuffle(d_perm[[fname]])
    loss_perm = loss(y, pred$predict(d_perm)[[1]])
    data.frame(fname = fname,
               group = g,
               n = nrow(d),
               importance = mean(loss_perm - loss_orig))
  })
  rbindlist(res)
}

# Makes sure that a value does not randomly get assigned the same value again
myshuffle = function(x){
  if(is.factor(x)){
    return(sample(x))
  }
  n = length(x)
  ord = order(x)
  # assures that order is different each time
  x = sample(x)
  ord2 = order(x)
  while(sum(ord == ord2) > 0) {
    x = sample(x)
    ord2 = order(x)
  }
  x
}
#' cPFI computation
#'
#' @param pred iml::Predictor object
#' @param loss loss function loss(ytrue, ypred)
#' @param conds list of Conditional
#' @param repetitions Number of repeated permutations
#' @return data.frame with cPFI
grouped_pfi = function(pred, loss, conds, repetitions = 1) {
  res <- lapply(1:repetitions, function(i) {
    res <- lapply(names(conds), function(fname) {
      cond <- conds[[fname]]
      group <- cond$cnode(pred$data$X)$node
      res <- grouped_pfi_feat(pred, group, loss, fname)
      imp <- weighted.mean(res$importance, w = res$n)
      condition = paste(get_split_vars(cond$model), collapse = ",")
      data.frame(feature = fname, importance = imp, condition = condition) 
    })
    res <- rbindlist(res)
    res$i <- i
    res
  })
  result <- rbindlist(res)
  result[, list(
                "importance" = median(importance),
                "importance.05" = quantile(importance, probs = 0.05),
                "importance.95" = quantile(importance, probs = 0.95),
                "condition" = unique(condition)
                ), by = list(feature)]
}


#' subgroup PFIs computation
#'
#' @param pred iml::Predictor object
#' @param loss loss function loss(ytrue, ypred)
#' @param conds list of Conditional
#' @param repetitions Number of repeated permutations
#' @return data.frame with all subgroup PFIs
grouped_pfi_groupwise = function(pred, loss, conds, repetitions = 1) {
  res <- lapply(1:repetitions, function(i) {
    res <- lapply(names(conds), function(fname) {
      cond <- conds[[fname]]
      group <- cond$cnode(pred$data$X)$.path
      res <- grouped_pfi_feat(pred, group, loss, fname)
      res$condition = paste(get_split_vars(cond$model), collapse = ",")
      res$feature = fname
      res
    })
    res <- rbindlist(res)
    res$i <- i
    res
  })
  result <- rbindlist(res)
  result[, list(
                "importance" = median(importance),
                "importance.05" = quantile(importance, probs = 0.05),
                "importance.95" = quantile(importance, probs = 0.95),
                "condition" = unique(condition)
                ), by = list(feature, group)]
}

#' Extract split features
#' 
#' @param ctree_model a partykit ctree
#' @return list of features
get_split_vars = function(ctree_model) {
  ids = partykit::nodeids(ctree_model, terminal = FALSE)
  nids = partykit::nodeids(ctree_model, terminal = TRUE)
  ids = setdiff(ids, nids)
  svars = sapply(ids, function(i) {
    node  = partykit::node_party(ctree_model[[i]])
    split <- partykit::split_node(node)
    ivar <- partykit::varid_split(split)
    svar <- unlist(names(ctree_model$data)[ivar])
})
  unique(svars)
}

