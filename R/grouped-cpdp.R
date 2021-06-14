#' Grouped cPDP
#'@param pred iml::Predictor
#'@param cond conditional model
#'@param fname character with name of feature
#'@return list with nodes data.frame and the cgPDP-plot 
grouped_pdp = function(pred, fname, cond, grid_size = 30, bplot = TRUE) {
  assert_numeric(grid_size)
  dat = pred$data$X
  nodes_df = cond$cnode(dat)
  nodes_df$feature = dat[[fname]]
  nodes_df = nodes_df %>%
    group_by(node, .path) %>%
    mutate(xmin= min(feature), xmax = max(feature))

  nodes = unique(nodes_df$node)
  res = lapply(nodes, function(node) {
    subs_ids = which(nodes_df$node == node)
    if (length(subs_ids) <= 1) return(data.frame())
    subs = dat[subs_ids, ]
    pred_sub = iml::Predictor$new(mod = pred$model, data = data.frame(subs))
    fe = FeatureEffect$new(pred_sub, feature = fname, grid.size = grid_size, method = "pdp")
    ret = fe$results
    ret$node = rep(nodes_df[subs_ids[1], ".path"][[1]], times = nrow(ret))
    bpstats = data.frame(stats = boxplot.stats(subs[[fname]])$stats,
                         feature = fname)
    bpstats$.value = fe$predict(bpstats$stats)
    outliers = subs[subs[[fname]] < bpstats$stats[1] |
                    subs[[fname]] > bpstats$stats[5],]
    outliers$node = ret$node[1]
    outliers$.value = fe$predict(outliers[[fname]])
    ret$in_hinges = ret[[fname]] >= bpstats$stats[1] &
                    ret[[fname]] <= bpstats$stats[5]
    ret$in_box = ret[[fname]] >= bpstats$stats[2] &
                 ret[[fname]] <= bpstats$stats[4]
    return(list("pdp" = ret, "outliers" = outliers, "bpstats" = bpstats))
  })
  pdp = rbindlist(lapply(res, function(x) x$pdp))
  outliers = rbindlist(lapply(res, function(x) x$outliers))
  bpstats = rbindlist(lapply(res, function(x) x$bpstats))
  if(bplot) {
    p = ggplot(data = pdp,
               aes_string(x = fname, y = ".value",
                          group = "node", color = "node")) + 
          geom_line(data = pdp[pdp$in_hinges,]) +
          geom_line(data = pdp[pdp$in_box,], size = 2) + 
          geom_point(data = outliers)
  } else {
    p = ggplot(data.frame(pdp),
               aes_string(x = fname, y = ".value",
                          group = "node", color = "node")) + 
          geom_line()
  }
 list(pdp = pdp, outliers = outliers, bpstats = bpstats, p = p)
}


grouped_pdp_factor = function(pred, fname, cond) {
  dat = pred$data$X
  nodes_df = cond$cnode(dat)
  nodes_df$feature = dat[[fname]]

  nodes = unique(nodes_df$node)
  res = lapply(nodes, function(node) {
    subs_ids = which(nodes_df$node == node)
    subs = dat[subs_ids, ]
    pred_sub = iml::Predictor$new(mod = pred$model, data = data.frame(subs))
    if (length(unique(subs[[fname]])) == 1) {
      ret = data.frame(pred$predict(subs)[[1]],
                       subs[[fname]])
      colnames(ret) = c(".value", fname) 
    } else {
      fe = FeatureEffect$new(pred_sub, feature = fname, method = "ice")
      ret = fe$results
    }
      ret$node = rep(nodes_df[subs_ids[1], ".path"][[1]], times = nrow(ret))
      ret$node_id = rep(nodes_df[subs_ids[1], "node"][[1]], times = nrow(ret))
      ret
  })
  res = rbindlist(res, fill = TRUE)
  p_group = ggplot(data.frame(res)) + 
    geom_boxplot(aes_string(x = fname, y = ".value")) + 
    facet_wrap("node", scales = "free_x")
  list(df = res, p = p_group)
}
