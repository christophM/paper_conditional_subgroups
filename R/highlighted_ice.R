#' Highlighted ICE curves
#'
#' @param pred iml::Predictor
#' @param feature Name of the feature
#' @return list of plot and results
highlighted_ice = function(pred, feature, cond){
  eff = FeatureEffect$new(pred, feature = feature, method = "pdp")
  xgrid = eff$results[[feature]]
  dens = cond$cdens(pred$data$X, xgrid = xgrid)
  ice = FeatureEffect$new(pred, feature = feature, method = "ice")$results
  ice = merge(ice, dens, by.x = c(".id", fname), by.y = c(".id.dist", fname)) 
  p = ggplot(ice) + geom_line(aes_string(x = feature, group = ".id", alpha = ".dens", y = ".value"))
  list(results = ice, dens = dens,  p = p)
}
