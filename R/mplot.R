#' M Plot
#'
#' Compute the M Plot as described in the ALE paper
#' 
#' @param pred iml::Predictor
#' @param feature character of feature name
#' @param xgrid numerical vector with the grid values
#' @param window_width width of window around each grid point for neighbourhood
mplot = function(pred, feature, xgrid, window_width) {
  assert_numeric(xgrid)
  assert_character(feature)
  assert_numeric(window_width)
  dat = pred$data$X
  feat = dat[[feature]]
  # for each grid point
  res  = lapply(xgrid, function(x) {
    datx = dat[((x - window_width/2) <= feat) &
               ((x + window_width/2) >= feat), , drop = FALSE]
    if (nrow(datx) == 0) {
      data.frame(.value = NA,
                 .x = x)
    } else {
      predictions = pred$predict(datx)[[1]]
      data.frame(.value = mean(predictions),
                 .x = x)
    }
  })
  res = rbindlist(res)
  colnames(res) = c(".value", feature)
  res$.type = "m-plot"
  p = ggplot(res) + geom_line(aes_string(x = feature, y = ".value"))
  list(results = res, p = p)
}
