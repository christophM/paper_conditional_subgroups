#' Create grid along feature
#'
#' @param feature vector with feature values
#' @param grid.size numeric size of grid
#' @param feature.type either "numerical" or "categorical"
#' @param type either "equidistant" or "quantile" 
#' @return grid values
get_grid = function(feature, grid.size,  feature.type = NULL, type = "equidistant") {
  checkmate::assert_vector(feature, all.missing = FALSE, min.len = 2)
  checkmate::assert_choice(feature.type, c("numerical", "categorical"), null.ok = TRUE)
  checkmate::assert_numeric(grid.size)
  checkmate::assert_choice(type, c("equidistant", "quantile"))
  
  if(is.null(feature.type)) feature.type = get.feature.type(class(feature))
  
  if (feature.type == "numerical") {
    # remove NaN NA and inf
    feature = feature[is.finite(feature)]
    if (length(feature) == 0) stop("Feature does not contain any finite values.")
    
    if(type == "equidistant") {
      grid = seq(from = min(feature), 
        to = max(feature), 
        length.out = grid.size)
    } else if(type == "quantile") {
      probs = seq(from = 0, to = 1, length.out = grid.size)
      grid = quantile(feature, probs = probs, names = FALSE, type = 1)
    }
  } else if (feature.type == "categorical") {
    grid = unique(feature)
  }
  grid
}

