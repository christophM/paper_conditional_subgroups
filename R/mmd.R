#' Radial basis function kernel
#'
#' Computes kernel matrix for one or two datasets
#' Mainly wraps kernlab::kernelMatrix
#'
#' @param d1 data.frame, one instance per row, one feature per column
#' @param d2 data.frame, same columns as d1
#' @param sigma The sigma of the kernel computation
#' @return Mean of similarity matrix
cross.kernel = function(d1, d2 = NULL, sigma) {
  checkmate::assert_data_frame(d1, any.missing = FALSE)
  checkmate::assert_data_frame(d2, null.ok = TRUE, ncols = ncol(d1), any.missing = FALSE)
  checkmate::assert_number(sigma, lower = 0)
  mradial = kernlab::rbfdot(sigma = sigma)
  mm = kernlab::kernelMatrix(mradial, as.matrix(d1), as.matrix(d2))
  mean(mm)
}

#' Maximum mean discrepancy (squared) 
#'
#' @param d1 data.frame, one instance per row, one feature per column
#' @param d2 data.frame, same columns as d1
#' @param sigma The sigma of the kernel computation
#' @return Maximum mean discrepancy
mmd2 = function(d1, d2, sigma = NULL) {
  checkmate::assert_data_frame(d1, any.missing = FALSE)
  checkmate::assert_data_frame(d2, null.ok = TRUE, ncols = ncol(d1), any.missing = FALSE)
  checkmate::assert_number(sigma, lower = 0, null.ok = TRUE)
  d1 = data.frame(model.matrix(~ . -1, data = d1))
  d2 = data.frame(model.matrix(~ . -1, data = d2))
  if(is.null(sigma)) {
    sigma = get_median_dist(rbind(d1, d2))
    # Confusingly, the sigma in rbfdot is acutally the gamma param
    sigma = 1/(2 * sigma^2)
  }
  cross.kernel(d1, d1, sigma = sigma) - 2 * cross.kernel(d1, d2, sigma = sigma) + cross.kernel(d2, d2, sigma = sigma)
}

#' MMD2 for PDP intervention
#'
#' If xgrid is NULL, xj will be drawn from marginal distribution, which implies that each xj value is a grid value.
#' If xgrid is set, only xgrid values will be sampled.
#'
#' @param dat data.frame
#' @param fname character with name of column
#' @param xgrid (optional) vector with xj values
#' @return MMD2 between full data and intervened data
	
mmd2_pdp = function(dat, dat_ref, fname, xgrid = NULL, grid_type = "quantile") {
  checkmate::assert_data_frame(dat_ref)
  checkmate::assert_data_frame(dat)
  checkmate::assert_character(fname)
  checkmate::assert_true(fname %in% colnames(dat))
  checkmate::assert_numeric(xgrid, null.ok = TRUE)
	checkmate::assert_choice(grid_type, c("quantile", "equidistant"))
  dat_pdp = dat
  if(is.null(xgrid)) {
		if(grid_type == "quantile") {
			dat_pdp[fname] = sample(dat_pdp[[fname]])
		} else {
			dat_pdp[fname] = runif(1, min(dat_pdp[[fname]]), max(dat_pdp[[fname]]))
		}
  } else {
    dat_pdp[fname] = sample(xgrid, size = nrow(dat), replace = TRUE)
  }
  mmd2(dat_ref, dat_pdp)
}

#' MMD2 for cgPDP intervention
#'
#' If xgrid is NULL, xj will be drawn from marginal distribution per group, which implies that each xj value is a grid value.
#' If xgrid is set, the sample values will be rounded to the closes xgrid value.
#'
#' @param dat data.frame
#' @param fname character with name of column
#' @param xgrid (optional) vector with xj values
#' @return MMD2 between full data and intervened data
mmd2_cond = function(dat, dat_ref, cond, fname, xgrid = NULL, grid_type = "quantile") {
  checkmate::assert_data_frame(dat)
  checkmate::assert_data_frame(dat_ref)
  checkmate::assert_character(fname)
  checkmate::assert_true(fname %in% colnames(dat))
  checkmate::assert_numeric(xgrid, null.ok = TRUE)
  checkmate::assert_class(cond, "R6")
	checkmate::assert_choice(grid_type, c("quantile", "equidistant"))
  if(is.null(xgrid)) {
		nodes_df = cond$cnode(X = data.table(dat))
		nodes_df$feature = dat[[fname]]
		if(grid_type == "quantile") {
			xnew = nodes_df %>%
				group_by(node) %>%
				mutate(xnew = sample(feature)) %>%
				pull(xnew)
		} else {
			xnew = nodes_df %>%
				group_by(node) %>%
				mutate(xnew = runif(1, min(feature), max(feature))) %>%
				pull(xnew)
		}
  } else  {
    nodes_df = cond$cnode(data.table(dat))
    nodes_df$feature = dat[[fname]]
    nodes_s = nodes_df %>%
      group_by(node, .path) %>%
      dplyr::summarize(xmin= min(feature), xmax = max(feature), n = n())
    xnew = rep(NA, times = nrow(dat))
    for (node in nodes_s$node) {
      node_info = nodes_s[nodes_s$node == node,]
      xgrid_node = xgrid[(xgrid >= node_info$xmin) & (xgrid <= node_info$xmax)]
      if (length(xgrid_node) == 0) {
        # This is the case without intersection of grid point.
        # Node not drawn in this case and therefore NA
        grid_sample = NA
      } else {
        grid_sample = sample(nodes_df$feature[nodes_df$node == node])
      }
      xnew[nodes_df$node == node] = grid_sample
    }
  }
  dat_cond = dat
  dat_cond[fname] = xnew
  n1 = nrow(dat_cond)
  dat_cond = na.omit(dat_cond)
  if(nrow(dat_cond) < n1) print(sprintf("%i/%i data points omitted because grid too sparse", n1 - nrow(dat_cond), n1))
  mmd2(dat_ref, dat_cond)
}

#' MMD2 for ALE intervention
#'
#' Computed as the average MMD2 by shifting to the left interval and to the right interval
#'
#' @param dat data.frame
#' @param fname character with name of column
#' @param xgrid (optional) vector with xj values
#' @return MMD2 between full data and intervened data

mmd2_ale = function(dat, dat_ref, fname, xgrid = NULL) {
  checkmate::assert_data_frame(dat)
  checkmate::assert_data_frame(dat_ref)
  checkmate::assert_character(fname)
  checkmate::assert_true(fname %in% colnames(dat))
  checkmate::assert_numeric(xgrid, null.ok = TRUE)
  if (is.null(xgrid)) return(NA)
  interval.index = findInterval(dat[[fname]], xgrid, all.inside = TRUE)
  # Data point in the left most interval should be in interval 1, not zero
  ddat1 = ddat2 = dat
  ddat1[fname] = xgrid[interval.index] 
  ddat2[fname] = xgrid[interval.index + 1] 
  mmd2(dat_ref, rbind(ddat1, ddat2))
}


#' Compute the median distance
#'
#' To be used for setting a default sigma param in MMD computation
#' Based on https://papers.nips.cc/paper/3110-a-kernel-method-for-the-two-sample-problem.pdf
#'
#' @param d data.frame with data points
#' @return median distance 
get_median_dist = function(d){
  checkmate::assert_data_frame(d)
  dists = dist(d, diag = FALSE, upper = FALSE, method = "euclidean")
  median(dists)
}
