# For a grid (grid.dat) of features (param features) creates a blown up
# dataset with the marginals of features not in 'features'.
# The samples (n.sample.dist number of samples) for the marginals are drawn from dist.dat.
#            If n.sample.dist is not set, the whole cartesian product between grid.dat and dist.dat is built
# grid.dat only needs to contain the columns which are fixed. Decide here which grid points should be used.
# dist.dat needs to contain all columns
Conditional = R6Class(
  public = list(
    feature = NULL,
    data = NULL,
    model = NULL,
    ctrl = NULL,
    initialize = function(data, feature, ctrl = partykit::ctree_control()) {
      assert_class(data, "Data")
      self$data = data
      self$feature = feature
      self$ctrl = ctrl
      private$fit_conditional()
    },
    csample_data = function(X, size){
      cmodel = self$models[[self$feature]]
      data_nodes = self$cnode(self$data$X)
      X_nodes = self$cnode(X)
      xj_samples = lapply(1:nrow(X), function(i) {
        node = X_nodes[i, "node"]
        data_ids = which(data_nodes$node == node)
        data_ids = setdiff(data_ids, i)
        data_ids_sample = sample(data_ids, size = size, replace = TRUE)
        xj = self$data$X[data_ids_sample, self$feature, with = FALSE]
        data.frame(t(xj))
      })
      rbindlist(xj_samples)

    },
    csample_parametric = function(X, size){
      cmodel = self$models[[self$feature]]
      if (self$data$feature.types[[self$feature]] == "categorical") {
        xgrid = unique(self$data$X[[self$feature]])
      } else {
        x = self$data$X[[self$feature]]
        xgrid = seq(from = min(x), to = max(x), length.out = 100)
      }
      dens = self$cdens(X, xgrid)
      xj_samples = lapply(1:nrow(X), function(i) {
        dens_i = dens[dens$.id.dist == i,]
        xj = sample(dens_i[[self$feature]], size = size, prob = dens_i[[".dens"]], replace = TRUE)
        data.frame(t(xj))
      })
      rbindlist(xj_samples)
    },
    csample = function(X, size, type = "parametric"){
      assert_number(size, lower = 1)
      assert_character(self$feature)
      assert_data_table(X)
      assert_choice(type, c("data", "parametric"))
      if (type == 'parametric') {
        self$csample_parametric(X, size)
      } else {
        self$csample_data(X, size)
      }
      },
    cdens = function(X, xgrid = NULL){
      cmodel = self$model
      if(inherits(cmodel, "trafotree")) {
        conditionals = predict(cmodel, newdata = X, type = "density", q = xgrid)
        densities = melt(conditionals)$value
        densities = data.table(.dens = densities, .id.dist = rep(1:nrow(X), each = length(xgrid)),
                   feature = rep(xgrid, times = nrow(X)))
      } else if (self$data$feature.types[self$feature] == "categorical") {
        probs = predict(cmodel, newdata = X, type = "prob")
        probs.m = melt(probs)$value
        densities = data.table(.dens = probs.m, .id.dist = rep(1:nrow(X), each = ncol(probs)),
                   feature = factor(rep(colnames(probs), times = nrow(X)), levels = levels(self$data[[self$feature]])))
      } else {
        pr = predict(cmodel, newdata = X, type = "density")
        at = unique(X[[self$feature]])
        res = sapply(pr, function(pr) pr(at) / sum(pr(at)))
        res = data.table(t(res))
        colnames(res) = as.character(at)
        res.m = melt(res, measure.vars = as.character(at))
        densities = data.table(.dens = res.m$value, .id.dist = rep(1:nrow(X), times = ncol(X)),
                               feature = rep(at, each = nrow(X)))
      }
      colnames(densities) = c(".dens", ".id.dist", self$feature)
      densities
  },
  cnode = function(X,  prob = NA) {
    cmodel = self$model
    node = predict(cmodel, newdata = X, type = "node")
    node_df = data.frame(node = (node), .id = names(node), .path = pathpred(cmodel, X))
    if (!is.na(prob)){
      if(inherits(cmodel, "trafotree")) {
        # case of numerical feature
        quants = predict(cmodel, newdata = X, type = "quantile", prob = prob)
        quants = data.frame(t(quants))
        colnames(quants) = paste0("q", prob)
      } else if (self$data$feature.types[[self$feature]] == "numerical") {
        # case of numerical features with few unique values
        quants = predict(cmodel, newdata = X, type = "quantile", at = prob)
        colnames(quants) = paste0("q", prob)
      } else {
        # case of categorical feature
        quants = predict(cmodel, newdata = X, type = "prob")
        names(quants) = levels(X[[self$feature]])
      }
      node_df = cbind(node_df, quants)
    }
    node_df
  }
  ),
  private = list(
  fit_conditional = function() {
    require("trtf")
    y = self$data$X[[self$feature]]
    if ((self$data$feature.types[self$feature] == "numerical") & (length(unique(y)) > 2)) {
      yvar = numeric_var(self$feature, support = c(min(y), max(y)))
      By  =  Bernstein_basis(yvar, order = 5, ui = "incr")
      m = ctm(response = By,  todistr = "Normal", data = self$data$X )
      form = as.formula(sprintf("%s ~ 1 | .", self$feature))
      part_cmod = trafotree(m, formula = form,  data = self$data$X, control = self$ctrl)
    } else {
      form = as.formula(sprintf("%s ~ .", self$feature))
      part_cmod = partykit::ctree(form, data = self$data$X, control = self$ctrl)
    }
      self$model = part_cmod
    }
  )
)



#' Fit conditional models
#'
#' Needed for conditional PDP and Feature Importance.
#'
#' @param data data.frame with data for which to fit the conditional models
#' @param ctrl ctree_control object
#' @param type Either trtf or cart. Decide which type of tree to use
#' @param features vector of feature names. fits conditionals for all feature if NULL
#' @return list of Conditional R6 objects
#' @importFrom partykit ctree_control
#' @export
fit_conditionals = function(data,  ctrl = NULL,  type = "trtf", features = NULL){
  assert_choice(type, c("trtf", "cart"))
  assert_data_frame(data)
  if (is.null(features)) features = colnames(data)
  cmods = lapply(features, function(fname){
    if (type == "trtf") {
      if (is.null(ctrl)) ctrl = ctree_control(minbucket = 30)
      Conditional$new(iml:::Data$new(data.frame(data)), fname, ctrl = ctrl)
    } else {
      if (is.null(ctrl)) ctrl = rpart::rpart.control(minbucket = 30)
      ConditionalCART$new(iml:::Data$new(data.frame(data)), fname, ctrl = ctrl)
    }
  })
  names(cmods) = features
  cmods
}


ConditionalCART = R6Class(
  inherit = Conditional,
  public = list(
    feature = NULL,
    data = NULL,
    model = NULL,
    ctrl = NULL,
    initialize = function(data, feature, ctrl = rpart::rpart.control(minbucket = 30)) {
      assert_class(data, "Data")
      self$data = data
      self$feature = feature
      self$ctrl = ctrl
      private$fit_conditional()
    }
    ),
  private = list(
  fit_conditional = function() {
    require("rpart")
    form = as.formula(sprintf("%s ~ .", self$feature))
    cart = rpart(form, data = self$data$X, control = self$ctrl)
    self$model = partykit::as.party(cart)
    }
  )
)

ConditionalStrobl= R6Class(
  inherit = Conditional,
  public = list(
    feature = NULL,
    data = NULL,
    model = NULL,
    ctrl = NULL,
    initialize = function(data, feature, ctrl = partykit::ctree_control()) {
      assert_class(data, "Data")
      self$data = data
      self$feature = feature
      self$ctrl = ctrl
      private$fit_conditional()
    },
		cnode = function(X) {
			cmodel = self$model
			node = predict(cmodel, newdata = X, type = "node")[[1]]
			data.frame(node = node, .id = node, .path = NA)
		}),
  private = list(
  fit_conditional = function() {
    require("partykit")
		form = as.formula(sprintf("%s ~ .", self$feature))
		self$model = partykit::cforest(form, data = self$data$X, control = self$ctrl, ntree = 1)
  })
)


