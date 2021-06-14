library("data.table")
library("MASS")
library("here")
library("iml")
library("parallel")
library("batchtools")

devtools::load_all("../")
data_dir = sprintf("%s/data/", here())
set.seed(1)

N_EXPERIMENTS = 100
NC = 7

#' compute true PFI given f and g
#'
#' @param X data.frame with feature
#' @param fun function f(X) = y
#' @param g function g(X_{-1}) = X1
#' @return conditional feature importance
tpfi = function(X, fun, g){
  n = nrow(X)
  y = fun(X) + rnorm(n, 0, 1)
  L = (y - fun(X))^2
  Xt = X
  Xt$X1 = g(X)
  Lt = (y - fun(Xt))^2
  mean(Lt - L)
}
# =============================================================================
# Setup registry for batchtools
# =============================================================================
if(file.exists("registry")) {
  reg = loadRegistry("registry", writeable = TRUE)
} else {
  reg = makeExperimentRegistry(file.dir = "registry")
}
#clearRegistry(reg)


# =============================================================================
# Setup Experiment 1 with true g and f
# =============================================================================

ex1 = function(data, job, n, p, ...){
  ntotal = 1000 + n 
  g = gs[[job$prob.pars$g]]
  f = fs[[job$prob.pars$f]]
  X = data.frame(matrix(rnorm(n = p * ntotal), ncol = p))
  # Overwrite X1 as function of other features
  X$X1 = g(X)
  pred.fun = function(mod, newdata) f(newdata)
  y = f(X) + rnorm(n = ntotal, mean = 0, sd = 1)
  tpfi = tpfi(X[1:1000,], f, g)
  train_index = sample(1:n, size = (2/3) * n, replace = FALSE)
  test_index = setdiff(1:n, train_index)
  ytrain = y[train_index]
  ytest = y[test_index]
  Xtrain = X[train_index, ]
  Xtest = X[test_index,]
  pred = Predictor$new(predict.fun = pred.fun, data = Xtest, y = ytest) 
  # todo train and test
  list("ytrain" = ytrain, "Xtrain" = Xtrain,
       "ytest" = ytest, "Xtest" = Xtest,
       pred = pred, tpfi = tpfi)
}


# =============================================================================
# Setup Experiment 2 with intermediate cforest
# =============================================================================
ex2 = function(data, job, n, p, ...){
  # Use stuff from experiment before
  d = ex1(data, job, n, p, ...)
  # Update the predictor to be a conditional random forest
  dat = cbind(y = d$ytrain, d$Xtrain)
  d$dat = dat
  cf = randomForest::randomForest(y ~ ., data = dat, 
                                  keep.forest = TRUE, keep.inbag = TRUE)
  pred_fun = function(mod, newdata) {predict(mod, newdata = newdata)}
  d$pred = Predictor$new(mod = cf, predict.fun = pred_fun, data = d$Xtest, y = d$ytest)
  # TODO: Use TPFI here to generate groundtruth and compare against this for setting 2
  fun = function(x) predict(cf, newdata =  x)
  X = data.frame(matrix(rnorm(n = p * 1000), ncol = p))
  # Overwrite X1 as function of other features
  g = gs[[job$prob.pars$g]]
  X$X1 = g(X)
  # Compute true conditional PFI
  n = nrow(X)
  f = fs[[job$prob.pars$f]]
  y = f(X) + rnorm(n = nrow(X), mean = 0, sd = 1)
  L = (y - fun(X))^2
  Xt = X
  Xt$X1 = g(X)
  Lt = (y - fun(Xt))^2
  d$tpfi_rf = mean(Lt - L)
  d$cf = cf
  d
}


addProblem(name = "e1", data = data.frame(), fun = ex1, seed = 1)
addProblem(name = "e2", data = data.frame(), fun = ex2, seed = 2)
# =============================================================================
# Setup cs-PFI (subgroups)
# =============================================================================
#' Subgroup PFI
#'
#' @param data batchtools data
#' @param job batchtools job
#' @param instance batchtool instance
#' @param return_nterminal Should the number of terminal nodes be returned? Means that data.frame is returned instead of scalar
#' @param alpha The alpha value for ctree_control
#' @param cp The complexity parameter for rpart.control
imp_tree = function(data, job, instance, return_nterminal = FALSE, ...){
 dat = iml:::Data$new(cbind(instance$Xtrain), y = instance$ytrain)
 maxdepth = job$algo.pars$maxdepth
 minbucket = job$algo.pars$minbucket

 if (job$algo.pars$type == "cart") {
   ctrl1 = rpart::rpart.control(maxdepth = maxdepth, minbucket = minbucket,
                                cp = job$algo.pars$cp)
   cond = ConditionalCART$new(dat, "X1", ctrl = ctrl1)
 } else {
   ctrl1 = partykit::ctree_control(maxdepth = maxdepth, minbucket = minbucket,
                                   alpha = job$algo.pars$alpha) 
   cond = Conditional$new(dat, "X1", ctrl = ctrl1)
 }
 group = cond$cnode(instance$Xtest)$.path
 gpfi =  grouped_pfi_feat(instance$pred, group = group, loss = mse, fname = "X1")

 res = weighted.mean(gpfi$importance, w = gpfi$n)
 if (return_nterminal){
   nterminal = length(nodeids(cond$model, terminal = TRUE))
   res = data.frame(pfi = res, nterminal = nterminal)
 } 
 res
}

# =============================================================================
# Setup true PFI
# =============================================================================
imp_true = function(data, job, instance, ...){
  instance$tpfi
}

# True PFI, but for random forest from setting 2
imp_true_rf = function(data, job, instance, ...){
  instance$tpfi_rf
}


# =============================================================================
# Setup knockoff PFI
# =============================================================================
imp_ko = function(data, job, instance, ...){
  knockoff_imp2(Xtrain = instance$Xtrain, ytrain = instance$ytrain, pred = instance$pred, features = "X1")$importance
}

# =============================================================================
# Setup imputation PFI
# =============================================================================
imp_imp = function(data, job, instance, ...){
  # Imputation Approach with Random Forest
  RF = ranger::ranger(X1 ~ ., data = instance$Xtrain)
  # We use approach from fisher here with residuals
  resids = instance$Xtest$X1 - predict(RF, instance$Xtest)$predictions
  X2c = instance$Xtest
  X2c$X1 = predict(RF, X2c)$predictions + sample(resids, size = nrow(X2c), replace = FALSE)
  yhat = instance$pred$predict(X2c)[[1]]
	mse_orig = mean((instance$pred$predict(instance$Xtest)[[1]] - instance$ytest)^2)
  mean((yhat - instance$ytest)^2) - mse_orig
}


# =============================================================================
# Setup marginal PFI
# =============================================================================
imp_marg = function(data, job, instance, ...){
	mse_orig = mean((instance$pred$predict(instance$Xtest)[[1]] - instance$ytest)^2)
  X2 = instance$Xtest
  X2$X1 = sample(X2$X1)
  yhat = instance$pred$predict(X2)[[1]]
  mean((yhat - instance$ytest)^2) - mse_orig
}

# =============================================================================
# Setup CVIRF (Stroble)
# =============================================================================
imp_cvirf = function(data, job, instance, ...){
  # the training data for the random forest must be available in
  # the same environment as the permimp call, that's how permimp
  # is implemented.
  if ("cf" %in% names(instance)) {
    dat = instance$dat
    th = job$algo.pars$permimp_threshold
    cf = instance$cf
    permimp::permimp(cf, conditional = TRUE, do_check = FALSE, threshold = th,
                         progressBar = FALSE)$values[["X1"]]
  } else {
    NA
  }
}

# see file R/ci-experiment.R
addAlgorithm(name = "tree",  imp_tree)
addAlgorithm(name = "true",  imp_true)
addAlgorithm(name = "true_rf",  imp_true_rf)
addAlgorithm(name = "imp",  imp_imp)
addAlgorithm(name = "ko",  imp_ko)
addAlgorithm(name = "marg", imp_marg)
addAlgorithm(name = "cvirf", imp_cvirf)

# =============================================================================
# Define scenarios for g(X_{-1}) = X_1
# =============================================================================
g0 = function(X) {rnorm(n = nrow(X), mean = 0, sd = 1)}
g1 = function(X) {rnorm(n = nrow(X), mean = X$X2, sd = 1)}
g2 = function(X) {
  mu = rep(3, times = nrow(X))
  # automatically: mu[X$X2 > 0] = 3
  mu[(X$X2 <= 0) & (X$X3 > 0)] = -3
  mu[(X$X2 <= 0) & (X$X3 <= 0)] = 0

  sigm = rep(1, times = nrow(X))
  # automatically: sigm[X$X2 > 0] = 1
  sigm[(X$X2 <= 0) & (X$X3 > 0)] = 2
  sigm[(X$X2 <= 0) & (X$X3 <= 0)] = 5

  rnorm(nrow(X), mu, sigm)
}

g3 = function(X) {
  mu = rowSums(X[sprintf("X%i", 2:10)])
  rnorm(nrow(X), mu, 5)
}

gs = list("independent" = g0,
          "linear_g" = g1,
          "nlinear_g" = g2, 
          "high_dim_g" = g3)
# =============================================================================
# Define y = f(X)
# =============================================================================
fn1 = function(X) {rowSums(X[,1:10]) + X$X1 * X$X2}
fs = list("linear_f" = fn1)



# =============================================================================
# Define scenarios and run simulation
# =============================================================================
s1 = expand.grid(p = c(10,  90),
								 n = c(300, 3000),
								 g = names(gs),
								 f = names(fs))

tree_ades = data.frame(maxdepth = 30, minbucket = 30, type = c("cart", "trtf"),
                       alpha = 0.05, cp = 0.01)

ades = list("tree" = tree_ades,
            "true" = data.frame(),
            "imp" = data.frame(),
            "ko" = data.frame(),
            "marg" = data.frame())

# Experiments for setting 1 (true f)
#addExperiments(list(e1 = s1), ades, repls = N_EXPERIMENTS)

# Experiments for setting 2 (intermediate cforest)
ades_rf = list("true_rf" = data.frame(), 
               "cvirf" = data.frame(permimp_threshold = c(0.95)))
ades2 =  append(ades, ades_rf)
# DELME
ades2 = list("true_rf" = data.frame())
#addExperiments(list(e2 = s1[s1$n == 1000,]), ades2, repls = N_EXPERIMENTS)
addExperiments(list(e2 = s1), ades2, repls = N_EXPERIMENTS)
# Fewer runs when n=3000, for now, as its very time consuming (intermediate cforest)
#addExperiments(list(e1 = s1, e2 = s1[s1$n == 3000,]), ades2, repls = 100)
summarizeExperiments()
reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = NC, fs.latency = 0)
#testJob(70)
submitJobs()
waitForJobs()
getErrorMessages()
# =============================================================================
# Save results
# =============================================================================
pars = unwrap(getJobPars())
pfis = reduceResultsDataTable()
res = merge(pars, pfis)
saveRDS(res, file = sprintf("%s/results/true-experiment.RDS", here()))

if (FALSE) {
# ============================================================================= 
# Tree Depth Experiments
# =============================================================================
if(file.exists("registry-nterminal")) {
  reg = loadRegistry("registry-nterminal", writeable = TRUE)
} else {
  reg = makeExperimentRegistry(file.dir = "registry-nterminal")
}
tree_ades = expand.grid(maxdepth = c(1:5),
												minbucket = c(30),
												type = c("cart", "trtf"),
                        return_nterminal = TRUE,
                        alpha = 1,
                        cp = 0)
ades = list("tree" = tree_ades)
addAlgorithm(name = "tree",  imp_tree)
addProblem(name = "e1", data = data.frame(), fun = ex1, seed = 1)
addExperiments(list(e1 = s1), ades, repls = N_EXPERIMENTS, reg = reg)

summarizeExperiments()
reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = NC, fs.latency = 0)
#testJob(1)
submitJobs()
waitForJobs()
# =============================================================================
# Save results
# =============================================================================
pars = unwrap(getJobPars())
pfis = reduceResults(function(a,b){ rbind(a, b)})
pfis$job.id = 1:nrow(pfis)
res = merge(pars, pfis)
saveRDS(res, file = sprintf("%s/results/true-experiment-depth.RDS", here()))
}
