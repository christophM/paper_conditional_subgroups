knockoff_imp = function(X, y, pred, type, features = NULL) {
  assert_choice(type, c("fixed", "gaussian"))
	if (type == "fixed") {
		ko = knockoff::create.fixed(X)$Xk
	} else {
		ko = create.second_order(as.matrix(X))
	}
	mse_orig = mean((pred$predict(X)[[1]] - y)^2)
  if (is.null(features)) features = colnames(X)
	mses = lapply(features, function(fname) {
	 if(type == "fixed") ko[,fname] = ko[,fname] * sqrt(sum(ko[,fname]^2)) 
	 Xnew = X
	 Xnew[[fname]] = ko[,fname]
	 yhat = pred$predict(Xnew)[[1]]
	 # MSE
	 ms = mean((yhat - y)^2)
  })
	mses = unlist(mses)
	data.frame(feature = features, importance = mses - mse_orig)
}


knockoff_imp2 = function(Xtrain, ytrain, pred, features = NULL) {
  Xtest = pred$data$X
  ytest = pred$data$y[[1]]
  # Code mostly taken from knockoff:::create.second_order
  # Not using create.second_order directly, as I want to split training and testing
  mu = colMeans(Xtrain)
  Sigma = cov(Xtrain)
  if (!knockoff:::is_posdef(Sigma)) stop("not positive definite")
  ko = create.gaussian(as.matrix(Xtest), mu, Sigma, method = "asdp")
	mse_orig = mean((pred$predict(Xtest)[[1]] - ytest)^2)
  if (is.null(features)) features = colnames(Xtest)
	mses = lapply(features, function(fname) {
	 Xnew = Xtest
	 Xnew[[fname]] = ko[,fname]
	 yhat = pred$predict(Xnew)[[1]]
	 # MSE
	 mean((yhat - ytest)^2)
  })
	mses = unlist(mses)
	data.frame(feature = features, importance = mses - mse_orig)
}

