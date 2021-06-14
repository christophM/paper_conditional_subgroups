# =============================================================================
# Data fidelity (MMD) experiment
# =============================================================================

library("foreach")
library("dplyr")
library("ggplot2")
library("tidyr")
library("here")
library("parallel")
library("rjson")
library("mlr3")
library("mlr3pipelines")


MAXROWS = 10000
N_SIMULATIONS = 30
N_FEATURES = 10
NCORES = 31

# Loads custom functions from R/ folder
devtools::load_all()

clean_names = function(x) {
	x = gsub("[^0-9a-zA-Z]+", "", x)
	# R adds a "." suffix when in, for, if are used as colnames in data.frame 
	x = gsub("for", "FOR", x, fixed = TRUE)
  x = gsub("if", "IF", x, fixed = TRUE)
  gsub("in", "IN", x, fixed = TRUE)
}


# Load data sets
data_dir = here("data")
# In processed are all datasets with right amount of features and shuffled feature order
openml_dir = here("data/cc18/processed/")
# In raw are all datasets + the meta-information
openml_dir_raw = here("data/cc18/raw/")

set.seed(10)
test_size = 0.6

data_ids = unique(gsub("\\.(arff|json)", "", list.files(openml_dir)))

set.seed(42)
experiments = expand.grid(data_id = rev(data_ids),
                          sim = 1:N_SIMULATIONS,
                          grid_type = "quantile",
                          stringsAsFactors = FALSE)
# Shuffles order of experiments
experiments = experiments[sample(1:nrow(experiments)),]
print(sprintf("Doing %i experiments", nrow(experiments)))
results = mclapply(1:nrow(experiments),
                   mc.cores = NCORES,
                   mc.preschedule = FALSE,
                   function(i) {

  ex = experiments[i,]
  print(ex)
  data_id = as.numeric(as.character(experiments[i, "data_id"]))
  j = fromJSON(file = sprintf("%s/%s.json", openml_dir_raw, data_id))
  data_name = j$name


  filename = paste0(res_dir, "/data-fidelity/", i, "-", data_name,  ".csv")
  if(file.exists(filename)){
    print(sprintf("Skipping %s, already computed", filename))
    return(TRUE)
  }

  # Load data
  dat = farff::readARFF(sprintf("%s/%s.arff", openml_dir, data_id), show.info = FALSE)
  print(sprintf("Processing data %i (n=%i, p=%i)", data_id, nrow(dat), ncol(dat)))
  # Some column names caused troubles
  colnames(dat) = clean_names(colnames(dat))
  # Also replace special characters in target
  # We need target, bc. of Strobl approach and to remove for others
  target =  clean_names(j$default_target_attribute)
  if(is.logical(dat[[target]])) dat[,target] = as.numeric(dat[,target])
  # Remove rows with missing values
  dat = na.omit(dat)
  # Downsample to keep it manageable
  if(nrow(dat) > MAXROWS) dat = dat[sample(1:nrow(dat), size = MAXROWS), ]

  train_ids = sample(1:nrow(dat), size = (1 - test_size) * nrow(dat))
  train_dat = dat[train_ids,]
  test_ids = setdiff(seq_len(nrow(dat)), train_ids)
  # For creating a task, bc. often no default target in openmml. Bit ugly, I know.
  dat$fake_target = rep(1, times = nrow(dat))
  test_dat = dat[test_ids,]
  task = TaskRegr$new(id = "dat", backend = dat, target = "fake_target")
  pos = po("scale")
  train_dat= task$clone()$filter(train_ids)
  test_dat = task$clone()$filter(test_ids)
  pos$train(list(train_dat))[[1]]
  train_dat = pos$predict(list(train_dat))[[1]]$data()
  train_dat = data.frame(train_dat)
  test_dat = pos$predict(list(test_dat))[[1]]$data()
  test_dat = data.frame(test_dat)
	# For some reason, this operations adds a point to keywords such as 'if' and 'in'
  train_dat$fake_target = test_dat$fake_target = NULL

  fnames = colnames(dat)[!(colnames(dat) %in% c("fake_target", target))] 
  train_dat = train_dat[fnames]
  test_dat = test_dat[fnames]
  # Further split test_dat for MMD computation
  test_dat_ref_index = sample(c(TRUE, FALSE), size = nrow(test_dat), replace = TRUE)
  test_dat_ref = test_dat[test_dat_ref_index,]
  test_dat = test_dat[!test_dat_ref_index,]

  # Cols are shuffled in load-cc18.R, so using first 10 is random sample
  fnames_sample = fnames[1:min(10, length(fnames))]
  res = lapply(fnames_sample, function(fname) {
    print(fname)
    idat = iml:::Data$new(train_dat)
		compute_mmd_tree = function(maxdepth, treetype){
			assert_choice(treetype, c("cart", "trtf"))
			if(treetype == "cart"){
				ctrl1 = rpart::rpart.control(maxdepth = maxdepth)
				cond = ConditionalCART$new(idat, fname, ctrl = ctrl1)
			} else {
				ctrl1 = ctrl(maxdepth)
				cond = Conditional$new(idat, fname, ctrl = ctrl1)
			}
			mmdc= mmd2_cond(test_dat, dat_ref = test_dat_ref, fname = fname, cond = cond, grid_type = ex$grid_type)
			rm(cond)
			mmdc
		}

		# Building up the results data.frame
    mmds = data.frame(repetition = ex$sim,
											feature = fname,
											grid_type = ex$grid_type,
											data_name = data_name)

    mmds$trtf1 = compute_mmd_tree(1, "trtf")
    mmds$trtf2 = compute_mmd_tree(2, "trtf")
    mmds$trtf3 = compute_mmd_tree(3, "trtf")
    mmds$trtf4 = compute_mmd_tree(4, "trtf")
    mmds$trtf5 = compute_mmd_tree(5, "trtf")
    mmds$trtf30 = compute_mmd_tree(30, "trtf")

    mmds$cart1 = compute_mmd_tree(1, "cart")
    mmds$cart2 = compute_mmd_tree(2, "cart")
    mmds$cart3 = compute_mmd_tree(3, "cart")
    mmds$cart4 = compute_mmd_tree(4, "cart")
    mmds$cart5 = compute_mmd_tree(5, "cart")
    mmds$cart30 = compute_mmd_tree(30, "cart")


    # Strobl a-like approach
    train_dat2 = train_dat
    train_dat2[target] = dat[train_ids, target]
    idat = iml:::Data$new(data.frame(train_dat2))
    cond = ConditionalStrobl$new(idat, target, ctrl(30))
    mmds$strobl = mmd2_cond(test_dat, dat_ref = test_dat_ref, cond = cond, fname = fname)

    # Missing: Only use for splitting all features that have high dependence with feature of interest

    # Imputation Approach with Random Forest
    form = as.formula(sprintf("%s ~ .", fname))
    RF = ranger::ranger(form, data = train_dat)
		# We use approach from fisher here with residuals
    resids = train_dat[fname] - predict(RF, train_dat)$predictions
    test_dat_copy = test_dat
		test_dat_copy[fname] = predict(RF, test_dat)$predictions + sample(resids[[1]], size = nrow(test_dat))
		mmds$imputation = mmd2(test_dat_copy, test_dat_ref)

    # Compute mmd2 for marg
    mmds$perm = mmd2_pdp(test_dat, dat_ref = test_dat_ref, fname = fname, grid_type = ex$grid_type)

    # For ale
    xgrid = get_grid(train_dat[[fname]], grid.size = 30,
                     type = ex$grid_type, feature.type = "numerical")
    xgrid = unique(xgrid)

    mmds$ale = mmd2_ale(test_dat, dat_ref = test_dat_ref, fname, xgrid = xgrid)
    mmds$none = mmd2(test_dat, test_dat_ref)

    # Gaussian Knock-Offs
    suppressWarnings({gaussian_ko = create.second_order(as.matrix(test_dat))})
    dat_gko = test_dat
    dat_gko[fname] = gaussian_ko[,fname]
    mmds$ko = mmd2(dat_gko, test_dat_ref)
    mmds
  })
  print(sprintf("Done with data %i", data_id))
  res = rbindlist(res, fill = TRUE)
  print(res)
  write.csv(file = filename, res, row.names = FALSE)
})

