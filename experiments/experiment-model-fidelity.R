# =============================================================================
# Model Fidelity experiment
# =============================================================================

library("mlr3")
library("mlr3learners")
library("mlr3pipelines")
library("here")
library("ranger")
library("parallel")
library("OpenML")
library("farff")
library("rjson")
library("rpart")
devtools::load_all()

# number of cores to be used
NC = 7

mserror = function(x,y){
  dif = (x - y)
  l = length(dif)
  dif = na.omit(dif)
  if(l > length(dif)) print(sprintf("Omitted %i/%i NA diffs", l - length(dif), l))
  mean(dif^2)
}

get_pdp_loss = function(mod, train_dat, test_dat, target, fname){
  pred = Predictor$new(mod, data = train_dat, y = target)
  nunique = length(unique(train_dat[[fname]]))
  # Case when some node has only one unique xj value
  if (nunique == 1){
    print(sprintf("Feature %s has only one unique value", fname))
    prediction = mean(pred$predict(train_dat)[[1]])
  } else {
    pdp = FeatureEffect$new(pred, feature = fname, method = "pdp")
    prediction = pdp$predict(test_dat)
  }
  mserror(prediction, test_dat[[target]])
}

get_gcpdp_loss = function(conds, mod, dat, fname, train_ids, test_ids, target){
  cond = conds[[fname]]
  train_nodes = cond$cnode(dat[train_ids,])$node
  test_nodes = cond$cnode(dat[test_ids,])$node
  mses = sapply(unique(test_nodes), function(wi) {
   train_dat = dat[train_ids,][train_nodes %in% wi,]
   test_dat =  dat[test_ids,][test_nodes %in% wi,]
   get_pdp_loss(mod, train_dat = train_dat, test_dat = test_dat,
                target = target, fname = fname) * nrow(test_dat)
  })
  gcpdp_loss = sum(mses) / nrow(dat[test_ids,])
}

data_dir = here("data")
openml_dir = here("data/openml")
set.seed(10)
test_size = 0.3


data_ids = unique(gsub("\\.(arff|json)", "", list.files(openml_dir)))
learners = list(lrn("regr.ranger"), lrn("regr.lm"), lrn("regr.kknn"))
lrn_ids = seq_along(learners)

experiments = expand.grid(lrn_id = lrn_ids, data_id = data_ids)
# shuffle experiments
experiments = experiments[sample(1:nrow(experiments)), ]
experiments$data_id = as.character(experiments$data_id)


res = mclapply(1:nrow(experiments),
               mc.cores = NC,
               function(i) {
  data_id = experiments[i, "data_id"]
  lrn_id = experiments[i, "lrn_id"]
  filename = sprintf("%s/model-fidelity/mod-fid-%s-%s.csv", res_dir, data_id, lrn_id)
  if(file.exists(filename)) {
    print(sprintf("Skipping %s", filename))
    return()
  }
  print(sprintf("Processing data %s and learner %s", data_id, lrn_id))
  dat = farff::readARFF(sprintf("%s/%s.arff", openml_dir, data_id), show.info = FALSE)
  # Some column names caused troubles
  colnames(dat) = gsub("[^0-9a-zA-Z]+", "", colnames(dat))
  print(dim(dat))
  j = fromJSON(file = sprintf("%s/%s.json", openml_dir, data_id))
  target = j$default_target_attribute
  target = gsub("[^0-9a-zA-Z]+", "", target)
  task = TaskRegr$new(id = as.character(data_id), backend = dat, target = target)
  # Split for training+tuning and for evaluation+iml
  train_ids = sample(1:nrow(dat), size = (1 - test_size) * nrow(dat))
  test_ids = setdiff(seq_len(nrow(dat)), train_ids)

  learner = learners[[lrn_id]]$clone()
  mod = learner$train(task, row_ids = train_ids)
  dat = data.frame(dat)
  dat_sub = dat[train_ids, task$feature_names]
  conds1 = fit_conditionals(dat_sub, ctrl = ctrl(1))
  conds2 = fit_conditionals(dat_sub, ctrl = ctrl(2))
  conds5 = fit_conditionals(dat_sub, ctrl = ctrl(5))
  conds10 = fit_conditionals(dat_sub, ctrl = ctrl(10))

  conds1c = fit_conditionals(dat_sub, ctrl = rpart.control(maxdepth=1), type = "cart")
  conds2c = fit_conditionals(dat_sub, ctrl = rpart.control(maxdepth=2), type = "cart")
  conds5c = fit_conditionals(dat_sub, ctrl = rpart.control(maxdepth=5), type = "cart")
  conds10c = fit_conditionals(dat_sub, ctrl = rpart.control(maxdepth=10), type = "cart")


  pred = Predictor$new(mod, data = dat[train_ids, ], y = target)
  mean_pred = mean(pred$predict(dat[train_ids,])[[1]])
  res = lapply(task$feature_names,  function(fname){
    print(fname)
    pdp_loss = get_pdp_loss(mod, train_dat = dat[train_ids, ], test_dat = dat[test_ids,], target = target, fname = fname)
    ale = FeatureEffect$new(pred, feature = fname, method = "ale")
    prediction = ale$predict(dat[test_ids,]) + mean_pred
    ale_loss = mserror(prediction, dat[test_ids, target])

    gcpdp_loss1  =  get_gcpdp_loss(conds1, mod, dat, fname, train_ids, test_ids, target)
    gcpdp_loss2  =  get_gcpdp_loss(conds2, mod, dat, fname, train_ids, test_ids, target)
    gcpdp_loss5  =  get_gcpdp_loss(conds5, mod, dat, fname, train_ids, test_ids, target)
    gcpdp_loss10 =  get_gcpdp_loss(conds10, mod, dat, fname, train_ids, test_ids, target)

    gcpdp_loss1c  =  get_gcpdp_loss(conds1c, mod, dat, fname, train_ids, test_ids, target)
    gcpdp_loss2c  =  get_gcpdp_loss(conds2c, mod, dat, fname, train_ids, test_ids, target)
    gcpdp_loss5c  =  get_gcpdp_loss(conds5c, mod, dat, fname, train_ids, test_ids, target)
    gcpdp_loss10c =  get_gcpdp_loss(conds10c, mod, dat, fname, train_ids, test_ids, target)


    data.frame(pdp_loss, ale_loss,
               feature = fname,
               gcpdp_loss1,
               gcpdp_loss2,
               gcpdp_loss5,
               gcpdp_loss10,
               gcpdp_loss1c,
               gcpdp_loss2c,
               gcpdp_loss5c,
               gcpdp_loss10c,
               dat_id = task$id,
               dat_name = j$name, learner = learner$id)
})
  res = rbindlist(res)
  write.csv(file = filename, res)
})
