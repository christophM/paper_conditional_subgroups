# =============================================================================
# Download datasets from OpenML
# =============================================================================

library("OpenML")
library("rjson")
library("here")
devtools::load_all()
# Random seed is needed bc. random feature selection for data fidelity experiment happens here
set.seed(42)

# All datasets are downloaded into "raw"
cc18dir_raw = sprintf("%s/cc18/raw/", data_dir)
# All processed datasets are stored in "processed"
cc18dir_processed= sprintf("%s/cc18/processed/", data_dir)

cc18list = jsonlite::fromJSON(sprintf("%s/cc18.json", data_dir))
data_ids = cc18list[[1]]$data$data_id

targets = lapply(data_ids, function(data_id) {
  j = jsonlite::fromJSON(sprintf("https://www.openml.org/d/%s/json", data_id))
  obj = toJSON(j)
	filename = sprintf("%s/%s.json", cc18dir_raw, data_id)
	if (!file.exists(filename)) {
		write(obj, file = filename)
	}
	fn = sprintf("%s/%s.arff", cc18dir_raw, data_id)
	if (!file.exists(fn)) {
		download.file(j$url, fn)
	}
})


# Remove factors
data_ids = unique(gsub("\\.(arff|json)", "", list.files(cc18dir_raw)))


datasets = lapply(data_ids, function(data_id) {
  json_file = sprintf("%s/%s.json", cc18dir_raw, data_id)
  j = fromJSON(file = json_file)
  data_name = j$name
  print(sprintf("Processing data %s", data_id))
  # Load data
  dat_file = sprintf("%s/%s.arff", cc18dir_raw, data_id)
  dat = farff::readARFF(dat_file, show.info = FALSE)
  ncol_orig = ncol(dat)
  # Some column names caused troubles
  colnames(dat) = gsub("[^0-9a-zA-Z]+", "", colnames(dat))
  print(dim(dat))
  target =  gsub("[^0-9a-zA-Z]+", "", j$default_target_attribute)

  # remove factors
  is_num = which(sapply(dat, is.numeric))
  tindex = which(colnames(dat) == target)
  keep = unique(c(is_num, tindex))
  print(sprintf("%i/%i columns left",length(keep), ncol(dat)))
  dat = dat[keep]

  # remove features that have no variability
  if (ncol(dat) > 1) {
    constant_feats =  sapply(dat, function(x) length(unique(x)) < 2)
    dat = dat[!constant_feats]
  }

  dat = na.omit(dat)

  # if no features left, skip
  if((ncol(dat) <= 7) | (ncol(dat) > 500)) {
    print(sprintf("Data with id %s has only non-numerical featurs", data_id))
    print(sprintf("Deleted data with id %s", data_id))
    return(data.frame())
    dat = NULL
  } else {
    filename = sprintf("%s/%s.arff", cc18dir_processed, data_id)
    # Mix feature order and add target at the end
    # This is so that in data-fidelity we can always use the first ten features
		# Last column is the target
    ord = c(sample(1:(ncol(dat) - 1)), ncol(dat))
    dat = dat[ord]
    farff::writeARFF(path = filename, dat, overwrite = TRUE)
  }


  data.frame(openmlid = data_id,
             name = data_name,
             nrows = nrow(dat),
             ncols = ncol(dat),
             ncols_before = ncol_orig)
})


dataset_info = data.table::rbindlist(datasets)
write.csv(file = sprintf("%s/cc18-infos.csv", fig_dir), dataset_info, row.names = FALSE)

