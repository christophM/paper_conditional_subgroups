# =============================================================================
# Download datasets from OpenML 
# =============================================================================

library("OpenML")
library("rjson")
library("here")
devtools::load_all()
data_dir = here("data")
openml_dir = here("data/openml")

# Load all datasets with numerical features, more than 3 features, less than 50 features
ds = listOMLTasks()
ds = ds[ds$task.type == "Supervised Regression",]
ds = ds[!is.na(ds$number.of.instances),]
ds = ds[ds$number.of.instances > 1000,]
ds = ds[ds$number.of.instances < 8000,]
ds = ds[ds$number.of.numeric.features  < 50,]
ds = ds[ds$number.of.numeric.features  >= 4,]
# Causes Problems
ds = ds[ds$format != "Sparse_ARFF",]
# ALE cannot handle categorical features
ds = ds[ds$number.of.symbolic.features == 0,]
data_ids = unique(ds$data.id)


# data that already exists with higher version (= duplicated) 
data_ids = setdiff(data_ids, 209)
# data without target
data_ids = setdiff(data_ids, c(4553, 504))

targets = lapply(data_ids, function(data_id) {
  j = jsonlite::fromJSON(sprintf("https://www.openml.org/d/%i/json", data_id))
  obj = toJSON(j)
	filename = sprintf("%s/%i.json", openml_dir, data_id)
	if (!file.exists(filename)) {
		write(obj, file = filename) 
	}
	fn = sprintf("%s/%s.arff", openml_dir, data_id)
	if (!file.exists(fn)) {
		download.file(j$url, fn)
	}
})


