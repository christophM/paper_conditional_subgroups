# =============================================================================
# Create the data overview table for LaTeX 
# =============================================================================

library("rjson")
library("here")
devtools::load_all()
data_dir = here("data")
openml_dir = here("data/openml")


data_ids = unique(gsub("\\.(arff|json)", "", list.files(openml_dir))) 
data_info = lapply(data_ids, function(data_id){
  dat = farff::readARFF(sprintf("%s/%s.arff", openml_dir, data_id), show.info = FALSE)
  j = fromJSON(file = sprintf("%s/%s.json", openml_dir, data_id))
  data_name = j$name
  data.frame(name = data_name, n = nrow(dat), p = ncol(dat))
})

data_info = rbindlist(unique(data_info))
data_info = t(data_info)
cap = "We selected data sets from OpenML \\cite{vanschoren2014openml, Casalicchio2017} having 1000 to 8000 instances and a maximum of 50 numerical features. We excluded data sets with categorical features, since ALE cannot handle them."
lab = "tab:datasets"
colnames(data_info) = data_info[1,]
rownames(data_info) = c("", "No. of rows", "No. of features") 
rnm_dat = c("space_ga" = "space",
            "wind" = "wind",
            "satellite_image" = "satellite",
            "pollen" = "pollen",
            "wine_quality" = "wine",
            "quake" = "quake")
colnames(data_info) = rnm_dat[colnames(data_info)]


xtab = xtable(data_info[2:3,], caption = cap, label = lab, align = c("l", rep("r", times = ncol(data_info))))

print(xtab, file = sprintf("%s/data-info.tex", fig_dir))

