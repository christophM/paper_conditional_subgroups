library(xtable)
devtools::load_all()

infos = read.csv(sprintf("%s/cc18-infos.csv", fig_dir))
colnames(infos) = c("OpenML ID", "Name", "No. Obs.", "No. numerical feat.", "No. feat.")
lab = "tab:cc18"
cap = "Overview of OpenML CC18 data sets used for the data fidelity experiment."
tab = xtable(infos, label = lab, caption = cap)
print(tab, include.rownames = FALSE, file = sprintf("%s/cc18-info.tex", fig_dir), size = "small")
