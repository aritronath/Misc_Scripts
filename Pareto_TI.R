### Install ParetoTI and python dependencies
BiocManager::install("vitkl/ParetoTI", dependencies = c("Depends", "Imports", "LinkingTo"))
library(ParetoTI)
ParetoTI::install_py_pcha(method = "conda", 
                          extra_packages = c("tensorflow", "tensorflow-probability",
                                             "pandas", "keras", "h5py",
                                             "geosketch", "pydot", "scikit-learn==0.20",
                                             "umap-learn"))
#check everything was insatlled and discoverable
reticulate::py_discover_config("py_pcha")



#### Example ----------------------------------------------
library(ParetoTI)
library(ggplot2)

# Generate random data that fits into the triangle (3D)
set.seed(4355)
archetypes = generate_arc(arc_coord = list(c(5, 0, 4), c(-10, 15, 0), c(-30, -20, -5)),
                          mean = 0, sd = 1)
data = generate_data(archetypes$XC, N_examples = 1e4, jiiter = 0.04, size = 0.9)
# Fit polytope to those data
arc_data = fit_pch(data, noc = as.integer(3), delta = 0)

# Show results as interactive 3D scatterplot using plotly
plot_arc(arc_data = arc_data, data = data,
         which_dimensions = 1:3)
# Plot static 2D scatterplot using ggplot2
plot_arc(arc_data = arc_data, data = data,
         which_dimensions = 1:2) +
  theme_bw()

# Plot data as 2D density rather than scatterplot
plot_arc(arc_data = arc_data, data = data,
         which_dimensions = 1:2, geom = ggplot2::geom_bin2d) +
  theme_bw()

# Project to UMAP coordinates (3D -> 2D)
arc_umap = arch_to_umap(arc_data, data, which_dimensions = 1:2,
                        method = c("naive", # implemented in R and slow
                                   "umap-learn")) # requires python module
plot_arc(arc_data = arc_umap$arc_data, data = arc_umap$data,
         which_dimensions = 1:2) +
  theme_bw()
# Project to tSNE coordinates (3D -> 2D, requires Rtsne package)
arc_tsne = arch_to_tsne(arc_data, data, which_dimensions = 1:2)
plot_arc(arc_data = arc_tsne$arc_data, data = arc_tsne$data,
         which_dimensions = 1:2) +
  theme_bw()

## --------------------------------------
dat <- raw2
arc_data = fit_pch(dat, noc = as.integer(3), delta = 0)

plot_arc(arc_data = arc_data, data = dat,
         which_dimensions = 1:3)

# Plot static 2D scatterplot using ggplot2
plot_arc(arc_data = arc_data, data = dat,
         which_dimensions = 1:2) +
  theme_bw()

# Plot data as 2D density rather than scatterplot
plot_arc(arc_data = arc_data, data = data,
         which_dimensions = 1:2, geom = ggplot2::geom_bin2d) +
  theme_bw()
