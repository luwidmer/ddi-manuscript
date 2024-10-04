renv::restore()

required_packages <- c(
  "OncoBayes2",
  "abind",
  "tidyverse",
  "ggplot2",
  "egg",
  "RColorBrewer",
  "here",
  "grid",
  "brms",
  "tinytex",
  "knitr",
  "posterior",
  "bayesplot",
  "tidybayes"
)

if (!all(required_packages %in% installed.packages()))
{
  install.packages(required_packages[!(required_packages %in% installed.packages())])
}

for (pkg in required_packages)
{
  library(package = pkg, character.only = TRUE)
}

if (!tinytex::is_tinytex() && length(Sys.which('tlmgr')[[1]]) == 0)
{
  tinytex::install_tinytex()
}
