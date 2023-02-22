source("renv/activate.R")

source("src/install_packages.R")
source("src/run_data_scenarios.R")
source("compile_manuscript.R")

compile_manuscript()

open_pdf("manuscript.pdf")
