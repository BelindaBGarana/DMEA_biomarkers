#installing packages for biomarker analysis
#BG 20220429; last edit: BG 20220429

install.packages("reshape2", repos = "http://cran.us.r-project.org")

#install rlang
#install.packages("rlang", repos = "http://cran.us.r-project.org")

#remove BiocManager, BiocVersion and re-install
#remove.packages("BiocVersion")
#remove.packages("BiocManager")

#if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager", repos = "http://cran.us.r-project.org")}
#BiocManager::install(version='3.14') #using '3.14' instead of 'devel' because USC HPC doesn't have R v4.2
#BiocManager::install("depmap")

