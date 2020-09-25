Code to recreate main results from https://www.medrxiv.org/content/10.1101/2020.09.08.20190918v1

All code is implemented in R and was tested on R version 4.0.2 (2020-06-22).
All packages required to run code are listed in Functions > loadpackages .
Install/run time is approximately 10 minutes. 

"Main.R" recreates Fig. 1

"NYepiestim.R" estimates R effective for New York using nytimes county-level data. NB/ this analysis was performed on the 21st July and results are recorded in "NYestimR.RData"

"NYq.RData" is weekly specific humidity data for New York summarized from ERA5
