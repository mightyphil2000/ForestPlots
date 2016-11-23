
source("./forest_plots_helper_functions.R")

## 1exp:Mout
source("./one_exp_many_out.R")


##########################################################################################
newtestdata <- read.delim("../inst/data/lungcancer.txt", as.is = TRUE)
newtestplot <- one_exp_many_out(inpt_Df = newtestdata, eff_Col = "Beta", exposure_Name = "Phenotype", outcome_Name = "subtype", forest_Title = 'Effect size in log units, presented as a 95% confidence interval', outfile_Name = "weighted_lung.wmf", left_Col_Names = c("SNPs"), left_Col_Titles = c("# SNPs"), right_Col_Names = c("Pvalue"), right_Col_Titles = NULL, se_Col = "SE", summary_List = c("IVW"), exp_ES = TRUE, weight = TRUE ) 
