
source("./forest_plots_helper_functions.R")

## 1exp:Mout
source("./one_exp_many_out.R")


##########################################################################################
newtestdata <- read.delim("dis_lpa.txt", as.is = TRUE)
newtestplot <- one_exp_many_out(inpt_Df = newtestdata, eff_Col = "b", exposure_Name = "exposure", outcome_Name = "subcat", forest_Title = 'Effect size in log units, presented as a 95% confidence interval', outfile_Name = "weighted_lung.png", left_Col_Names = c("trait_strict", "ncase"), left_Col_Titles = c("Outcome","# Cases"), right_Col_Names = NULL, right_Col_Titles = NULL, se_Col = "se", summary_List = c("Wald ratio"), exp_ES = TRUE, weight = TRUE ) 

