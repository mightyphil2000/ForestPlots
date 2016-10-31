## Mexp:1out
source("./many_exp_one_out.R")
meoodata <- read.delim("./resultsfilesforplots/Mexp-1Out - Sheet1.tsv", as.is = TRUE)
meooplot <- many_exp_one_out( inpt_Df = meoodata,  eff_Col = "b", exposure_Name="exposure", outcome_Name="id.outcome", forest_Title = 'Effect size estimates on log scale with 95% confidence', outfile_Name = 'annot_meoo.png', left_Col_Names=c("method", "b", "se"), left_Col_Titles = c("method", "b", "se"), right_Col_Names = NULL, right_Col_Titles = NULL, log_ES = FALSE, exp_ES = FALSE,  se_Col = "se", ub_Col =  NULL, lb_Col = NULL, summary_list = c("Wald ratio", "Maximum likelihood"))
    
## 1exp:1out - SNP; current status on this is that 
source("./one_exp_one_out.R") ## there really isn't an example of this, this is merely for summaries
oeoodataSNPs <- read.delim("./resultsfilesforplots/1exp1out - snps.tsv", as.is = TRUE)
oeooSNPplot <- one_exp_one_out( inpt_Df = oeoodataSNPs,  eff_Col = "b", exposure_Name="SNP", outcome_Name="outcome", instrument_Size = 'nsnp',forest_Title = 'Effect size in log units, presented as a 95% confidence interval', outfile_Name = 'annot_oeoo.png', left_Col_Names=NULL, left_Col_Titles = NULL, right_Col_Names = NULL, right_Col_Titles = NULL, log_ES = FALSE, exp_ES = FALSE, se_Col = "se", ub_Col =  NULL, lb_Col = NULL, onetype_Only = TRUE, summary_List = c("Wald ratio"))

## 1exp:1out - Summaries 
## Summaries done separately means that 
oeoodataSUMs <- read.delim("./resultsfilesforplots/1exp1out - summaries.tsv", as.is = TRUE)
oeooSUMplot <- one_exp_one_out( inpt_Df = oeoodataSUMs,  eff_Col = "b", exposure_Name="method", outcome_Name="outcome", instrument_Size = 'nsnp',forest_Title = 'Effect size in log units, presented as a 95% confidence interval', outfile_Name = 'annot_oeooS.png', left_Col_Names=NULL, left_Col_Titles = NULL, right_Col_Names = NULL, right_Col_Titles = NULL, log_ES = FALSE, exp_ES = FALSE, se_Col = "se", ub_Col =  NULL, lb_Col = NULL, onetype_Only = TRUE, summary_List = c("MR Egger", "Maximum likelihood", "Inverse variance weighted", "Weighted median"))

  ## 1exp:Mout
  source("./one_exp_many_out.R")
  oemedata <- read.delim("./resultsfilesforplots/1exp-Mout - Sheet1.tsv", as.is = TRUE)
  oemeplot <- one_exp_many_out(inpt_Df = oemedata , eff_Col = "b", exposure_Name = "id.exposure", outcome_Name = "outcome", forest_Title = 'Effect size in log units, presented as a 95% confidence interval', outfile_Name = 'annot_oemo.png', left_Col_Names = c('nsnp', 'ncontrol'), left_Col_Titles = NULL, right_Col_Names = c("Q", "Q_df", "Q_pval"), right_Col_Titles = c("Q", "Q_sd", "p"), log_ES = FALSE, exp_ES = FALSE, se_Col = "se", ub_Col =  NULL, lb_Col = NULL, sort_at_All = TRUE, N_Breaks = 7)
  
  ## Mexp:Mout
source("many_exp_many_out.R")
memodata <- read.delim("./resultsfilesforplots/MExp-Mout-Mmethod - Sheet1.tsv", as.is = TRUE)
memotestplot <- many_exp_many_out(inpt_Df = memodata, eff_Col = "b", se_Col = "se", ub_Col = NULL, lb_Col = NULL, exposure_Name = "exposure", outcome_Name = "outcome",forest_Title = 'Effect size in log units, presented as a 95% confidence interval', outfile_Name = "annot_memo.png", left_Col_Names = c('exposure', 'method'), left_Col_Titles = NULL, right_Col_Names = NULL, right_Col_Titles = NULL, log_ES = FALSE, exp_ES = FALSE, summary_List = c('Inverse variance weighted', 'Maximum likelihood', 'MR Egger', 'Weighted median'), weight = TRUE, sort_at_All = FALSE, N_Breaks = 3)



## Method4 Exp1
## "Space"
## Method1 Exp2
## Method2 Exp2
## Method3 Exp2
## Method4 Exp2

## I'm not sure if this requires a new function or perhaps just requires
## an extension to the existing "1 to many" functions with an argument
## for a

