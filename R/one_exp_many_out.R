

##### Begin one_exp_many_out()
one_exp_many_out <- function( inpt_Df = NULL , eff_Col = "effect.outcome", exposure_Name="snps.test", outcome_Name="effect.estimate",forest_Title = 'Effect', outfile_Name = 'annot_FP.pdf', left_Col_Names=c("outcome", "snps.test", "chr.names"), left_Col_Titles = c("Outcome", "SNP", "Chr"), right_Col_Names = NULL, right_Col_Titles = NULL, log_ES = FALSE, exp_ES = FALSE, se_Col = "se.outcome", ub_Col =  NULL, lb_Col = NULL, summary_List = c("Inverse variance weighted"), weight = TRUE, sort_at_All = TRUE, N_Breaks = 3) {
  
  library(grid)
  library(ggplot2)
  library(plyr)
  library(reshape2)
  library(scales)
  
  source("./forest_plots_helper_functions.R")
  # inpt_Df = (character) filename of the data
  # eff_Col = (character) character giving the MR effect column
  # exposure_Name = (character) name of the SNP Name/Summary Name column
  # outcome_Name = (character) name of the column giving the outcome
  # forest_Title = (character) title of the forest plot
  # outfile_Name = (character) name of the file to save the png of the forest plot to
  # left_Col_Names = (character) vector of the column names of the LHS annotations to be used from the original data file
  # left_Col_Titles = (character) vector of headings for the forest plot columns of the LHS
  # right_Col_Names = (character) vector of the column names of the RHS annotations to be used from the original data file
  # right_Col_Titles = (character) vector of headings for the forest plot columns of the RHS
  # log_ES = (logical) log-transform the effect size and confidence interval limits Y/N?
  # exp_ES = (logical) exponential-transform the effect size and confidence interval limits Y/N?
  # se_Col = (character) column name of the standard error column in the original data
  # ub_Col =  (character) column name of the upper bound of CI column in the original data
  # lb_Col = (character) column name of the lower bound of CI column in the original data
  # summary_List = (character) names of the effect estimate methods to be plotted (to select which summaries to plot)
  # weight = (logical) column name of column containing instrument sizes (to eliminate single-snp results) 
  
  options(warn=-1)
  keep_rows <- NULL
  for (i in summary_List){
    keep_rows <- c(which(i == inpt_Df[,"method"]), keep_rows)
  }
  number_Exp <- length(unique(inpt_Df[, exposure_Name])) # Number of exposures
  number_Out <- length(unique(inpt_Df[, outcome_Name])) # Number of outcomes
  total_Rows <- nrow(inpt_Df)
  
  inpt_Df <- inpt_Df[keep_rows,]
  
  if(length(lb_Col)  == 0){
    inpt_Df$lb <- inpt_Df[,eff_Col] - qnorm(p = 0.975) * inpt_Df[,se_Col]
    inpt_Df$ub <- inpt_Df[,eff_Col] + qnorm(p = 0.975) * inpt_Df[,se_Col]
    lb_Col <- "lb"
    ub_Col <- "ub"
  }
  ## inpt_Df[,"orci"] <- paste0(" (", signif(inpt_Df[,lb_Col],3),", ", signif(inpt_Df[,ub_Col],3), ")")
  
  #
  space1 <- spacer( major = outcome_Name, sort_var = eff_Col, data_Fm = inpt_Df, sort_at_all = sort_at_All)
  
  expand_data  <- space_Out(data_Fm = inpt_Df, space_Fm = space1)
  
  
  ## Make the forest plot
  fo1  <- makeforest( data_Fm = expand_data, eff_col = eff_Col, lb_col = lb_Col, ub_col = ub_Col, log_ES = log_ES, title_text = forest_Title, exp_ES = exp_ES, se_col = se_Col, se_Weight = weight, n_Breaks = N_Breaks)
  
  ## Construct left-hand-side annotations
  left <- anot_side( data_Fm = expand_data, col_names = left_Col_Names, title_list = left_Col_Titles, notitles = FALSE, content_Name = outcome_Name)
  
  ## Construct right-hand-side annotations
  ## right <- anot_side(data_Fm = expand_data, col_names = c("orci", right_Col_Names), title_list = c("CI Limits", right_Col_Titles), notitles = FALSE, left = FALSE )
  right <- NULL
  if (length(right_Col_Names) > 0){
    right <- anot_side(data_Fm = expand_data, col_names = right_Col_Names, title_list = right_Col_Titles, notitles = FALSE, left = FALSE )
  }
  
  options(warn=0)
  output_forestplot(forest = fo1, left =  left, right = right, outfile_Name = outfile_Name, nrows = nrow(expand_data), ncols = length(right_Col_Names) + length(left_Col_Names), plotType = "oEmO")
  return(list(forest = fo1, left = left, right = right))
  
}
##### END one_exp_many_out()
