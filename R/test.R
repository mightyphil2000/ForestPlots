rm(list=ls())
setwd("~/ForestPlots/R/")
source("~/ForestPlots/R/forest_plots_helper_functions.R")
source("~/ForestPlots/R/many_exp_many_out.R")
source("~/ForestPlots/R/one_exp_many_out.R")

memodata <- read.delim("~/ForestPlots/inst/data/lungcancer.txt", as.is = TRUE)


# memodata<-memodata[memodata$Phenotype!="Smoking",]

memodata$Phenotype<-gsub("^\\s+|\\s+$", "", memodata$Phenotype) #remove leading and trailing white space from phenotype column

memodata$subcat<-"Other risk factors"  
memodata$subcat[memodata$Phenotype %in% c("Alpha-tocopherol",
"Beta-Carotene" ,
"Bilirubin",
"Calcium",
"Homocysteine",
"Iron",
"kynurenine",
"Magnesium",
"methionine",
"Retinol",
"Selenium",
"tryptophan",
"Vitamin B12",
"Vitamin B6",
"Vitamin D",
"Vitamin E")]<-"Vitamins and minerals"
memodata$subcat[memodata$Phenotype=="Smoking"] <-"Smoking" 
memodata$subcat[memodata$Phenotype %in% c("Birth weight","body mass index","height","hip circumference","waist circumference","weight","waist-to-hip ratio")]<-"Anthropometrics" 
memodata$subcat[memodata$Phenotype %in% c("FEV","Chronic bronchitis and chronic obstructive pulmonary disease","asthma")]<-"lung disease/function"
memodata$subcat[memodata$Phenotype %in% c("2hr glucose","Total cholesterol","CRP","fasting glucose","HDL cholesterol","MUFA","PUFA","Tot.FA")]<-"Metabolic biomarkers"

memodata$Phenotype[memodata$Phenotype=="Chronic bronchitis and chronic obstructive pulmonary disease"]<-"CB/COPD"
#Order the risk factors; plot function will use the order in the data
memodata$order<-2
memodata$order[memodata$Phenotype=="Smoking"]<-1
memodata$order[memodata$subcat=="Other risk factors"]<-3
memodata<-memodata[order(memodata$order),]

Res<-memodata[memodata$subtype=="adenocarcinoma",]

# # Res$colour_col=1
# one_exp_many_out(inpt_Df = Res, eff_Col = "Beta", se_Col = "SE", exposure_Name = "Phenotype", 
# 	outcome_Name = "subcat", 
# 	forest_Title = 'OR (95% CI) for lung cancer \nper SD/log odds change in risk factor', 
# 	outfile_Name = "~/Google Drive/MR_base/papers/lung_cancer/plots/plot.png", 
# 	left_Col_Names = c("Phenotype", "SNPs"), 
# 	left_Col_Titles = c("\n","No. of \nSNPs"), 
# 	right_Col_Names = c("Pvalue"),  
# 	right_Col_Titles = c("\nP-value"), 
# 	summary_List = c("IVW"), 
# 	exp_ES = TRUE, #exponentiate 
# 	weight = TRUE, #inverse variance weight
# 	sort_at_All=F,
# 	Labels=c(0.25,0.5,1,2,4,8,16)) #weight marker size by inverse variance of log odds ratio


# memodata<-memodata[order(memodata$Phenotype,decreasing=T),]

memotestplot <- many_exp_many_out(inpt_Df = memodata, eff_Col = "Beta", 
	se_Col = "SE", ub_Col = NULL, lb_Col = NULL, 
	exposure_Name = "subtype", # the results are coloured by this variable, e.g. the lung cancer subtypes 
	outcome_Name = "subcat", # this stratifies the results by phenotype subcategory, e.g. anthropometric , glycemic, etc
	forest_Title = 'Effect size in log units, presented as a 95% confidence interval', 
	outfile_Name = "annot_memo.png", 
	left_Col_Names = c("Phenotype","subtype"), 
	left_Col_Titles = c("",""), 
	right_Col_Names = NULL, 
	right_Col_Titles = NULL, 
	log_ES = FALSE, 
	exp_ES = TRUE, #exponentiate log odds ratio
	summary_List = c('IVW'), #name of method
	weight = TRUE, #weight the effect sizes by the inverse of the variance of the log odds ratio
	sort_at_All=FALSE, #when FALSE uses the order in the data; TRUE means automatically sort the data by effect size 
	Labels=c(0.06,0.12,0.25,0.5,1.0,2.0,4.0,8.0,16,32,64)) #labels for the x axis
