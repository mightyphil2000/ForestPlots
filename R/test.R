rm(list=ls())
setwd("~/ForestPlots/R/")
source("~/ForestPlots/R/forest_plots_helper_functions.R")
source("~/ForestPlots/R/many_exp_many_out.R")
source("~/ForestPlots/R/one_exp_many_out.R")

memodata <- read.delim("~/ForestPlots/inst/data/lungcancer.txt", as.is = TRUE)

Res<-memodata[memodata$subtype=="adenocarcinoma",]
# Res<-Res[Res$Phenotype!="Smoking",]

Res$Phenotype<-gsub("^\\s+|\\s+$", "", Res$Phenotype) #remove leading and trailing white space from phenotype column

Res$subcat<-"Other risk factors"  
Res$subcat[Res$Phenotype %in% c("Alpha-tocopherol",
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
Res$subcat[Res$Phenotype=="Smoking"] <-"Smoking" 
Res$subcat[Res$Phenotype %in% c("Birth weight","body mass index","height","hip circumference","waist circumference","weight","waist-to-hip ratio")]<-"Anthropometrics" 
Res$subcat[Res$Phenotype %in% c("FEV","Chronic bronchitis and chronic obstructive pulmonary disease","asthma")]<-"lung disease/function"
Res$subcat[Res$Phenotype %in% c("2hr glucose","Total cholesterol","CRP","fasting glucose","HDL cholesterol","MUFA","PUFA","Tot.FA")]<-"Metabolic biomarkers"

#Order the risk factors; plot function will use the order in the data
Res$order<-2
Res$order[Res$Phenotype=="Smoking"]<-1
Res$order[Res$subcat=="Other risk factors"]<-3
Res<-Res[order(Res$order),]


# Res$colour_col=1
one_exp_many_out(inpt_Df = Res, eff_Col = "Beta", se_Col = "SE", exposure_Name = "Phenotype", 
	outcome_Name = "subcat", 
	forest_Title = 'OR (95% CI) for lung cancer \nper SD/log odds change in risk factor', 
	outfile_Name = "~/Google Drive/MR_base/papers/lung_cancer/plots/plot.png", 
	left_Col_Names = c("Phenotype", "SNPs"), 
	left_Col_Titles = c("\n","No. of \nSNPs"), 
	right_Col_Names = c("Pvalue"),  
	right_Col_Titles = c("\nP-value"), 
	summary_List = c("IVW"), 
	exp_ES = TRUE, #exponentiate 
	weight = TRUE, #inverse variance weight
	sort_at_All=F,
	Labels=c(0.25,0.5,1,2,4,8,16)) #weight marker size by inverse variance of log odds ratio


# memodata <- read.delim("MExp-Mout-Mmethod - Sheet1.tsv", as.is = TRUE)
memodata <- read.delim("~/ForestPlots/inst/data/lungcancer.txt", as.is = TRUE)
memodata$method<-"IVW"
memodata<-memodata[order(memodata$Phenotype,decreasing=T),]

lb<-memodata$Beta-1.96*memodata$SE
ub<-memodata$Beta+1.96*memodata$SE

Limits<-c(min(lb),max(ub))
Limits[2]-Limits[1]
From=Limits[1]
To=Limits[2]
exp(seq(From,To,by=0.6931472))

exp(seq(0.06, o,by=0.6931472))
0.25/2/2
log(2)
log(4)
log(8)-log(4)
log(4)-log(2)
# png("png.png",width=3200,height=3200)

memotestplot <- many_exp_many_out(inpt_Df = memodata, eff_Col = "Beta", 
	se_Col = "SE", ub_Col = NULL, lb_Col = NULL, exposure_Name = "subtype", 
	outcome_Name = "Phenotype",forest_Title = 'Effect size in log units, presented as a 95% confidence interval', 
	outfile_Name = "~/Google Drive/MR_base/papers/lung_cancer/plots/annot_memo.png", 
	left_Col_Names = c('subtype'), 
	left_Col_Titles = c(""), 
	right_Col_Names = NULL, 
	right_Col_Titles = NULL, 
	log_ES = FALSE, 
	exp_ES = TRUE, 
	summary_List = c('IVW'), 
	weight = TRUE,
	sort_at_All=T,
	Labels=c(0.06,0.12,0.25,0.5,1.0,2.0,4.0,8.0,16,32,64))
