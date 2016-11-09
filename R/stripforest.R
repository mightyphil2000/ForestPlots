### Functions to generate forest plots, strip forest plots, and  (soon) scatterplots
### Version: 0.53
### Date: 20160401
### Authors: cl14022, ph14916
### Latest change: Adding strip plot function

### This script does:
### 1. Read in a set of meta-analysis data, for moderate numbers of SNPs and exposures
### 2. Count the numbers of exposures (SNPs) and of outcomes
### 3. Detect the presence of summary data
### 4. Determine the appropriate plot according to the following table
###     Scatterplot: if the number of exposures and outcomes together exceeds 50 ( or go straight to summaries if it's 1:many)
###     Grouped forest plot: if there are multiple outcomes (regardless of the number of exposures)
###     Forest plot with summaries: if there are multiple exposures and a single outcome ( or multiple outcomes and one exposure)
### 5. Produces the appropriate plot in the R environment
### 6. If a filename is given, saves the appropriate plot as a png


### To do this, it uses several sub-functions
### 1. one_exp_one_out(...,summaryOnly=FALSE): produces a forest plot with multiple SNP results and multiple summaries, or optionally only the summaries
### 2. one_exp_many_out(): produces a forest plot  of summary results based on one exposure and several outcomes
### 3. many_exp_one_out(): produces a forest plot  of summary results based on several exposures and one outcome
### 4. many_exp_many_out(...,snpsOrsummary="snps"): [NOT IMPLEMENTED YET]--produces a scatterplot of effect estimates, either SNPs or summaries
### 5. output_forestplot(): produces output files and R plots

##### BEGIN forest_stripplot()


forest_stripplot <- function( inpt_df = NULL ,  eff_Col = "effect.outcome", exposure_Name="snps.test", outcome_Name="effect.estimate", point_Size = 'n.snps',  header_Col = "subcat", rowtext_Col = "trait_strict", log_ES = FALSE, exp_ES = FALSE, decrease = TRUE, se_Col = "se.outcome", ub_Col =  NULL, lb_Col = NULL) {

        # The stripplot function does not distinguish between individual-SNP estimates and summary estimates, and I am not sure that
library(ggplot2)
library(plyr)

library(reshape)
library(scales)



options(warn=-1)


if(length(lb_Col)  == 0){ ### Check if a lower_bound column has been specified or not
    inpt_df$lb <- inpt_df[,eff_Col] - qnorm(p = 0.975) * inpt_df[,se_Col]
    inpt_df$ub <- inpt_df[,eff_Col] + qnorm(p = 0.975) * inpt_df[,se_Col]
    lb_Col <- "lb"
    ub_Col <- "ub"
}


stripspacer <- function(head_col, eff_col, rowtext_col, Data_Fm, decrease = TRUE) {

    ## This function fills in a 'spacing vector', to determine placement of rows (that contain pointranges and their descriptions) in a study plot, an attribute list of textual attributes
    ## which ensures that subgroups have bold face head_cols, a content list, which gives the content for head_cols, and a row list, of the rows to be filled in with this head_col data. Further, it arranges the spacing of strips into 40-row columns of results.

    # head_col = (character), column name of the column containing the name of the particular header to use (for example, the outcome )
    # eff_col = (character), column name of the column giving the effect size of interest
    # rowtext_col = (character), column name giving the row text that distinguishes the estimates (for example, instruments used, exposures, or methods of computing the association)
    # Data_Fm = (name) name of the data frame containing the meta-analytic data
    # decrease = (logical) order the instruments by decreasing effect size y/n?


    N_study <- nrow(Data_Fm)  ## Baseline number of associations: could be any number of instruments/estimates associated with any number of outcomes
    N_cols <- ceiling(N_study / 40) ## How many 40-item columns to generate?

    # Spacing the pointranges: one row per instrument, then, buffering room added for each header
    N_space <- N_study + 1 * length(unique(Data_Fm[,head_col]))

    # Preassigning the different columns of the spacing matrix
    spacing_vec <- vector(mode = "numeric")
    content_list <- vector(mode = "character")
    attr_list <- vector(mode = "character")
    row_list <- vector(mode = "numeric")
    column_list <- vector(mode = "numeric")

    ## A shameful for loop is used to populate the matrix
    idx <- 1
    for (header in unique(Data_Fm[,head_col])) {

        # loop over headers to fill the 'spacing vector,' which determines placement of the rows in the forest plot

        inst_idxs <- which(Data_Fm[,head_col] == header) ## where are the effect sizes of interest located in the initial data frame?
        n_inst <- length(inst_idxs) ## How many instruments per header??

        columns <- ceiling((idx:(idx + 1 + n_inst))/40)
        if (length(unique(columns)) > 1){

            # if we're breaking across columns, start a new one!
            idx <- ceiling(idx/40)*40 + 1 ## restart at next multiple of 40 ( + 1 )
        }

        spacing_vec[[idx]] <- idx
        spacing_vec[[(idx + 1 + n_inst)]] <- (idx + 1 + n_inst)

        spacing_vec[((idx + 1):(idx + n_inst))] <- idx + 1:n_inst

        # In the same loop, fill a list with study names but also with head_cols to serve as the primary annotation of the forest plot
        content_list[[idx]] <- ""
        content_list[[(idx + 1 + n_inst)]] <- header

        inst_sort <- order(as.numeric(Data_Fm[inst_idxs, eff_col]), decreasing = !decrease) ## reversing direction is necessary because ascending order will naturally go with ascending rows

        content_list[((idx + 1):(idx + n_inst))] <- Data_Fm[inst_idxs[inst_sort], rowtext_col]

        # A list of attributes for the typeface so that head_cols are in bold face
        attr_list[[idx]] <- "plain" ## a space
        attr_list[[(idx + 1 + n_inst)]] <- "bold" ## header row
        attr_list[((idx + 1):(idx + n_inst))] <- "plain"  ## rowss with text

        # The most important list, contains mapping between spaces on the forest plot and rows in the meta-analytic data frame
        row_list[[idx]] <- NA ## space row, ignore
        row_list[[(idx + 1 + n_inst)]] <- NA ## Header row, ignore
        row_list[((idx + 1):(idx + n_inst))] <- inst_idxs[inst_sort] ##

        ## List of which column
        column_list[idx:(idx + 1 + n_inst)] <- rep(max(columns), (n_inst + 2)) ## Assign to highest-numbered column in the data (shorten columns you would break)

        idx <- idx + 2 + n_inst ## next index


    }

    ## repeat for summary_rows
    # returns data frame giving spacings, primary annotation (content_list), typeface attributes for the primary annotation (attr_list), and  mapping between rows in forest plot and rows in meta-analytic data frame
    return( data.frame( spacing_vec = spacing_vec%%40 + 40*(spacing_vec%%40 == 0) , content_list = content_list, row_list = row_list, attr_list = attr_list, column_list = column_list ) )
}



### 4. Generate a properly-spaced forest plot (without a summary effect) for all studies
columnforest <- function(data_Fm, column = 1, eff_col, lb_col, ub_col, se_col,  point_size = NULL , log_ES = FALSE, exp_ES = FALSE) {
    # data_Fm = (name) spaced data frame (must be from output from space_Out()) containing meta-analytic data and appropriate annotation, spacing cols etc

    # eff_col = (character) column name of effect-size column
    # lb_col = (character) column name of the confidence interval lower bound column
    # ub_col = (character) column name of the confidence interval upper bound column
    # title_text = (character) title text for main part of forest plot (will use column name if no text is supplied)
    # log_ES = (logical) convert effect sizes TO natural log scale y/n?
    # exp_ES = (logical) convert effect sizes FROM natural log scale y/n?
    # summary_rows = (integer) vector listing rows with summaries in them
    # instrument_size = (integer) vector listing the number of SNPs to have gone into each instrument


    data_Fm$lb_col <- as.numeric(data_Fm[,lb_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,lb_col]))*(exp_ES)

    data_Fm$ub_col <- as.numeric(data_Fm[,ub_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,ub_col]))*(exp_ES)

    data_Fm$eff_col <- as.numeric(data_Fm[,eff_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,eff_col]))*(exp_ES)
    data_Fm$se_col <- as.numeric(data_Fm[,se_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,se_col]))*(exp_ES)

    if(log_ES == TRUE){
        stopifnot(all(as.numeric(data_Fm[,lb_col])) >= 0) # throw error if there are negative limits

        data_Fm$lb_col <-  log(as.numeric(data_Fm[,lb_col]))

        data_Fm$ub_col <-  log(as.numeric(data_Fm[,ub_col]))

        data_Fm$eff_col <- log(as.numeric(data_Fm[,eff_col]))

        data_Fm$se_col <- log(as.numeric(data_Fm[,se_col]))
    }

   # data_Fm <- data_Fm[data_Fm[,"column_list"] == column,]
    # ggplot code to generate the forest plot using geom_segments and geom_points and to make a relatively minimal theme
    #
    #
    segment_width <- max(data_Fm[,ub_col] - data_Fm[,lb_col], na.rm = TRUE)
    label_width <- max(nchar(as.character(data_Fm[,"content_list"])))
    column_lb <- min(data_Fm[,lb_col], na.rm = TRUE) - 0.1*segment_width - 1.5*label_width
    column_ub <- max(data_Fm[,ub_col], na.rm = TRUE) + 0.05*segment_width
    column_width <- column_ub - column_lb
    bg_centre <-  mean(c(column_ub, column_lb))
    text_centre <- column_lb - 0.06*column_width

    raw_forest  <- ggplot(data = data_Fm) + geom_tile(aes(x = bg_centre, y = spacing_vec, fill= factor((spacing_vec + 1 )%%2), width = 1.3*column_width, height = 1)) + scale_fill_grey(start = 1,end = 0.85) +geom_text(data = data_Fm,aes(y = spacing_vec, x = text_centre, na.rm = TRUE , label = content_list, vjust = "bottom", fontface = attr_list, hjust = "inward")) + geom_segment(aes( y = spacing_vec, yend = spacing_vec, x = as.numeric(lb_col), xend = as.numeric(ub_col) ))  + geom_point(aes( y = spacing_vec,  x = eff_col))  + theme_minimal()+ facet_grid(~column_list) + theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),axis.line=element_line(), axis.line.y=element_blank(), legend.position = 'none', axis.line.x = element_line(size = 1), strip.background = element_blank(), strip.text = element_blank())
    # returns ggplot2 object with the (un-annotated) forest plot
    if(exp_ES == TRUE){
        raw_forest <- raw_forest  + scale_x_continuous(trans = log2_trans()) + geom_vline(xintercept = 1)
    } else {
        raw_forest <- raw_forest + geom_vline(xintercept = 0)
    }
    return(raw_forest)
}






## order and structure effect sizes and CIs for forest plot
space1 <- stripspacer( head_col = header_Col, eff_col = eff_Col, rowtext_col  = rowtext_Col, Data_Fm = inpt_df, decrease = decrease)

expand_data  <- space_Out(data_Fm = inpt_df, space_Fm = space1)


## Make the forest plot
stripfo1  <- columnforest( data_Fm = expand_data, eff_col = eff_Col, lb_col = lb_Col, ub_col = ub_Col,se_col = se_Col )
## data_Fm = expand_data; space_col = 'spacing_vec'; eff_col = eff_Col; lb_col = lb_Col; ub_col = ub_Col; log_ES = log_ES; title_text = forest_Title; exp_ES = exp_ES;summary_rows = summary_Rows;instrument_size = instrument_Size; se_col = se_Col
return(stripfo1)
}
