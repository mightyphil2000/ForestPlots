### Functions to generate forest plots, strip forest plots, and  (soon) scatterplots
### Version: 0.55
### Date: 20160906
### Authors: cl14022, ph14916

### Given a user's classification of columns in a data frame, where 
#### each column is either: effect size, ordering, LHS, RHS, do not plot;
#### and the data frame contains effects of one exposure on one outcome,
#### possibly with different association and summary measures

## This script does:
### 1. Order rows by a particular column
### 2. Restrict  to only the appropriate columns to be plotted 
### 3. Given an ordering of rows, determine an ordering of columns to plot
### 4. Determine row and column headers (row headers from columns, column headers user-supplied)
### 5. Provide spacing for rows
### 6. Provide spacing for columns 
### 7. Size canvas for plot
### 8. Plot LHS, centre-forest, RHS, with given spacing and headers
### 9. Output to file

one_exp_one_out <- 
    function( inpt_Df = NULL ,  
              eff_Col = "effect.outcome", 
              exposure_Name="snps.test", 
              outcome_Name="effect.estimate", 
              instrument_Size = 'n.snps',
              forest_Title = 'Effect on a scale given with 95% Confidence Interval', 
              outfile_Name = 'annot_FP.pdf', 
              left_Col_Names=c("outcome", "snps.test", "chr.names"), 
              left_Col_Titles = c("Outcome", "SNP", "Chr"), 
              right_Col_Names = NULL, 
              right_Col_Titles = NULL, 
              log_ES = FALSE, 
              exp_ES = FALSE, 
              se_Col = "se.outcome", 
              ub_Col =  NULL, 
              lb_Col = NULL, 
              onetype_Only=FALSE, 
              summary_List = c("Inverse variance weighted"),
              sort_at_All=F)
              ){
        
        library(ggplot2)
        library(plyr)
        library(reshape2)
        library(scales)
        
        ## inpt_Df = (character) filename of the data
        ## eff_Col = (character) character giving the MR effect column (used for sorting and required for forest plot)
        ## exposure_Name = (character) name of the SNP Name/Summary Name column
        ## outcome_Name = (character) name of the column giving the outcome
        ## instrument_Size = (character) name of the column giving the number of SNPs in each instrument
        ## forest_Title = (character) title of the forest plot
        ## outfile_Name = (character) name of the file to save the png of the forest plot to
        ## left_Col_Names = (character) vector of the column names of the LHS annotations to be used from the original data file
        ## left_Col_Titles = (character) vector of headings for the forest plot columns of the LHS
        ## right_Col_Names = (character) vector of the column names of the RHS annotations to be used from the original data file
        ## right_Col_Titles = (character) vector of headings for the forest plot columns of the RHS
        ## log_ES = (logical) log-transform the effect size and confidence interval limits Y/N?
        ## exp_ES = (logical) exponential-transform the effect size and confidence interval limits Y/N?
        ## se_Col = (character) column name of the standard error column in the original data (either SE or CI bounds are needed for the forest plot)
        ## ub_Col =  (character) column name of the upper bound of CI column in the original data
        ## lb_Col = (character) column name of the lower bound of CI column in the original data
        ## onetype_Only = (logical) plot only summary data or single-SNPs and summaries?
        ## summary_List = (character) names of summary statistics to include (a subset of the unique values in exposure_Name)
        
        options(warn=-1)
        
        summary_Rows <- which(inpt_Df[, instrument_Size] > 1)
        
        if (onetype_Only == FALSE){
            keep_rows <- NULL
            for (i in summary_List){
                keep_rows <- c(which(i == inpt_Df[,"method"]), keep_rows)
            }
            
            inpt_Df <- rbind(inpt_Df[-summary_Rows,], inpt_Df[intersect(summary_Rows, keep_rows),]) ## put summaries at the bottom of the plot
            summary_Rows <- (nrow(inpt_Df) - length(intersect(summary_Rows, keep_rows)) + 1) : nrow(inpt_Df)
        }
        
        total_Rows <- nrow(inpt_Df)
        
        if(length(lb_Col)  == 0){
            inpt_Df$lb <- inpt_Df[,eff_Col] - qnorm(p = 0.975) * inpt_Df[,se_Col]
            inpt_Df$ub <- inpt_Df[,eff_Col] + qnorm(p = 0.975) * inpt_Df[,se_Col]
            lb_Col <- "lb"
            ub_Col <- "ub"
        }
        
       
        if(onetype_Only == TRUE){
            summary_Rows <- NULL
        }
        
        
        ### 4. Generate a properly-spaced forest plot (without a summary effect) for all studies
        ggforest_oeoo <- function(data_Fm, space_col, eff_col, lb_col, ub_col,se_col, title_text = '', log_ES = FALSE, exp_ES = FALSE, summary_rows, instrument_size) {
            # data_Fm = (name) spaced data frame (must be from output from space_Out()) containing meta-analytic data and appropriate annotation, spacing cols etc
            # space_col = (character) column name giving the column listing spacing of effects in the forest plot
            # eff_col = (character) column name of effect-size column
            # lb_col = (character) column name of the confidence interval lower bound column
            # ub_col = (character) column name of the confidence interval upper bound column
            # title_text = (character) title text for main part of forest plot (will use column name if no text is supplied)
            # log_ES = (logical) convert effect sizes TO natural log scale y/n?
            # exp_ES = (logical) convert effect sizes FROM natural log scale y/n?
            # summary_rows = (integer) vector listing rows with summaries in them
            # instrument_size = (integer) vector listing the number of SNPs to have gone into each instrument
            
            if (nchar(title_text) <= 1) {
                title_text <- eff_col
            }
            ## reassignment of column names to avoid ggplot scope problems
            data_Fm$space_col <- data_Fm[,space_col]
            
            data_Fm$lb_col <- as.numeric(data_Fm[,lb_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,lb_col]))*(exp_ES)
            
            data_Fm$ub_col <- as.numeric(data_Fm[,ub_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,ub_col]))*(exp_ES)
            
            data_Fm$eff_col <- as.numeric(data_Fm[,eff_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,eff_col]))*(exp_ES)
            data_Fm$se_col <- as.numeric(data_Fm[,se_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,se_col]))*(exp_ES)
            
            data_Fm$nsnps <- data_Fm[,instrument_size]
            
            data_Fm$is_summ <- 0
            if(length(summary_rows) > 0 ){ # variable showing whether something is a summary
                
                data_Fm$is_summ <- data_Fm$is_summ + data_Fm$nsnps > 1
            }
            
            data_Fm$is_summ <- factor(data_Fm$is_summ)
            if(log_ES == TRUE){
                stopifnot(all(as.numeric(data_Fm[,lb_col])) >= 0) # throw error if there are negative limits
                
                data_Fm$lb_col <-  log(as.numeric(data_Fm[,lb_col]))
                
                data_Fm$ub_col <-  log(as.numeric(data_Fm[,ub_col]))
                
                data_Fm$eff_col <- log(as.numeric(data_Fm[,eff_col]))
                
                data_Fm$se_col <- log(as.numeric(data_Fm[,se_col]))
            }
            
            
            # ggplot code to generate the forest plot using geom_segments and geom_points and to make a relatively minimal theme
            raw_forest  <- ggplot(data = data_Fm, aes( y = space_col, yend = space_col, x = as.numeric(lb_col), xend = as.numeric(ub_col) )) + geom_segment(aes(colour = is_summ))  + geom_point(aes( y = space_col,  x = as.numeric(eff_col), size = (se_col)^(-2), shape = is_summ, fill = is_summ)) + scale_shape_manual(values = c(15,18)) + scale_size_continuous(range = c(1.25, 8.5)) + scale_colour_manual(values = c('gray40', 'black'))  + theme_bw() + theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),axis.line=element_line(), axis.line.y=element_blank(), title = element_text(size = 30), legend.position = 'none', axis.line.x = element_line(size = 1)) + expand_limits(y = c(data_Fm[,space_col] - 1, data_Fm[,space_col] + 2),x = c(min(as.numeric(lb_col),na.rm = TRUE), max(as.numeric(ub_col), na.rm = TRUE))) + labs(title = title_text) # returns ggplot2 object with the (un-annotated) forest plot
            if(exp_ES == TRUE){
                raw_forest <- raw_forest  + scale_x_continuous(trans = log2_trans()) + geom_vline(xintercept = 1)
            } else {
                raw_forest <- raw_forest + geom_vline(xintercept = 0)
            }
            return(raw_forest)
        }
        ## order and structure effect sizes and CIs for forest plot
        space1 <- spacer(data_Fm = inpt_Df, major = exposure_Name, sort_var = eff_Col,sort_at_all=sort_at_All)
        
        expand_data  <- space_Out(data_Fm = inpt_Df, space_Fm = space1)
        ## Make the forest plot
        if (onetype_Only == TRUE) {
            fo1  <- makeforest( data_Fm = expand_data,  eff_col = eff_Col, lb_col = lb_Col, ub_col = ub_Col, log_ES = log_ES, title_text = forest_Title, exp_ES = exp_ES, se_col = se_Col )
            } else {
            fo1  <- ggforest_oeoo( data_Fm = expand_data, space_col = 'spacing_vec', eff_col = eff_Col, lb_col = lb_Col, ub_col = ub_Col, log_ES = log_ES, title_text = forest_Title, exp_ES = exp_ES,summary_rows = summary_Rows,instrument_size = instrument_Size, se_col = se_Col )
        }
        ## Construct left-hand-side annotations
        left <- anot_side( data_Fm = expand_data, col_names = left_Col_Names, title_list = left_Col_Titles, notitles = FALSE, content_Name = exposure_Name)
        
        ## Construct right-hand-side annotations
        right <- NULL
        if (length(right_Col_Names) > 0){
            right <- anot_side(data_Fm = expand_data, col_names = right_Col_Names, title_list = right_Col_Titles, notitles = FALSE, left = FALSE )
        }
         
        options(warn=0)
        output_forestplot(forest = fo1, left =  left, right = right, outfile_Name = outfile_Name, nrows = nrow(expand_data), ncols = length(right_Col_Names) + length(left_Col_Names), plotType = ifelse(onetype_Only,"oEoOo","oEoO"))
        return(list(forest = fo1, left = left, right = right))
        
    }

