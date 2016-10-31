### Functions to help with the generation of forest plots; no logic to actually generate forest plots 
### Version: 0.60
### Date: 20160818
### Authors: cl14022, ph14916
### Latest change: Split into separate files

### This script does:




###### BEGIN Shared Helper Functions
###### 

spacer <- function(major, sort_var, data_Fm) {
    
    ## This function fills in a 'spacing vector', to determine placement of rows in a forest plot, grouped by the major variable (i.e. exposure) and sorted by the sort variable (i.e. effect size)
    ## and also to ensure that subgroups have bold face headers, a content list, which gives the content for headers, and a row list, of the rows to be filled in with this header data
    
    # major = (character), column name of the column containing the name of the column most important for ordering rows, e.g. exposure (SNP for it)
    # Data_Fm = (name) name of the data frame containing the data of interest
    N_entries <- nrow(data_Fm)  ## includes summary data
    
    # The layout of rows is one per entry, with one per grouping and an overall buffer at to
    N_space <- N_entries + length(unique(data_Fm[, major])) ## 
    
    spacing_vec <- vector(mode = "numeric")
    content_list <- list()
    attr_vec <- vector(mode = "character")
    row_vec <- vector(mode = "numeric")

    idx <- 1 ## have to do it in this shameful way
    for (entry in unique(data_Fm[order(data_Fm[,sort_var], decreasing = TRUE),major])) {
        
        # loop over entry ids to fill the 'spacing vector,' which determines placement of the rows in the forest plot
        entry_idxs <- which(data_Fm[,major] == entry)

        n_entries <- length(entry_idxs)
        spacing_vec[(idx :(idx + n_entries + 1))] <- idx : (idx + n_entries + 1)
        
        # In the same loop, fill a list with entry names but also with headers to serve as the primary annotation of the forest plot
        content_list[[idx]] <- entry ## The header is the value of the major sort
        content_list[((idx + 1):(idx + n_entries + 1))] <- "" ## spacer
        
        # A list of attributes for the typeface so that headers are in bold face
        attr_vec[[idx]] <- "bold"
        attr_vec[((idx + 1):(idx + n_entries + 1))] <- "plain"
        
        # The most important list, contains mapping between spaces on the forest plot and rows in the meta-analytic data frame, sorted
        row_vec[[idx]] <- NA
        row_vec[[(idx + n_entries + 1)]] <- NA
    
        if( n_entries > 1 ){
          row_vec[((idx + 1):(idx + n_entries))] <- entry_idxs[order(data_Fm[entry_idxs,sort_var], decreasing = TRUE)]
        } else {
          row_vec[((idx + 1):(idx + n_entries))] <- entry_idxs
          }
        
        idx <- idx + n_entries + 2 ## add 2 so that the next entry begins in the correct place
    }
   ##  browser()   
    ## repeat for summary_rows
    # returns data frame giving spacings, primary annotation (content_list), typeface attributes for the primary annotation (attr_vec), and  mapping between rows in forest plot and rows in meta-analytic data frame
    result <- data.frame( spacing_vec = (1 + length(spacing_vec) - spacing_vec), content_list = unlist(content_list), row_vec = row_vec, attr_vec = attr_vec )
    return( result )
}

space_Out <- function(data_Fm, space_Fm) {
    # a function to expand the initial data frame and add the spacing vector in a proper way
    # data_Fm = (name) data frame containing the meta-analytic data
    # space_Fm = (name) data frame containing the spacing vector and the mapping between spaces in the forest plot and rows in the meta-analytic data
    exp_data_Fm  <- data.frame(space_Fm, stringsAsFactors = FALSE)
    for (jj in colnames(data_Fm)) {
        exp_data_Fm[,jj] <- data_Fm[space_Fm$row_vec,jj]
    }
    # Returns a data frame containing spacing, primary annotations, typeface attributes for the primary annotation, and correctly-mapped rows in the meta-analytic data frame
    return(exp_data_Fm)
}


### 5. Generate a properly-spaced set of annotations for the left-hand side (lhs) of the figure
### and
### 6. Generate a properly-spaced set of annotations for the right-hand side (rhs) of the figure
## 5.1. Function to generate a single column of annotations
anot_col <- function(data_Fm, text_col, title_text = '', notitles) {
    # data_Fm = (name) space_Out data frame
    # text_col = (character) name of column containing annotations
    # title_text = (character) title of annotation column, defaults to the column name if nothing is entered
    
    data_Fm[is.na(data_Fm[, text_col]),text_col] <- ""
    data_Fm$labs <-  format(data_Fm[,text_col],digits = 3, width = 16)

    
    # A hard rule to set the width of the annotation column, which sometimes truncates very wide columns (complex disease names, numbers with 16 digits, etc)
    text_widths <- c(-1, max(10,0.5 * max(sapply( as.character(data_Fm[,text_col]),nchar ))))

    # GGplot rendering of the annotation column
    coltext  <- ggplot(data = data_Fm, aes( y = spacing_vec, x = 0, label = labs, fontface = attr_vec )) + geom_text(hjust = "inward") + theme_bw() + theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(),axis.text.x = element_text(colour = "white"),axis.ticks.x = element_line(colour = "white"), axis.title = element_blank(), rect = element_blank(), panel.grid = element_blank(), title = element_text(size = 19) ) + expand_limits(x = text_widths,  y = c(data_Fm[,'spacing_vec'] - 1, data_Fm[,'spacing_vec'] + 2)) + labs(title = title_text, size = 18)  # returns two-item list with left_text, the GGplot annotations, and text_widths, the x-axis limits of the plot
    if(notitles == TRUE){
        coltext <- coltext + theme(title = element_blank())
    }
    return(list(col_text = coltext, text_widths = text_widths))
}

## 5.2. Function to aggregate all of the left-hand (all of the right-hand) notation columns into a single subplot, which works by applying the anot_col function to each item in the list of annotation columns
anot_side <- function(data_Fm, col_names, title_list = '', notitles = FALSE, content_Name = "content_list", left = TRUE) {
    library(gtable)
    library(grid)
    # data_Fm = (name) space_Out data frame containing effect sizes, all annotations, and spacings
    # space_col = (character) name of the column containing the spacing vector
    # col_names = (character vector) vector of the column names for the columns to be the left-hand-side or right-hand-side annotations
    # title_list = (character vector) vector of the titles to be given to each of the annotation columns, defaults to the column names
   
  
    if(length(title_list) < 1){
      title_list <- col_names
    }
  
    if(left == TRUE){
      col_names <- c("content_list", col_names)
      title_list <- c(content_Name, title_list)
      } 

    relative_widths <- vector(mode = "numeric", length = length(col_names))
    output <- vector(mode = "list")


    for (i in 1:length(col_names)) {
        # loop to get the widths of each annotation column and to group the annotation objects together
        col <- anot_col( data_Fm = data_Fm, text_col = col_names[i], title_text = title_list[[i]], notitles = notitles)
        relative_widths[i] <- col$text_widths[2] - col$text_widths[1]
        output[[col_names[i]]] <- ggplotGrob(col$col_text)
    }
    # calculate the widths of the annotation columns, relative to each other, and then multiply it by a third so that it fits onto one side of the plot

    output$relative_widths <- relative_widths / (sum(relative_widths)) * 0.29
    # returns an output object containing: each annotation plot, listed under its own name, and the relative_widths of the columns
    return(output)
}


### 4. Generate a properly-spaced forest plot (without a summary effect) for all studies
makeforest <- function(data_Fm, eff_col, lb_col, ub_col, se_col, title_text = '', log_ES = FALSE, exp_ES = FALSE, se_Weight = FALSE ) {
    # data_Fm = (name) spaced data frame (must be from output from space_Out()) containing meta-analytic data and appropriate annotation, spacing cols etc
    # space_col = (character) column name giving the column listing spacing of effects in the forest plot
    # eff_col = (character) column name of effect-size column
    # lb_col = (character) column name of the confidence interval lower bound column
    # ub_col = (character) column name of the confidence interval upper bound column
    # title_text = (character) title text for main part of forest plot (will use column name if no text is supplied)
    # log_ES = (logical) convert effect sizes TO natural log scale y/n?
    # exp_ES = (logical) convert effect sizes FROM natural log scale y/n?
    # summary_rows = (integer) vector listing rows with summaries in them
    
    
    if (nchar(title_text) <= 1) {
        title_text <- eff_col
    }
    ## reassignment of column names to avoid ggplot scope problems

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
    
    
    # ggplot code to generate the forest plot using geom_segments and geom_points and to make a relatively minimal theme
    raw_forest  <- ggplot(data = data_Fm, aes( y = spacing_vec, yend = spacing_vec, x = as.numeric(lb_col), xend = as.numeric(ub_col) )) + geom_segment()  + geom_point(aes( y = spacing_vec,  x = as.numeric(eff_col), size = (se_col)^(-2 * as.numeric(se_Weight))), shape = 15) + scale_size_continuous(range = c(1.25, 8.5))  + theme_bw() + theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),axis.line=element_line(), axis.line.y=element_blank(), title = element_text(size = 20), legend.position = 'none', axis.line.x = element_line(size = 1)) + expand_limits(y = c(data_Fm[,'spacing_vec'] - 1, data_Fm[,'spacing_vec'] + 2),x = c(min(as.numeric(lb_col),na.rm = TRUE), max(as.numeric(ub_col), na.rm = TRUE))) + labs(title = title_text) # returns ggplot2 object with the (un-annotated) forest plot
                                                                                                                                                                        
    if (exp_ES == TRUE) {
        raw_forest <- raw_forest  + scale_x_continuous(trans = log2_trans(), breaks = scales::pretty_breaks(n = 5)) + geom_vline(xintercept = 1) 
    } else {
        raw_forest <- raw_forest + geom_vline(xintercept = 0) + scale_x_continuous(breaks = scales::pretty_breaks(n = 5))
    }
    
    return(raw_forest)
}

###### END Shared Helper Functions

##### BEGIN output_forestplot()
output_forestplot <- function(forest, left, right, outfile_Name, nrows, ncols, plotType){
    library(gtable)
    library(grid)
    group_top_Plots <- function(forst_Pt, left_Hs, right_Hs = NULL) {
        # forst_Pt = (name) the name of the forest plot object
        # left_Hs = (name) the name of the aggregate annotation object for the left hand side of the plot
        # right_Hs  = (name) the name of the aggregate annotation object for the right hand side of the plot

        # Aggregate all of the plots into a single grid object for plotting
        grob_Bag <- vector(mode = 'list')

        left_RW <- left_Hs$relative_widths
        left_Grobs <- left_Hs
        left_Grobs$relative_widths <- NULL

        for (i in 1:length(left_Grobs)) {
            grob_Bag[paste('l',names(left_Grobs)[i],sep = '')] <- left_Grobs[i]
        }

        grob_Bag$m_forest <- ggplotGrob(forst_Pt)
        right_RW <- NULL
        if(!is.null(right_Hs)){
            right_RW <- right_Hs$relative_widths
            right_Grobs <- right_Hs
            right_Grobs$relative_widths <- NULL
            for (i in 1:length(right_Grobs)) {
                grob_Bag[paste('r',names(right_Grobs)[i], sep = '')] <- right_Grobs[i]
            }
        }


        width_vec <- c(left_RW,0.42,right_RW)
        width_vec <- width_vec / sum(width_vec)
        # convert the grid objects (now grouped) into a table of grid objects that can be plotted using grid.draw
        grp_FP <- gtable_matrix( name = "groupplot", grobs = matrix(grob_Bag, nrow = 1), widths = unit(width_vec, "npc"), heights = unit(1,"npc") )

        # return the grid object table, to be plotted
        return(grp_FP)
    }
    group <- group_top_Plots(forst_Pt = forest, left_Hs = left, right_Hs = right)
    ## draw and export the annotated forest plot
    grid.newpage()
    ### Problem occurs here, which means it's a problem in defining group
    grid.draw(group)
    ## browser()

    typeSpec <- data.frame(c('oEmO', 'mEoO', 'oEoO', 'oEoOo'),1:4)
    colnames(typeSpec) <- c('type', 'spec')
    widthSpec <- switch(typeSpec$spec[match(plotType,typeSpec$type)], 6500 + 480*(ncols%%30), 6000 + 480*(ncols%%30), 6500 + 480*(ncols%%30), 6500 + 480*(ncols%%30))
    heightSpec <- switch(typeSpec$spec[match(plotType,typeSpec$type)],3000+ 250*(nrows%%13), 3000+ 250*(nrows%%13), 6000+ 500*(nrows%%13), 3000+ 125*(nrows%%13))
    if(length(outfile_Name) > 0 ){
        png(file = outfile_Name,width = widthSpec, height = heightSpec, res = 300)
        grid.newpage()
        grid.draw(group)
        dev.off()
    }
    return(group)
}
##### END output_forestplot()
