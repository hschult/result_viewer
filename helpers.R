#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#Functions


create_complexheatmap=function(m, mode="raw", unitlabel='auto', rowlabel=T, collabel=T, clustering='none', clustdist='auto', clustmethod='auto', distribution='auto', color_vector_onesided=NULL, color_vector_twosided=NULL, reverse_coloring=FALSE, optimize=T)
{
  ### subset of parameters available for complexheatmap:
  #name="legend"					#title of legend
  #cluster_columns=TRUE				#skip clustering for columns = no reordering of columns
  #cluster_rows=TRUE				#skip clustering for rows = no reordering of rows
  #show_column_dend=TRUE				#show column dendrogram
  #show_row_dend=TRUE				#show row dendrogram
  #column_dend_height=unit(2,"cm")		#height of dendrogram
  #clustering_distance_rows="euclidean"		#euclidean, pearson, spearman, kendall, maximum, manhattan, canberra, binary, minkowski
  #clustering_method_rows="complete"		#ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
  #show_column_dend=T				#show dendrogram
  #show_row_dend=T				#show dendrogram
  #show_row_names=TRUE				#show row names
  #row_names_gp=gpar(fontsize=20))		#increase row names font size
  #km=2						#apply kmeans clustering on rows
  #split=...					#manually split rows into cluster
  #row_names_max_width=unit(3,"cm")		#size for row names
  #column_names_max_height=unit(3,"cm")		#size for column names
  #column_title					#headline for each heatmap
  #width						#width for each heatmap
  #na_col="grey"					#color for na values
  #column_names_gp = gpar(fontsize=8)		#set column label font size
  #row_names_gp = gpar(fontsize=8)		#set row label font size
  
  ### transform matrix according to mode (log2, zscore), set special distance function for mode=zscore, create labels
  if (mode == "zscore"){
    m <- t(scale(t(m)))
    
  }else if (mode == "log2"){
    m <- m + 1
    m <- log2(m)
  }
  
  #replace inf->NA->0
  is.na(m) <- sapply(m, is.infinite)
  m[is.na(m)] <- 0
  
  ### mapping of colors to min/max value depending on distribution type (one- or two-sided)
  min=min(m, na.rm = TRUE)
  absmin=abs(min(m, na.rm = TRUE))
  max=max(m, na.rm = TRUE)
  absmax=abs(max(m, na.rm = TRUE))
  totalmax=max(absmax, absmin)
  message(paste("value min: ", min, sep="\t"))
  message(paste("value max: ", max, sep="\t"))
  if (distribution=='auto') {
    if (min<0 && max>0) {						#two-sided distribution: eg. zscore, log2fc
      distribution='two-sided'
    } else {							#one-sided distribution: eg. count, log2 count
      distribution='one-sided'
    }
  }
  if (distribution=='two-sided') {
    color_vector=heat_color(color_vector_twosided)
  } else {
    color_vector=heat_color(color_vector_onesided)
  }

  if (reverse_coloring == TRUE){                  #reverse colorset
    color_vector = reverse_coloring(color_vector)
  }
  
  message("distribution: ", distribution)
  color_vector_n=length(color_vector)					#number of colors in vector
  maxlimit=max
  minlimit=min
  if (optimize==T && distribution=='two-sided'){				#try to create better color breaks (identical min/max for two-sided, reasonable rounding)
    maxlimit=max(absmax, absmin)
    minlimit= -maxlimit
  }
  
  #catch error minlimit == maxlimit -> no break-vlaues
  if(minlimit == maxlimit){
    maxlimit <- maxlimit + 1
  }
  
  breaks=seq(minlimit,maxlimit,length=color_vector_n)#one break point for each color; even distribution between min/max
  message(paste("breaks min: ", minlimit, sep="\t"))
  message(paste("breaks max: ", maxlimit, sep="\t"))
  col_fun=colorRamp2(breaks,color_vector)					#create color mapping function to fix limits/breaks
  
  ### clustering
  if (clustering=='none') {
    cluster_rows=F
    cluster_columns=F
  } else if (clustering=='row') {
    cluster_rows=T
    cluster_columns=F
  } else if (clustering=='column') {
    cluster_rows=F
    cluster_columns=T
  } else if (clustering=='both') {
    cluster_rows=T
    cluster_columns=T
  }
  
  ### plot
  ht1 = Heatmap(m,
                name=unitlabel,
                column_title=NULL,
                col=col_fun,
                cluster_rows=cluster_rows,
                cluster_columns=cluster_columns,
                clustering_distance_rows=clustdist,
                clustering_distance_columns=clustdist,
                clustering_method_rows=clustmethod,
                clustering_method_columns=clustmethod,
                show_row_names=rowlabel,
                show_column_names=collabel,
                row_names_side="right",
                row_dend_side="left",
                row_dend_width= unit(1,"inches"),
                column_dend_height= unit(1,"inches"),
                row_names_max_width= unit(8,"inches"),
                column_names_max_height= unit(4,"inches"),
                row_names_gp=gpar(fontsize=12),
                column_names_gp=gpar(fontsize=12),
                #show_heatmap_legend = F,
                #raster_device = "png",
                #width=unit({
                #  0.5*ncol(m)
                #}, "cm"),
                #rect_gp = gpar(lineheight = 50, lwd = NA),
                #cell_fun = function(j, i, x, y, width, height, fill){
                #  grid.rect(x = x, y = y, width = width, height = height*0.5, gp = gpar(fill = fill, col = NA))
                #  }, 
                heatmap_legend_param=list(
                  color_bar="continuous", 			#continuous, discrete
                  legend_direction="vertical"			#horizontal, vertical
                )
  )
  
  return(ht1)
}


# heat_color --------------------------------------------------------------

heat_color <- function(palette){
  #get color palette for heatmap
  
  switch(palette,
         #one-sided
        "reds" = reds,
        "viridis" = viridis,
        "plasma" = plasma,
        "inferno" = inferno,
        "magma" = magma,
        "blues" = blues,
        "heat" = heat,
        "cubehelix" = cubehelix,
        "ylgnbu" = ylgnbu,
        #two-sided
        "buwtrd" = buwtrd,
        "rdblgr" = rdblgr,
        "ylwtpu" = ylwtpu,
        "spectral" = spectral
       )
}


# reverse_coloring --------------------------------------------------------

reverse_coloring <- function(colors){
  rev(colors)
}


# heatmap_size ------------------------------------------------------------

heatmap_size <- function(data, row_label = T, column_label = T, clustering){
  col_names_maxlength_label_width=max(sapply(colnames(data),strwidth, units="in", font=12))	#longest column label when plotted in inches	
  col_names_maxlength_label_height=max(sapply(colnames(data),strheight, units="in", font=12))	#highest column label when plotted in inches	
  row_names_maxlength_label_width=max(sapply(rownames(data),strwidth, units="in", font=12))	#longest row label when plotted in inches	
  row_names_maxlength_label_height=max(sapply(rownames(data),strheight, units="in", font=12))	#highest row label when plotted in inches
  col_count <- ncol(data)
  row_count <- nrow(data)
  
  #message("colWidth: ", col_names_maxlength_label_width)
  #message("colheight: ", col_names_maxlength_label_height)
  
  #unit: inches
  width <- 0
  height <- 0.15
  
  #legend
  width <- width + 1
  
  #labels
  if(row_label == TRUE){
    width <- width + row_names_maxlength_label_width + 0.09
  }
  if(column_label == TRUE){
    height <- height + col_names_maxlength_label_width
  }
  
  #clustering c("none", "row", "column", "both")
  if(clustering == "row" && row_count > 1){
    width <- width + 2
  }else if(clustering == "column" && col_count > 1){
    height <- height + 1
  }else if(clustering == "both"){
    if(row_count > 1){width <- width + 2}
    if(col_count > 1){height <- height + 1}
  }
  
  #entries
  
  width <- width + col_count * (col_names_maxlength_label_height + 0.06)
  height <- height + row_count * (row_names_maxlength_label_height + 0.08)
  
  #message("width: ", width)
  #message("height: ", height)
  return(c(width, height))
}

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#Functions

dynamic_matrixsplit <- function(data, reps, plot_type,facet_target,color_palette,facet_cols,widthfactor,heightfactor){
  if (missing(plot_type)){plot_type="box"}
  if (missing(facet_target)){facet_target="gene"}
  if (missing(color_palette)){color_palette = "Dark2"}
  if (missing(facet_cols)){facet_cols=NULL}
  if (missing(widthfactor)){widthfactor=1}
  if (missing(heightfactor)){heightfactor=2}
  
  options(bitmapType='cairo')
  options(width=10000)								#do not wrap lines after 80 chars
  options(max.print=1000000)							#do not wrap lines after 80 chars
  
  genes=nrow(data)												#number of genes (rows in matrix)
  genes_order=unique(as.character(rownames(data)))								#get genes in correct order
  conditions=length(unique(reps$V1))											#number of conditions (columns in matrix)
  conditions_order=unique(as.character(reps$V1))									#get conditions in correct order
  
  ###################
  # Combine and transform dataframes
  ###################
  #detach ids from data
  data_id <- data[[1]]
  data <- data[, sapply(data, is.numeric), with = FALSE]
  
  data_cols <- colnames(data)
  data = transpose(data) 								#switch columns <> rows
  data <- data +1									#add pseudocount
  data = log2(data)								#log transform
  #data <- merge(data_id, data, by=0, sort=F, all= T)
  
  #place former colnames in first column
  data$cols <- data_cols
  setcolorder(data, c("cols", colnames(data)[1:ncol(data)-1]))
  #reattach ids as colnames
  colnames(data)[2:ncol(data)] <- data_id
  setkey(data, cols)

  
  colnames(reps)[1]=c("condition") #add header for condition
  data <- data[reps]					#merge dataframes by rownames
  colnames(data)[1]="sample"							#change Row.names to sample
  data$sample=NULL								#completely remove sample column again
  data=transform(data,condition=factor(condition, levels=unique(condition)))	#order conditions in plot according to reps.txt (instead of alphabetic)
  
  data=melt(data, id.vars = "condition")
  
  ###################
  # Choose color palette
  ###################
  
  
  if (facet_target=="gene") {											#facet = gene
    num_colors=conditions
  }
  if (facet_target=="condition") {										#facet = condition
    num_colors=genes
  }

  
  if (color_palette=="None") {
    color_fill_grayscale="grey75"										#color to use for filling geoms in grayscale mode
    colour_palette=rep(color_fill_grayscale,num_colors)
  } else {
    options(warn=-1)
    colour_palette=colorRampPalette(brewer.pal(12, color_palette))(num_colors)
    options(warn=0)
  }
  
  ###################
  # Function to get standard error for error bars (box, bar, violin)
  ###################
  
  get.se = function(y) {
    se=sd(y)/sqrt(length(y))
    mu=mean(y)
    c(ymin=mu-se, ymax=mu+se)
  }
  
  ###################
  # Function to collapse the dataframe to the mean and the standard deviation/error before plotting (ONLY used for line plot)
  ###################
  
  # data : a data frame
  # varname : the name of a column containing the variable to be summarized
  # groupnames : vector of column names to be used as grouping variables
  data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
      c(
        mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE),
        se = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]]))
      )
    }
    data_sum <- ddply(data, groupnames, .fun=summary_func, varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
  }
  if (plot_type=="line") {
    data = data_summary(data, varname="value", groupnames=c("condition","variable"))			#collapse the dataframe to the mean and the standard deviation for line plot
    #head(data)
  }
  
  if (plot_type=="box" || plot_type=="violin" || plot_type=="bar" || plot_type=="line") {
    
    ###################
    # Calculate number of facets / row, pdf size for bar, box, violin plots
    ###################
    
    if (facet_target=="gene") {							#facet = gene
      facets=genes								#eg. number of genes
    }
    if (facet_target=="condition") {						#facet = condition
      facets=conditions							#eg. number of conditions
    }
    #cat("facets: ", facets, "\n")
    
    ##calculate number of columns/rows for facetting the plot
    if ( is.null(facet_cols) ) { 						#no fixed column number for facet grid -> rectangular
      min_cols=5								#rectangular; only start a second row if more than this number of columns
      if (facets <= min_cols) {						#grow the first row to 5 columns before starting a second
        facet_cols=facets
        facet_rows=1
      } else {								#auto width/height = 2:1
        facet_cols = ceiling(sqrt(facets))*1.5 					#quadratic
        facet_rows = ceiling(facets/facet_cols)
        #			facet_rows = facet_cols
        #facet_cols = ceiling((facets*(1/5)))
        #facet_rows = ceiling(facets/facet_cols)
      }
    } else {									#fixed number of facet columns				
      if (facet_cols>facets) {
        facet_cols=facets
      }
      facet_rows = ceiling(facets/facet_cols)
    }
    #cat("facet_cols: ", facet_cols, "\n")
    #cat("facet_rows: ", facet_rows, "\n")
    
    #if (facet_target=="gene") {							#facet = gene
    #  pdf_width = ((facet_cols*conditions*0.5)+1)*widthfactor
    #  pdf_height = ((facet_rows*3)+0.1)*heightfactor
    #}
    #if (facet_target=="condition") {						#facet = condition
    #  pdf_width = ((facet_cols*genes*0.4)+1)*widthfactor
    #  pdf_height = ((facet_rows*3.5)+0.1)*heightfactor
    #}
    
    ###################
    # Set common parameters for all plots
    ###################
    
    theme1 <- theme (															#no gray background or helper lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x=element_text(angle=45, hjust=1, vjust=1),										#x-axis sample lables = 45 degrees
      strip.background=element_blank(),
      panel.border = element_rect(colour="black"),
      legend.position="none",														#remove legend
      legend.title=element_blank(),
      axis.title.x = element_blank()
      
      #	axis.line.x = element_line(size=.3),
      #	axis.line.y = element_line(size=.3),
      #	panel.background = element_blank(),
      #	axis.title.y = element_text(face="bold", color="black", size=10),
      #	plot.title = element_text(face="bold", color="black", size=12),
      #	axis.text.x=element_text(angle=90, hjust=1)											#x-axis sample lables = vertical
    )
    
    matrixplot = ggplot(data, aes(y=value))
    
    matrixplot = matrixplot + 
      theme_bw() + theme1 +
      ylab("Log2 Expression") +
      xlab("") +
      scale_fill_manual(values=colour_palette) +
      scale_color_manual(values=colour_palette)
    
    ###################
    # Handle facetting and special parameters for line plot (no facetting, etc.)
    ###################
    
    if (facet_target=="gene") {														#facet = gene
      matrixplot = matrixplot + aes(x=condition, fill=condition)
      
      if (plot_type=="line") {													#line plot: no facetting, different size algorithm
        matrixplot = matrixplot + aes(x=variable, colour=condition, group=condition, fill=NULL)
        matrixplot = matrixplot + scale_x_discrete(expand=c(0.05,0.05))								#expand to reduce the whitespace inside the plot (left/right)
        #pdf_width = ((genes*1)+1) * widthfactor
        #pdf_height = ((conditions * 0.15) +2) * heightfactor
        #pdf_height = 3
        
      } else {
        ##scales="free_y"		#separate y scale for each facet (similar: free_x)
        ##ncol=3			#fix number of columns
        ##nrow=3			#fix number of rows
        matrixplot = matrixplot + facet_wrap(~variable, ncol=facet_cols, scales="free_x")
      }
    }
    if (facet_target=="condition") {													#facet = condition
      matrixplot = matrixplot + aes(x=variable, fill=variable)
      
      if (plot_type=="line") {													#line plot: no facetting, different size algorithm
        matrixplot = matrixplot + aes(x=condition, colour=variable, group=variable, fill=NULL)
        matrixplot = matrixplot + scale_x_discrete(expand=c(0.05,0.05))								#expand to reduce the whitespace inside the plot (left/right)
        #pdf_width = ((conditions*1)+1) * widthfactor
        #pdf_height = ((genes * 0.15) +2) * heightfactor
        #pdf_height = 3
      } else {
        ##scales="free_y"		#separate y scale for each facet (similar: free_x)
        ##ncol=3			#fix number of columns
        ##nrow=3			#fix number of rows
        matrixplot = matrixplot + facet_wrap(~condition, ncol=facet_cols, scales="free_x")
      }
    }
    
    ###################
    # Further handle plot types
    ###################
    
    if (plot_type=="box") {																#plot type: box
      matrixplot = matrixplot + stat_boxplot(geom='errorbar', width=0.2) 										#add horizontal line for errorbar
      matrixplot = matrixplot + geom_boxplot(position=position_dodge(1))
      #matrixplot = matrixplot + stat_summary(fun.data=get.se, geom="errorbar", width=0.2)					#error bar of standard error
    }
    if (plot_type=="violin") {																#plot type: violin
      matrixplot = matrixplot + geom_violin()
      #matrixplot = matrixplot + stat_summary(fun.y="median", geom="point")										#add median dot
      #matrixplot = matrixplot + stat_summary(fun.data=get.se, geom="errorbar", width=0.2, position=position_dodge())					#error bar of standard error
    }
    if (plot_type=="bar") {																#plot type: box
      matrixplot = matrixplot + stat_summary(fun.y=mean, geom="bar", position="dodge")							#bar plot of the mean (color=condition)
      matrixplot = matrixplot + stat_summary(fun.data=get.se, geom="errorbar", width=0.2, position=position_dodge())					#error bar of standard error
    }
    if (plot_type=="line") {
      matrixplot = matrixplot + theme(legend.position="right")
      #matrixplot = matrixplot + geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=0.05)								#error bar = standard deviation
      matrixplot = matrixplot + geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.05)								#error bar = standard error
      matrixplot = matrixplot + geom_line() + geom_point()											#bar plot of the mean (color=condition)
    }
  }
  
  return (matrixplot)
}



# unique columns ----------------------------------------------------------

uniqueColumns <- function(data){
  allNames <- colnames(data)
  uniqueNames <- make.names(allNames, unique = TRUE)
  names(data) <- uniqueNames
  
  return(data)
}


#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#Variables

#options(bitmapType='cairo')

# set some reasonable defaults for the options that are needed,
# but were not specified.

#if ( is.null(opt$unitlabel) ) { opt$unitlabel='auto'; }
#if ( is.null(opt$clustdist) ) { opt$clustdist='auto'; }
#if ( is.null(opt$clustmethod) ) { opt$clustmethod='auto'; }
#if ( is.null(opt$widthfactor) ) { opt$widthfactor = 1; }
#if ( is.null(opt$heightfactor) ) { opt$heightfactor = 1; }



##########################################################################################################################################################
########### setup color palettes available by cmd interface
##########################################################################################################################################################

num_colors=256

################## one-sided (expression, enrichment)
heat		=colorRampPalette(rev(brewer.pal(9,"YlOrRd")))(num_colors)
viridis		=viridis(num_colors)
magma		=magma(num_colors)
inferno		=inferno(num_colors)
plasma		=plasma(num_colors)
ylgnbu		=colorRampPalette(brewer.pal(9, "YlGnBu"))(num_colors)
blues		=colorRampPalette(brewer.pal(9,"Blues"))(num_colors)
reds		=colorRampPalette(brewer.pal(9,"Reds"))(num_colors)
cubehelix	=cubeHelix(num_colors)

col1=colorRampPalette(c("black", "orange", "yellow"))(num_colors)						#one-sided (0 .. x): go enrichment
col2=colorRampPalette(c("khaki1","yellow","orange","red","darkred"))(num_colors)				#one-sided (0 .. x): expression
col2a=heat.colors(num_colors)											#one-sided (0 .. x): expression 
col2c=colorRampPalette(rev(c("khaki1","yellow","orange","red","darkred")))(num_colors)				#one-sided (0 .. x): expression
col10=colorRampPalette(brewer.pal(9, "GnBu"))(num_colors)							#one-sided (0 .. x): expression
col17=colorRampPalette(brewer.pal(9, "PuBuGn"))(num_colors)							#one-sided (0 .. x): expression
col18=colorRampPalette(c("#000041", "#0000CB", "#0081FF", "#02DA81", "#80FE1A", "#FDEE02", "#FFAB00", "#FF3300"))(num_colors)	#one-sided (0 .. x): expression, ~=spectral
col18a=colorRampPalette(c("#781C81", "#3F4EA1", "#4683C1", "#57A3AD", "#6DB388", "#B1BE4E", "#DFA53A", "#E7742F", "#D92120"))(num_colors)	#one-sided (0 .. x): expression, ~=spectral

################## two-sided (log2fc, zscore)
buwtrd		=colorRampPalette(c("royalblue4","steelblue4","white","indianred","firebrick4"))(num_colors)
rdblgr		=redgreen(num_colors)
rdylgr		=colorRampPalette(brewer.pal(11, "RdYlGn"))(num_colors)
ylwtpu		=colorRampPalette(c("gold","white","white","mediumpurple4"))(num_colors)					#two-sided (-1 .. +1): correlation 
spectral	=colorRampPalette(brewer.pal(11, "Spectral"))(num_colors)

col4=colorRampPalette(c("dodgerblue4", "cadetblue1", "yellow", "darkolivegreen1", "darkgreen"))(num_colors)	#two-sided (-x .. +x): fold-change
col5=colorRampPalette(c("darkslategray","darkturquoise", "cornsilk", "indianred3", "red3"))(num_colors)		#two-sided (-x .. +x): fold-change
col6=colorRampPalette(c("yellow","grey25","red"))(num_colors)							#two-sided (-x .. +x): fold-change
col7a=colorRampPalette(brewer.pal(9, "RdBu"))(num_colors)							#two-sided (-x .. +x): fold-change
col9=colorRampPalette(c("chartreuse3","white","firebrick1"))(num_colors)					#two-sided (-x .. +x): fold-change
col12=colorRampPalette(brewer.pal(11, "RdYlBu"))(num_colors)							#two-sided (-x .. +x): fold-change
col13=colorRampPalette(brewer.pal(11, "RdGy"))(num_colors)							#two-sided (-x .. +x): fold-change
col14=colorRampPalette(brewer.pal(11, "PuOr"))(num_colors)							#two-sided (-x .. +x): fold-change
col15=colorRampPalette(c("royalblue3","steelblue3","white","indianred3","firebrick3"))(num_colors)		#two-sided (-x .. +x): fold-change

################## scatter colors
lightgoldenrod1 <- colorRampPalette(colors = c("lightgoldenrod1", "indianred2", "steelblue2"))(num_colors)
azure2 <- colorRampPalette(colors = c("azure2", "red3", "blue3"))(num_colors)	

#if ( opt$colorpalette=='auto' ) { 							#automatic mode
#	color_vector_onesided=heat
#	color_vector_twosided=buwtrd
#} else {
#	colorpalettemain=sub("^(.*)_r", "\\1", opt$colorpalette)				#remove _r (if it existed)
#	palette=get(colorpalettemain)								#get main palette
#	if (grepl('_r$', opt$colorpalette)) {							#reverse colors of param ended with _r
#		palette=rev(palette)
#	} 
#	color_vector_onesided=palette
#	color_vector_twosided=palette
#}
#}




##########################################################################################################################################################
########### read and transform input matrix
##########################################################################################################################################################

#setwd(".")
#args <- commandArgs(trailingOnly = TRUE)					#parse commandline
#data=read.table(file=opt$matrix,sep="\t",header=T,stringsAsFactors=F,row.names=NULL,check.names=FALSE, quote="")	#with header, column 0 = rownames, do not convert strings to character vectors
#fileprefix = sub("^(.*)[.].*", "\\1", opt$matrix)				#remove suffix from filename (everything after last ".") = output file prefix

#data_temp <- data[!duplicated(data[,1]),]					#remove rows where 1. column is duplicate
#data <- data_temp
#data_temp <- data[,-1]								#convert 1. column to rownames
#rownames(data_temp) <- data[,1]
#data <- data_temp

#data = na.omit(data)   							#remove rows with missing values.
#data[is.na(data)] <- 0								#replace NA with 0	

#if (!is.null(opt$categories)) {
#	categories=read.table(file=opt$categories,sep="\t",header=T,stringsAsFactors=T,row.names=1,check.names=FALSE,na.strings=c("","NA"))	#with header, column 0 = rownames, do not convert strings to character vectors
#}

##########################################################################################################################################################
########### create heatmap
##########################################################################################################################################################


###############
##  complexheatmap wrapper function
###############



###############
##  create heatmap
###############


##########################################################################################################################################################
########### categories = additional colour bars
##########################################################################################################################################################

categories_nr=0


# Scatterplot -------------------------------------------------------------

#data:
#   column 1: id
#   column 2, 3(, 4): x, y(, z) 

create_scatterplot <- function(data, round = F, log10 = F, transparency = 1, pointsize = 2, colors = NULL, maxaxis = NULL, x_label = "", y_label = "", z_label = "", density = T, line = T, categorized = F){
  #get intern columnnames
  x_head <- colnames(data)[2]
  y_head <- colnames(data)[3]
  if(ncol(data) >= 4){
    z_head <- colnames(data)[4]
  }
  
  #set labelnames if needed
  x_label <- ifelse(nchar(x_label), x_label, x_head)
  y_label <- ifelse(nchar(y_label), y_label, y_head)
  if(ncol(data) >= 4){
    z_label <- ifelse(nchar(z_label), z_label, z_head)
  }
  

  #delete rows where both 0 or at least one NA
  rows <- which(as.logical((data[[2]]!=0) + (data[[3]] != 0)))
  data <- data[rows]
  
  #round data to Integer
  if(round == TRUE){ 
    if(categorized == TRUE){
      data[,2:3] <- round(data[,2:3])
    }else{
      data[,2:ncol(data)] <- round(data[,2:ncol(data)])
    }
  }
  #log10
  if(log10 == TRUE){
    if(categorized == TRUE){
      data[,2:3] <- log10(data[,2:3])
    }else{
      data[,2:ncol(data)] <- log10(data[,2:ncol(data)])
    }
  }
  
  #autoscale axis
  if(is.null(maxaxis)){
    maxcolcounts <- apply(data[,2:3], 2, max)
    maxaxis <- ceiling(max(maxcolcounts, na.rm = T))
  }
  
  #replace inf->NA->0
  is.na(data) <- sapply(data, is.infinite)
  data[is.na(data)] <- 0

  ########## assemble plot ##########
  
  theme1 <- theme (											#no gray background or helper lines
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(size=.3),
    axis.line.y = element_line(size=.3),
    axis.title.x = element_text(face="bold", color="black", size=10),
    axis.title.y = element_text(face="bold", color="black", size=10),
    plot.title = element_text(face="bold", color="black", size=12)
    #		legend.background = element_rect(color = "red")			#border color
    #		legend.key = element_rect("green")						#not working!
  )
  
  ###scatter with color axis
  if(ncol(data) >= 4 && categorized == FALSE){
    plot <- ggplot(data = data, aes(x = data[[x_head]],y = data[[y_head]], color = data[[z_head]])) +
      ###color_gradient
      scale_color_gradientn(colors = scatter_color(colors), name = z_label) + 
      ### point options
      geom_point(size=pointsize, alpha=transparency)
  
  ###scatter with categories    
  }else if(ncol(data) >= 4 && categorized == TRUE){
    ###categorized plot
    color_vector <- scatter_color(colors)
    breaks <- seq(from = 1,to = nrow(data), length.out = length(color_vector))
    
    plot <- ggplot(data = data, aes(x = data[[x_head]],y = data[[y_head]])) +
      ### point options
      geom_point(size=pointsize, alpha=transparency, aes(color = factor(data[[z_head]]))) +
    
      scale_color_manual (
        #labels = data[,z_head],
        values = colorRampPalette(color_vector)(length(unique(data[[z_head]]))), #get color for each value,
        #breaks = ,
        drop=FALSE,								#to avoid dropping empty factors
        name = z_label
        #			guide=guide_legend(title="sdsds" )					#legend for points
      )
      
  }else{
    plot <- ggplot(data = data, aes(x = data[[x_head]],y = data[[y_head]])) +
      geom_point(size=pointsize, alpha=transparency)
  }
  
  plot <- plot +
    theme1 +
    
    ### binhex
    #		stat_binhex(bins=30) + 								
     
    
    ### smooth curve
    #		geom_smooth(method="loess", se=FALSE, color="black") +				#se=display confidence interval (shaded area)
    #		geom_smooth(method="loess", se=FALSE) +				#se=display confidence interval (shaded area)
    
    ### additional density plot at x and y axis
    #		geom_rug(col="darkred", alpha=.1) +						#density plot at x and y axis
    
    ### axis range and labels
    xlim(0, maxaxis) +								#set x axis limits
    ylim(0, maxaxis) +								#set y axis limits
    xlab(x_label) +								#axis labels
    ylab(y_label) 
  
  #		guides(fill =guide_legend(keywidth=3, keyheight=1))				#legend for density
  
  if(line == TRUE){
    ### diagonal line
    plot$layers <- c(geom_abline(intercept=0, slope=1), plot$layers) #plot$layers, so line is in background
  }  
  
  if(density == TRUE){
    ### kernel density
    #stat_density2d(geom="tile", aes(fill=..density..), n=200, contour=FALSE) +		#n=resolution; density more sparse
    plot$layers <- c(stat_density2d(geom="tile", aes(fill=..density..^0.25), n=200, contour=FALSE), plot$layers)#n=resolution; density less sparse 
    
    plot <- plot + scale_fill_gradient(low="white", high="black") +
    #guides(fill=FALSE) +		#remove density legend
    labs(fill="Density")
  }
  
  return(plot)
}

# scatter_color --------------------------------------------------------------

scatter_color <- function(palette){
  #get color palette for scatter
  
  colors <- switch(palette,
         "lightgoldenrod1" = lightgoldenrod1,
         "azure2" = azure2
         )
}
