#!/usr/bin/env Rscript
#
#create singular matrix for counts, log2fc, enrichment etc.
#remove rows with duplicate names
#further categorization / row can be optionally included (will be shown right of the heatmap)
#
#format: heatmap matrix (can be counts, log2c, pvalue, zscore, ...)
#name	zf_6h	zf_1d	md_6h	md_2d
#gene1	1	21	43	21
#gene2	4	53	33	24
#...
#
#format: categories matrix (for optional categorization shown right of heatmap)
#name	Pathway	Methylation
#gene1	immune response 	up
#gene2	immune response 	down
#...
#
#default modes:
#
#raw:
#data transform:	raw
#clustdist:		euclidean
#clustmethod:		average
#
#expression:
#data transorm:		+1, log2
#clustdist:		euclidean
#clustmethod:		average
#
#zscore:
#data transorm:		row zscore
#clustdist:		pearson
#clustmethod:		average
#
#known issues:
#clustering crashes for very large heatmaps (try -c none)
#
#needs r 3.2.2-local
#syntax:
#. switch_modules.sh -m r -n 3.2.2-local				#load this module version if not already loaded
#./complexheatmap_single.r -m matrix.txt
#. switch_modules.sh -m r						#switch back to old module version
###############################################################################################################

require('getopt');
options(bitmapType='cairo')
options(width=5000)

##0=no, 1=required, 2=optional
##logical,integer,double,complex,character

options = matrix(c(
  'matrix','m',1,'character','Tab-delimited input matrix bearing id and count columns (including headline).',
  'mode','z',2,'character','raw (show data as it is), expression (do log2 first), zscore (perform row-wise zscore). [default=expression]',
  'colorpalette','p',2,'character','Color palettes (add \"_r\" to the palette to reverse the order of colors; eg. heat_r).
  					one-sided distributions (* = default): 
						*heat		red, orange, yellow
						viridis		blue, green, yellow
						magma		black, violet, orange, yellow
						inferno		black, violet, orange, yellow (more saturated than magma)
						plasma		violet, orange, yellow
						ylgnbu		yellow, green, blue
						blues		white, blue
						reds		white, red
						cubehelix	black, green, pink, white
						...
  					two-sided distributions (* = default): 
						*buwtrd 	blue, white, red
						rdblgr		red, black, green
						rdylgr		red, yellow, green
						ylwtpu		yellow, white, purple
						spectral	rainbow
						...',
  'distribution','d',2,'character','type of data distribution. auto (determine automatically), one (one-sided, eg. expression, enrichment), two (two-sided, eg. log2fc, zscore.). [default=auto]',
  'clustering','c',2,'character','clustering and dendrogram: none, row, column, both. [default: row]',
  'clustdist','l',2,'character','clustering distance function: euclidean, pearson, spearman, kendall, maximum, manhattan, canberra, binary, minkowski. [default: determine automatically]',
  'clustmethod','n',2,'character','average, ward.D, ward.D2, single, complete, mcquitty, median, centroid. [default: average]',
  'norowlabel','r',2,'logical','Do not show row labels. This will automatically resize the pdf to roughly fit on one page. [default: show row labels]',
  'nocollabel','s',2,'logical','Do not show column labels. [default: show column labels]',
  'unitlabel','u',2,'character','Label for unit displayed above color key.',
  'title','v',2,'character','Title for plot. [default: input filename w/o suffix].',
  'subtitle','k',2,'character','Subtitle for plot. [default: mode, clustdist, clustmethod].',
  'categories','t',2,'character','Tab-delimited input matrix bearing id and category columns (including headline). [default: not included]',
  'category_colours','i',2,'character','How to select colors for the categories: preset, random, *Rcolorbrewer Set* [default: preset].
  				preset: pick colors from  pre-defined palettes: one for each category (up to 5 different categories).
				random: pick random color.
				*Rcolorbrewer Set*: eg. Dark2, Accent, Set1: Use the same palette for all categories.',
  'widthfactor','w',2,'double','Factor to multiply the width with. [default: 1].',
  'heightfactor','e',2,'double','Factor to multiply the height with. [default: 1].',
  'outfile','o',2,'character','Output file. [default: input_prefix.pdf].',
  'help','h',0,'logical','Provides command line help.'),ncol=5,byrow=T)

opt = getopt(options);

# help was asked for.
if ( !is.null(opt$help) ) {
  cat(getopt(options, usage=TRUE),file = stderr());
  q(status=1);
}

if ( is.null(opt$matrix) ) {
  cat(getopt(options, usage=TRUE),file = stderr());
  q(status=2);
}

# set some reasonable defaults for the options that are needed,
# but were not specified.
if ( is.null(opt$mode) ) { opt$mode = 'expression'; }
if ( is.null(opt$clustering) ) { opt$clustering = 'row'; }
if ( is.null(opt$distribution) ) { 
	opt$distribution = 'auto'
} else if ( opt$distribution=='one' ) { 
	opt$distribution = 'one-sided'
} else if ( opt$distribution=='two' ) { 
	opt$distribution = 'two-sided'
}
if ( is.null(opt$norowlabel) ) {
	rowlabel=T
} else {
	rowlabel=F
}								#do not show rowlabels = small pdf
if ( is.null(opt$nocollabel) ) {
	collabel=T
} else {
	collabel=F
}								#do not show col labels = small pdf
if ( is.null(opt$unitlabel) ) { opt$unitlabel='auto'; }
if ( is.null(opt$title) ) { opt$title = 'auto'; }
if ( is.null(opt$subtitle) ) { opt$subtitle = 'auto'; }
if ( is.null(opt$clustdist) ) { opt$clustdist='auto'; }
if ( is.null(opt$clustmethod) ) { opt$clustmethod='auto'; }
if ( is.null(opt$widthfactor) ) { opt$widthfactor = 1; }
if ( is.null(opt$heightfactor) ) { opt$heightfactor = 1; }
if ( is.null(opt$outfile) ) {
	fileprefix = sub("^(.*)[.].*", "\\1", opt$matrix)				#remove suffix from filename (everything after last ".") = output file prefix
	opt$outfile = paste(fileprefix, ".pdf", sep="")
}
if ( is.null(opt$colorpalette) ) { 							#automatic mode
	opt$colorpalette='auto'
} 
if ( is.null(opt$category_colours) ) { opt$category_colours = "preset"; }

##############################################################################################################


library(ComplexHeatmap)
library(gplots)
library(circlize)
library(RColorBrewer)
library(circlize)
library(viridis)	#for color palettes
library(rje)		#for color palettes


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

if ( opt$colorpalette=='auto' ) { 								#automatic mode
	color_vector_onesided=heat
	color_vector_twosided=buwtrd
} else {
	colorpalettemain=sub("^(.*)_r", "\\1", opt$colorpalette)				#remove _r (if it existed)
	palette=get(colorpalettemain)								#get main palette
	if (grepl('_r$', opt$colorpalette)) {							#reverse colors of param ended with _r
		palette=rev(palette)
	} 
	color_vector_onesided=palette
	color_vector_twosided=palette
}




##########################################################################################################################################################
########### read and transform input matrix
##########################################################################################################################################################

setwd(".")
args <- commandArgs(trailingOnly = TRUE)					#parse commandline
data=read.table(file=opt$matrix,sep="\t",header=T,stringsAsFactors=F,row.names=NULL,check.names=FALSE, quote="")	#with header, column 0 = rownames, do not convert strings to character vectors
fileprefix = sub("^(.*)[.].*", "\\1", opt$matrix)				#remove suffix from filename (everything after last ".") = output file prefix
#head(data)

data_temp <- data[!duplicated(data[,1]),]					#remove rows where 1. column is duplicate
data <- data_temp
headline=colnames(data)[-1]
data_temp <- data.frame(data[,-1])								#convert 1. column to rownames
rownames(data_temp) <- data[,1]
colnames(data_temp) <- headline
data <- data_temp

#data = na.omit(data)   							#remove rows with missing values.
data[is.na(data)] <- 0								#replace NA with 0	
#data

if (!is.null(opt$categories)) {
	categories=read.table(file=opt$categories,sep="\t",header=T,stringsAsFactors=T,row.names=1,check.names=FALSE,na.strings=c("","NA"))	#with header, column 0 = rownames, do not convert strings to character vectors
}

##########################################################################################################################################################
########### create heatmap
##########################################################################################################################################################


###############
##  complexheatmap wrapper function
###############

create_complexheatmap=function(m, mode="raw", unitlabel='auto', rowlabel=T, collabel=T, clustering='none', clustdist='auto', clustmethod='auto', distribution='auto', color_vector_onesided=NULL, color_vector_twosided=NULL, optimize=T, title='auto', subtitle='auto')
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
	if ( opt$mode=="raw" ) {						#raw
		if ( unitlabel=='auto' ) {
			unitlabel="Value"
		}
		if ( clustdist=='auto' ) {
			clustdist="euclidean"
		}
		if ( clustmethod=='auto' ) {
			clustmethod="average"
		}
	} else if ( mode=="expression" ) {					#do log2
		m <- m +1							#add pseudocount
		m = log2(m)							#log transform
		if ( unitlabel=='auto' ) {
			unitlabel="Log2 Expression"
		}
		if ( clustdist=='auto' ) {
			clustdist="euclidean"
		}
		if ( clustmethod=='auto' ) {
			clustmethod="average"
		}
	} else if ( mode=="zscore" ) {						#do zscore
		if ( unitlabel=='auto' ) {
			unitlabel="Z-Score"
		}
		m=t(scale(t(m)))										#row zscore
		#m=scale(m)											#column zscore
		if ( clustdist=='auto' ) {
			clustdist="pearson"
		}
		if ( clustmethod=='auto' ) {
			clustmethod="average"
		}
	} else {								#error
		stop("unknown mode:", mode)
	}

	if ( title=="auto" ) { 
		title = sub("^(.*)[.].*", "\\1", opt$matrix)						#remove suffix
		title = basename(title)									#remove path from filename
		opt$title<<-title									#to allow access outside of function (for size computation)
	}
	if ( subtitle=="auto" ) { 
		#subtitle = paste("Mode: ", mode, ", Distance: ", clustdist, ", Clustering: ", clustmethod, sep="")
		subtitle = paste(mode, ", ", clustdist, ", ", clustmethod, sep="")
		opt$subtitle<<-subtitle								#to allow access outside of function (for size computation)
	}

	### mapping of colors to min/max value depending on distribution type (one- or two-sided)
	min=min(m)
	absmin=abs(min(m))
	max=max(m)
	absmax=abs(max(m))
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
		color_vector=color_vector_twosided
	} else {
		color_vector=color_vector_onesided
	}
	message("distribution: ", distribution)
	color_vector_n=length(color_vector)					#number of colors in vector
	maxlimit=max
	minlimit=min
	if (optimize==T && distribution=='two-sided'){				#try to create better color breaks (identical min/max for two-sided, reasonable rounding)
		maxlimit=max(absmax, absmin)
		minlimit=-maxlimit
	}		
	breaks=seq(minlimit,maxlimit,length=color_vector_n)			#one break point for each color; even distribution between min/max
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
		#column_title=NULL,
		column_title=paste(title, subtitle, sep="\n"),
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
		row_dend_width=unit(1,"inches"),
		column_dend_height=unit(1,"inches"),
		row_names_max_width=unit(8,"inches"),
		column_names_max_height=unit(4,"inches"),
		row_names_gp=gpar(fontsize=12),
		column_names_gp=gpar(fontsize=12),
		column_title_gp=gpar(fontsize=10, units='in'),
		#width=unit(3,"inches"),
		heatmap_legend_param=list(
			color_bar="continuous", 			#continuous, discrete
			legend_direction="horizontal"			#horizontal, vertical
		)
	)

	return(ht1)
}


###############
##  create heatmap
###############


ht1=create_complexheatmap(
	m=data,
	mode=opt$mode,								#raw, expression, zscore
	clustering=opt$clustering,						#cluster for 'none', 'row', 'column', 'both'
	clustdist=opt$clustdist,						#auto, euclidean, pearson, spearman, kendall, maximum, manhattan, canberra, binary, minkowski
	clustmethod=opt$clustmethod,						#auto, average, ward.D, ward.D2, single, complete, mcquitty, median, centroid
	unitlabel=opt$unitlabel,						#auto, or string denoting unit label
	rowlabel=rowlabel,							#show row label: T/F
	collabel=collabel,							#show column label: T/F
	distribution=opt$distribution,						#'auto': pick automatically; 'one-sided': eg. expression, 'two-sided': eg. log2fc, zscore
	color_vector_onesided=color_vector_onesided,				#pick color palette for one-sided ditribution (NOT (min<0, max >0))
	color_vector_twosided=color_vector_twosided,				#pick color palette for two-sided ditribution (min<0, max >0)
	optimize=T,								#TRUE: try to create better color breaks/bar (identical min/max for two-sided)
	title=opt$title,
	subtitle=opt$subtitle
)


##########################################################################################################################################################
########### categories = additional colour bars
##########################################################################################################################################################

categories_nr=0
if (!is.null(opt$categories)) {

	#############
	# opt$category_colours=="preset": pre-defined color palettes: one for each category (works for up to 5 different categories, 8-12 elements each)
	#############
	colour_palettes_default=list()							#pre-define palettes used consecutively per category (=column in categories file)
	colour_palettes_default[[1]]=colorRampPalette(brewer.pal(8, "Dark2"))
	colour_palettes_default[[2]]=colorRampPalette(brewer.pal(8, "Accent"))
	colour_palettes_default[[3]]=colorRampPalette(brewer.pal(9, "Set1"))
	colour_palettes_default[[4]]=colorRampPalette(brewer.pal(12, "Set3"))
	colour_palettes_default[[5]]=colorRampPalette(brewer.pal(8, "Set2"))
	
	#colour_palettes_default[[xxx]]=colorRampPalette(c("black", "orange", "yellow"))(3)						#manual colors
	
	
	#head(categories)
	categories_nr=ncol(categories)							#number of categories (=columns)
	categories_names=names(categories)
	cat("\ncategories: ",categories_nr,"(",categories_names,")\n")
	
	cat_list=list()									#list of category column pre-plots
	for(i in 1:ncol(categories)) {							#loop category columns
		cat_name=as.character(names(categories)[i])				#name of category column
		cat_rows=as.character(categories[,i])					#rows of category column
		cat_levels=levels(categories[[cat_name]])
		
		cat_df=data.frame(cat_rows, stringsAsFactors=F)				#create dataframe as input for rowAnnotation (construct needed to permit dynamic legend names)
		colnames(cat_df)=cat_name
	
		#cat("i: ",i,"\n")
		cat("cat_name: ",cat_name,"\n")
		#cat("cat_rows: ",cat_rows,"\n")
		cat("cat_levels: ",cat_levels,"\n")
		cat_levels
		#message(length(cat_levels))
		
		#############
		# fix colours for one category manually here (must be HARDCODED for each plot if desired!)
		#############
		
		if (cat_name=="categoryXYZ") {								#fix color palette for this category
			#num_colors=length(cat_levels)
			#cat("numcolors: ", num_colors, "\n")
			options(warn=-1)
			colour_palette=colorRampPalette(brewer.pal(12, "Dark2"))((length(cat_levels)))	#Dark2, Accent, Set1, ...
			options(warn=0)
			
			colour_vector=c(								#setup colors here
					"immune response"="red",
					"proliferation"="blue",
					#"immune response"=colour_palette[1],
					#"proliferation"=colour_palette[2],
					"angiogenesis"=colour_palette[3],
					"others"=colour_palette[4],
					"blabla"=colour_palette[6]
			)
			colour_map=list()
			colour_map[[cat_name]]=colour_vector
			#print(colour_map)
			
			rowanno = rowAnnotation(df=cat_df, name=cat_name, col=colour_map)
		
		#############
		# random: pick random colors
		#############
	
		} else if (opt$category_colours=="random") {						#random colors
			rowanno = rowAnnotation(df=cat_df, name=cat_name)
	
		#############
		# preset: pick colors from  pre-defined palettes: one for each category (up to 5 different categories, 8-12 elements each)
		#############
	
		} else if (opt$category_colours=="preset") {						#random colors
			colour_palette=colour_palettes_default[[i]](length(cat_levels))					#pick palette for current category
	
			colour_map=list()
			colour_vector=c()
			for(j in 1:length(cat_levels)) {						#loop categories in this column
				cat_level=as.character(cat_levels[j])
				#cat("cat_level: ",cat_level, "\t -> ", colour_palette[j], "\n")
				colour_vector[[cat_level]]=colour_palette[j]
			}
			colour_map=list()
			colour_map[[cat_name]]=colour_vector
			
			rowanno = rowAnnotation(df=cat_df, name=cat_name, col=colour_map)
	
		#############
		# otherwise use selected colorbrewer palette for ALL categories (opt$category_colours)
		#############
	
		} else {										#preset color palette, but random mapping
	
			options(warn=-1)
			colour_palette=colorRampPalette(brewer.pal(12, opt$category_colours))(length(cat_levels))	#Dark2, Accent, Set1, ...
			options(warn=0)
			
			colour_map=list()
			colour_vector=c()
			for(j in 1:length(cat_levels)) {						#loop categories in this column
				cat_level=as.character(cat_levels[j])
				#cat("cat_level: ",cat_level, "\t -> ", colour_palette[j], "\n")
				colour_vector[[cat_level]]=colour_palette[j]
			}
			colour_map=list()
			colour_map[[cat_name]]=colour_vector
			
			rowanno = rowAnnotation(df=cat_df, name=cat_name, col=colour_map)
		}
	
		cat_list[[i]]=rowanno
	}
}


##########################################################################################################################################################
########### combine plots
##########################################################################################################################################################

outlist = ht1

if (!is.null(opt$categories)) {
	for(i in 1:length(cat_list)) {							#loop category pre-plots (one per column)
		rowanno=cat_list[[i]]							#get plot for this category
		outlist = outlist + rowanno						#add category annotation
	}
}

##########################################################################################################################################################
########### calculate necessary pdf size (automatically fits on one page if no rownames are desired) and other params + plot
##########################################################################################################################################################

col_nr=ncol(data)							#nr. of samples
row_nr=nrow(data)							#nr. of genes
col_names_maxlength_label_width=max(sapply(colnames(data),strwidth, units="in", font=12))	#longest column label when plotted in inches	
col_names_maxlength_label_height=max(sapply(colnames(data),strheight, units="in", font=12))	#highest column label when plotted in inches	
row_names_maxlength_label_width=max(sapply(rownames(data),strwidth, units="in", font=12))	#longest row label when plotted in inches	
row_names_maxlength_label_height=max(sapply(rownames(data),strheight, units="in", font=12))	#highest row label when plotted in inches	

#message(paste("col_labelwidth: ", col_names_maxlength_label_width, sep=""))
#message(paste("col_labelheight: ", col_names_maxlength_label_height, sep=""))
#message(paste("row_labelwidth: ", row_names_maxlength_label_width, sep=""))
#message(paste("row_labelheight: ", row_names_maxlength_label_height, sep=""))

#message(opt$title)
#xx=strwidth(opt$title, font=10, units='in')
#message(xx)
#message(opt$subtitle)
#xx=strwidth(opt$subtitle, font=10, units='in')
#message(xx)
title_maxlength_label_width=max(sapply(c(opt$title, opt$subtitle),strwidth, units="in", font=10))
#message(title_maxlength_label_width)

width_margin=0							#buffer for dendrogram + labels
height_margin=0							#buffer for dendrogram + labels + title
if ( !is.null(opt$norowlabel) ) {											#do not show rowlabels = small pdf
	### SMALL
	
	width_margin=width_margin+0.3											#width buffer: labels + small whitespaces
	if (opt$clustering=='both' || opt$clustering=='row') {								#width buffer: dendrogram + small whitespaces between viewports
		width_margin=width_margin +1	
	}
	#message(paste("widthmargin: ",width_margin, sep=""))

	height_margin=height_margin +0.2 +0.5 + (5 * row_names_maxlength_label_height)					#height buffer: labels + small whitespaces + color legend + 2 title rows(+whitespace)
	if (opt$clustering=='both' || opt$clustering=='column') {							#height buffer: dendrogram
		height_margin=height_margin+1
	}
	#message(paste("heightmargin: ",height_margin, sep=""))

	pdf_height=(9+height_margin)											#no row labels = small
	pdf_width=((col_nr*(row_names_maxlength_label_height+0.06))+width_margin+(categories_nr*0.2))			#no row labels = small, categories_nr=number of additional category columns
} else {
	### BIG

	width_margin=width_margin+row_names_maxlength_label_width + 0.3							#width buffer: labels + small whitespaces
	if (opt$clustering=='both' || opt$clustering=='row') {								#width buffer: dendrogram + small whitespaces between viewports
		width_margin=width_margin +1	
	}
	#message(paste("widthmargin: ",width_margin, sep=""))

	height_margin=height_margin +0.2 +0.5 + (5 * row_names_maxlength_label_height)					#height buffer: small whitespaces + color legend + 2 title rows(+whitespace)	
	if (collabel==T) {
		height_margin=height_margin +col_names_maxlength_label_width						#height buffer: labels
	} 
	if (opt$clustering=='both' || opt$clustering=='column') {							#height buffer: dendrogram
		height_margin=height_margin+1
	}
	#message(paste("heightmargin: ",height_margin, sep=""))

	pdf_height=((row_nr*(row_names_maxlength_label_height+0.06))+height_margin)					#readable row labels = big
	pdf_width=((col_nr*(row_names_maxlength_label_height+0.08))+width_margin+(categories_nr*0.2))			#readable row labels = big, categories_nr=number of additional category columns
}

#cat("pdf_width:",pdf_width,"\n")
#cat("title_maxlength_label_width:",title_maxlength_label_width,"\n")

if (pdf_width < title_maxlength_label_width) {										#if the title is longer than the current width make sure to fit the title
	pdf_width = title_maxlength_label_width
}

#cat("pdf_height:",pdf_height,"\n")
#cat("pdf_width:",pdf_width,"\n")

pdf_height=pdf_height * opt$heightfactor
pdf_width=pdf_width * opt$widthfactor

pdf(height=pdf_height, width=pdf_width, opt$outfile)
draw(outlist, heatmap_legend_side="bottom", annotation_legend_side="bottom")
dev.off()

##############################################################################################################
# CLEANUP

if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
