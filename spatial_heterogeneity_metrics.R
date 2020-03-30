#############################################################################################
#### Calcuate spatial heterogeneity metrics for CEAP treatment simulations
#### Author: Jeffery Cannon (Jeffery.Cannon@jonesctr.org) 
#### Institution: The Jones Center at Ichauway, Newton, GA
#### Date Created: 03/30/2020
#### Last Modified: 08/13/2018
#############################################################################################
# The goal of this script is to input canopy cover data from measured or simulated
# treatment outputs and calculate spatial heterogeneity metrics. Data inputs canoy cover raster,
# catchment boundaries, thresholds. Outputs include 

#[].. finish description
# Include a CC input, catchment input

#...contagion index, cover diversity (H),  coverage of canopy cover classes (p); and semi-
# variogram information for all three treatments. Outputs include statistical summaries,
# analysis, and figures for spatial heterogeneity portion of study
#############################################################################################

#######################################START SET UP##########################################
#---> Load libraries
packages <-
  c('sf', 'raster', 'units', 'plyr', 'doParallel', 'foreach', 'rgdal','rgeos')
for (package in packages) {
  if (suppressMessages(!require(package, character.only = T))) {
    install.packages(package,repos='https://cran.mtu.edu/')
    suppressMessages(library(package, character.only = T))
  }
}

#---> Load data (cite sources)
# Raster data of canopy cover
cc = raster('data/canopy_cover_raster.tif')

# Shapefile data containing all boundaries within which to calculate metrics
bound = st_read('data/boundary.shp')

# Shapefile data containing all potential gaps created using "generate_gaps.R"
in_gaps = st_read('data/gaps_initial.shp')

# Location for output data
output_location ='data/results.csv'

#---> Load custom functions with explanations

#########################
# FUNCTION: lgGapDelin()
#
# This function takes shapefile of potential gaps, and extracts those that meet a size requirement
# of greater than the specified dimension (60 m diameter by default)
#########################
lgGapDelin <- function(gaps, catches, minDimension_m = 60) {
  # For efficiency, delete all gaps that cannot meet minimum dimension based on size, if given
  cat('.Calculating area for all input gaps\n')
  gaps$total_gap_area_ha = st_area(gaps, by_element = TRUE)
  units(gaps$total_gap_area_ha) = 'ha'
  min_radius_m = minDimension_m / 2
  min_area = set_units(min_radius_m ^ 2 * pi,'m^2')
  
  cat('.Removing gaps smaller than minimum dimension\n')
  gaps = subset(gaps, gaps$total_gap_area_ha > min_area)
  
  cat('.Repairing any geometry errors\n')
  gaps = st_buffer(gaps, dist = 0)
  
  #For remaining gaps, complete negative then positive buffer for minimum dimension
  cat('.Executing negative buffer\n')
  gaps = st_buffer(gaps, dist = -min_radius_m)
  cat('.Executing positive buffer\n')
  gaps = st_buffer(gaps, dist = min_radius_m)
  
  cat('.Repairing any geometry errors\n')
  gaps = st_buffer(gaps, dist = 0)
  cat('.Removing gaps outside of study area\n')
  catches = st_union(catches)
  catches = st_buffer(catches, dist = 150)
  catches = st_buffer(catches, dist = -150)
  gaps = st_intersection(gaps, catches)
  
  #recalculate gap area after disaggregating
  cat('.Disaggregating gaps\n')
  gaps = st_cast(gaps, 'POLYGON')
  cat('.Recalculating gap area\n')
  gaps$total_gap_area = st_area(gaps, by_element = TRUE)
  units(gaps$total_gap_area) = 'ha'
  cat('.Removing small gaps\n')
  gaps = subset(gaps, gaps$total_gap_area > min_area)
  cat('.Delineated', nrow(gaps), 'gaps\n')
  return(gaps)
}

#########################
# FUNCTION: getHeterogeneityStats()
#
# This function takes a canopy cover raster, gap shapefile, and number of canopy cover category breaks (default = 5)
# and outputs spatial heterogeneity metrics for the area. Statistics calculated include canopy cover distribution,
# gap statistics, landscape diversity, and relative contagion index. See support functions below for detail.
#########################
getHeterogeneityStats <-
  function(catch,
           CC,
           gaps,
           breaks = 5,
           max = 100,
           min = 0) {
    a = getCC_stats(CC)
    b = getGapCover(catch, gaps)
    CC_cls = classifyCC(CC, breaks = breaks)
    c = getP(CC_cls,
             max = max,
             min = min,
             breaks = breaks)
    d = getDiversity(CC_cls, breaks = breaks)
    e = ContagionCalc(CC_cls)
    x = list(a, b, c, d, e)
    x = do.call(cbind, x)
    return(x)
  }

#########################
# FUNCTION: getCC_stats()
#
# This function takes a canopy cover raster, and calculcates mean, standard devation, and coefficient of variation
# of canopy cover within the area included.
#########################
getCC_stats = function(CC) {
  mean = cellStats(CC, stat = 'mean')
  sd = cellStats(CC, stat = 'sd')
  cv = sd / mean
  return(data.frame(
    cc_mean = mean,
    cc_sd = sd,
    cc_cv = cv
  ))
}

#########################
# FUNCTION: getGapCover()
#
# This function takes a boundary shapefile (catchment), and a gap feature and calculates the total percent cover
# for gaps within the boundary/catchment
#########################
getGapCover = function(catch, gaps) {
  g = crop(gaps, catch)
  if (length(g) == 0) {
    return(data.frame(gap_cover = 0))
  } else {
    g_area = gArea(g)
    c_area = gArea(catch)
    g_cover = g_area / c_area * 100
    return(data.frame(gap_cover = g_cover))
  }
}

#########################
# FUNCTION: classifyCC()
#
# This function takes a canopy cover raster, and number of breaks and reclassifies the raster into equal bins
# assuming canopy cover ranges from 0 to 100%.
#########################
classifyCC = function(CC, breaks) {
  breakpts = seq(0, 100, by = 100 / breaks)
  rcl = matrix(c(breakpts[1:breaks],
                 breakpts[2:(breaks + 1)],
                 0:(breaks - 1)),
               ncol = 3)
  ras_cls = reclassify(CC, rcl)
  return(ras_cls)
}

#########################
# FUNCTION: getP()
#
# This function takes a classified canopy cover raster (see classifyCC(), and its associated max, min, and breaks
# and calculautes the percent cover for each category, including 0 for any missing categories
#########################
getP = function(CC_cls, min, max, breaks) {
  breakpts = seq(min, max, by = max / breaks)
  breakpts = data.frame(start = breakpts[1:breaks], end = breakpts[2:(breaks +
                                                                        1)])
  colNames = paste('p', breakpts$start, '_', breakpts$end, sep = '')
  blank_freq = data.frame(value = c(0:(breaks - 1), NA))
  freq = as.data.frame(freq(CC_cls))
  freq = merge(blank_freq, freq, all.x = TRUE)
  freq$count[is.na(freq$count)] <- 0
  freq = subset(freq, !is.na(value))
  freq$p = freq$count / sum(freq$count) * 100
  rownames(freq) <- colNames
  freq = as.data.frame(t(freq))
  freq = freq['p', ]
  rownames(freq) = 1:dim(freq)[1]
  return(freq)
}

#########################
# FUNCTION: getDiversity()
#
# This function takes a classified canopy cover raster (see classifyCC()), and its associated breaks and calculates
# landscape/diversity using Shannon's diversity index H. See Cannon et al. 2020 for details.
#########################
getDiversity = function(CC_cls, breaks) {
  p = getP(CC_cls, max = 100, min = 0, breaks)
  p = as.vector(t(p)) # convert from df to vector
  p = p[p != 0]
  classes = length(p)
  p = p / 100 #conver to proportion
  num  = sum(p * log(p))
  Hmax = log(classes)
  H = -num / Hmax
  return(data.frame(H = H))
}

#########################
# FUNCTION: getAdjMatrix()
#
# This function takes a classified canopy cover raster (see classifyCC()), and calculates an adjacency matrix used
# as an input to the relative contagion calculation (see ContagionCalc()). The result is a table showing the
# proportion of pixels of class i that are adjacent to class j for all combinations of classes.
#########################
getAdjMatrix = function(CC_cls) {
  # Create class names and empty transition table
  classes = unique(CC_cls)
  n_classes = length(classes)
  full_table = expand.grid(as.factor(classes), as.factor(classes))
  colnames(full_table) <- c('from_cls', 'to_cls')
  # Get list of cells that != NA
  all_cells = Which(!is.na(CC_cls), cells = TRUE)
  # Create cell adjacency table
  adj_matrix = adjacent(CC_cls, all_cells, pairs = TRUE, target = all_cells)
  adj_matrix = as.data.frame(adj_matrix)
  colnames(adj_matrix) = c('from_cell', 'to_cell')
  # Add values from raster to adjacency table
  adj_matrix$from_cls = CC_cls[adj_matrix$from_cell]
  adj_matrix$to_cls = CC_cls[adj_matrix$to_cell]
  # Summarize transition probabilities
  adj_matrix = ddply(adj_matrix, .(from_cls, to_cls), summarize, n = length(to_cls))
  adj_matrix = merge(full_table, adj_matrix, all.x = TRUE)
  adj_matrix$n[is.na(adj_matrix$n)] = 0
  adj_totals = ddply(adj_matrix, .(from_cls), summarize, t = sum(n))
  adj_totals$Pi = with(adj_totals, t / sum(t))
  adj_matrix = merge(adj_matrix, adj_totals)
  adj_matrix$Pji = with(adj_matrix, n / t)
  adj_matrix$t <- NULL
  adj_matrix$n <- NULL
  return(adj_matrix)
}


#########################
# FUNCTION: ContagionCalc()
#
# This function takes a classified canopy cover raster (see classifyCC()), and calculates relative contagion
# Using the method described in Li and Reynolds 1993, Landscape Ecology (RC2). 
#########################
ContagionCalc = function(CC_cls) {
  #Following Li and Reynolds 1993, Landscape Ecology
  adjMatrix = getAdjMatrix(CC_cls)
  n = length(unique(adjMatrix$from_cls))
  Pij = with(adjMatrix, Pi * Pji)
  RC2 = 1 + sum(Pij * log(Pij), na.rm = TRUE) / (2 * log(n))
  return(data.frame(RC2 = RC2))
}


#########################
# FUNCTION: loopHeterogeneityStats()
#
# This function takes takes a list of feature boundaries (catches), with the id column specified, as well as a 
# canopy cover raster and gap shapefile, and breaks. The function loops through each catchment and, using 
# parallel processing, generates output on canopy cover statistics, gap statistics, landscape diversity and contagion
#########################
loopHeterogeneityStats <-
  function(catches,
           id = 'FEATURE',
           CC,
           gaps,
           breaks = 10,
           max = 100,
           min = 0) {
    #Set up blank data.frame for parallel processing
    out_df = data.frame(
      id = NA,
      cc_mean = NA,
      cc_sd = NA,
      gap_cover = NA
    )
    colnames(out_df)[1] = id
    break_names = seq(min, max, by = max / breaks)
    break_names = paste('p', break_names[1:breaks], '-', break_names[2:length(break_names)], sep = '')
    out_df[, break_names] <- NA
    out_df[, c('H', 'RC2')] <- NA
    out_df = subset(out_df, FALSE)
    for (i in catches@data[, id]) {
      tmp_catch = subset(catches, catches@data[, id] == i)
      tmp_CC = mask(crop(CC, tmp_catch), tmp_catch)
      tmp_df = getHeterogeneityStats(
        tmp_catch,
        tmp_CC,
        gaps,
        breaks = breaks,
        max = 100,
        min = 0
      )
      name_df = data.frame(id = tmp_catch@data[, id])
      colnames(name_df) = id
      tmp_df = cbind(name_df, tmp_df)
      out_df = rbind(out_df, tmp_df)
    }
    return(out_df)
  }


########################################END SET UP###########################################

######################################START ANALYSIS#########################################

#---> Delineate large gaps for each treatment
lg_gaps = lgGapDelin(in_gaps, bound, minDimension_m = 60)
lg_gaps = as(lg_gaps,'Spatial')
bound = as(bound, 'Spatial')
lg_gaps$total_gap_area = as.vector(lg_gaps$total_gap_area)
lg_gaps$total_gap_area_ha = as.vector(lg_gaps$total_gap_area_ha)

# Calculate gap cover, cc metrics, and heterogeneity metrics
stats = loopHeterogeneityStats(bound, id = 'FEATURE', cc, lg_gaps)

#######################################END ANALYSIS##########################################

######################################START OUTPUTS##########################################
### Write metrics to disk
write.csv(stats, output_location, row.names = FALSE)
print(stats)
cat('Output written to ', output_location, '\n')
#######################################END OUTPUTS###########################################

