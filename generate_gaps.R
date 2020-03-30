#############################################################################################
#### Script for delineating gaps from treatment simulations for RESTSIM
#### Author: Jeffery Cannon (Jeffery.Cannon@jonesctr.org)
#### Institution: The Jones Center at Ichauway, Newton, GA
#### Date Created: 06/26/2018
#### Last Modified: 03/30/2020
#############################################################################################
# The goal of this script is to take as an input a modified LANDFIRE canopy cover layer and
# a threshold canopy cover (e..g, 10%) and output a raster of all areas with canopy cover
# lower than the threshold for further analysis of spatial heterogeneity
# Estimated analysis time is approximately 4 minutes per raster for HUC 8 watershed
#############################################################################################

#######################################START SET UP##########################################
packages <- c('raster', 'sf')
for (package in packages) {
  if (suppressMessages(!require(package, character.only = T))) {
    install.packages(package, repos = 'https://cran.mtu.edu/')
    suppressMessages(library(package, character.only = T))
  }
}

#---> Load data (cite sources)
# Load treatment canopy cover raster output from simulation (from simulations)
CC <- raster('data/canopy_cover_raster.tif')

# Load LANDFIRE Existing Vegetation Type and attribute table
EVT <- raster('data/Landfire_EVT.tif')
EVT_attrib <- read.csv('data/EVT_09152016_attribute.csv')

# study boundaries to shapefile
bound <- st_read('data/boundary.shp')

# Layer to use as default projection
default_layer = CC
default_proj = proj4string(default_layer)

#---> Analysis parameters
gap_thresh = 10 # Mark as "low" canopy cover all areas <= 10% CC
rasterOptions(maxmemory = 1 * 10 ^ 9) # for faster processing

#########################
# FUNCTION: rasterToGapPolygon()
#
# This function takes a binary raster (higher than treshold = 0; below threshold = 1) and delineates
# continuous "meadows". The script requires installation of ArcGIS with spatial analyst tool,
# Python script that is called was written by Gannon and Cannon
#########################
rasterToGapPolygon <- function(in_raster)
{
  dir.create('data/scratch')
  writeRaster(in_raster,
              'data/scratch/tmp_raster',
              format = 'ascii',
              overwrite = TRUE) #write raster
  cwd <- getwd()
  cwd <- gsub('/', '\\', cwd, fixed = TRUE)
  py_exe_path <- 'C:/Python27/ArcGIS10.3/python.exe'
  py_script_path <- 'raster2polygon.py'
  sys_call <- paste(py_exe_path, py_script_path, cwd, sep = ' ')
  system(sys_call) #call python script with argument of pathname
  cat('Loading file into R memory...\n')
  gap_poly <- st_read('data/scratch/tmp_poly.shp') #load polygon
  gap_poly <- subset(gap_poly, GRIDCODE == 1) # include only gap areas
  cat('Deleting temporary shapefile\n')
  unlink('data/scratch', recursive = TRUE) #cleanup temp files
  return(gap_poly)
}

#########################
# FUNCTION: CreateGapPolys()
#
# This function takes a continous raster of canopy cover, a canopy cover threshold, study boundary,
# and mask and delineates gaps using the threshold. This is a wrapper function provides necessary
# pre- and post-processing for the rasterToGapPolygons() function which primarily preps input for
# Pyton. The pre- and post-processing includes (1) thresholding continous canopy cover data, (2)
# masking (e.g., only including forested areas), (3) disaggreating resulting polygons, and (4)
# defining projection of  the resulting polygons. See rasterToGapPolygons() for software and file
# structure requirements.
#########################

CreateGapPolys <-
  function(CC,
           gap_thresh,
           mask,
           default_proj,
           boundary) {
    cat('Masking binary gaps\n')
    CC_gap = (CC <= gap_thresh) * crop(mask, CC)
    cat('Finalizing gap raster\n')
    CC_gap[CC_gap == 0] <- NA
    cat('Running ArcGIS/Python raster to polygon tool...\n')
    CC_gap_poly = rasterToGapPolygon(CC_gap)
    cat(paste(
      'Feature currently contains',
      dim(CC_gap_poly)[1],
      'polygons\n',
      sep = ' '
    ))
    cat(paste('Defining polygon CRS:\n', default_proj, '\n'))
    st_crs(CC_gap_poly) <- CRS(default_proj)
    cat(paste('Disaggregating polygons...\n'))
    CC_gap_poly= st_cast(CC_gap_poly, 'POLYGON')
    cat(paste('Clipping away features not touching study boundary\n'))
    
    boundary = st_transform(boundary, CRS(default_proj))
    CC_gap_poly = CC_gap_poly[boundary,]
    cat(paste('Output feature contains', dim(CC_gap_poly[1]), 'polygons\n', sep = ' '))
    return(CC_gap_poly)
  }
########################################END SET UP###########################################

######################################START ANALYSIS#########################################
start_time = Sys.time()
#---> Create raster of "vegetated" lands, ommitting, "Unvegetated", and "No dominated lifeform" (i.e., developed)
EVT_rcl    <- data.frame(EVT_ORDER = levels(EVT_attrib$EVT_ORDER), EVT_VEGETATED = c(1, 1, 0, 0, 0, 0, 1, 1)
)
EVT_attrib <- merge(EVT_attrib, EVT_rcl)
EVT_attrib <- subset(EVT_attrib, EVT_VEGETATED == 1)
VEG <- EVT %in% EVT_attrib$VALUE

gap_poly = CreateGapPolys(
  CC = CC,
  gap_thresh = gap_thresh,
  mask = VEG,
  default_proj = default_proj,
  boundary = bound
)

#######################################END ANALYSIS##########################################

######################################START OUTPUTS##########################################
#---> Write all outputs to disk
st_write(gap_poly, dsn = 'data/gaps_initial.shp', driver = 'ESRI Shapefile', delete_layer =  TRUE)
end_time = Sys.time()
print(end_time - start_time)
#######################################END OUTPUTS###########################################

#####################################START GRAPHICS##########################################

#######################################END GRAPHICS##########################################
