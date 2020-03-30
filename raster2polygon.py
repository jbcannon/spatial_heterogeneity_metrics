# Created by: Jeff Cannon/Ben Gannon
# Created on: 07/26/2017
# Last Updated: 07/26/2017

'''
The purpose of this script is to convert a raster input file
into a polygon using ArcGIS
'''

#-> Import modules
import arcpy, os, sys, traceback, shutil
from arcpy.sa import *

# Check out spatial analyst extension
arcpy.CheckOutExtension('spatial')

path=sys.argv[1]
print(path)
arcpy.env.workspace = path + r'\data\scratch'
arcpy.env.scratchWorkspace = path + r'\data\scratch'
print(arcpy.env.workspace)

#-> Create variables

# Input
tmp_raster_asc = arcpy.env.workspace + r'\tmp_raster.asc'

# Output
tmp_poly_fc = r'\tmp_poly.shp'

###---> Main body of analysis

#-> Convert float input into integer raster
try:
    raster_int = Int(tmp_raster_asc)
    print('Converted float input into integer ')
except Exception as err:
    print(err.args[0])

# Process: Raster to Polygon
try:
    if arcpy.Exists(tmp_poly_fc):
        arcpy.Delete_management(tmp_poly_fc)
    arcpy.RasterToPolygon_conversion(raster_int, tmp_poly_fc, "SIMPLIFY", "VALUE")
    print('Converted raster to polygon')
except Exception as err:
    print(err.args[0])
print('Polygon shapefile created')
'''
'''
