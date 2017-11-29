#!/bin/bash
set -e

# Converts grib2 files to netcdf using a subset of variables
# Note: source activate ncl
# Note: Hardcoded extent to subset to below

source activate ncl

base_dir=$1 # This should be the dir that contains grib2 and netcdf
folderin=$base_dir'/grib2' #folder where the grib are
folderout=$base_dir'/netcdf' #folder to store netcdf

cd $folderin # cd to that folder

#now loop through all the grib files in there
shopt -s extglob
for f in HRRR*.grib2
do
    # Subset to region of interest
    filename=$(basename "$f")
    small_f="subROI_"$filename
    #echo $small_f
    wgrib2  $f -small_grib -121.48187300000001:-120.336557 46.44775:47.01475 $small_f
    # Clean up
    #rm -f $small_f
done

# Concatenate
echo "Starting concat"
cat subROI_* > ../temp.grib2
echo "Finished concat"

# to netcdf
ncl_convert2nc ../temp.grib2 -o "../"$folderout -itime
# Clean up
rm -f ../temp.grib2
echo "Finished converting to netcdf"



