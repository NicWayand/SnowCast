#!/bin/bash

# Input
indir=$1
var=$2

cd $indir
shopt -s nullglob
for f in *${var}*.tif
do
    echo $f
    # Create KML
    #echo /create_KML.sh $f col_snowdepthavg.txt
    # Zip to KMZ
    #f_kml="${f%.*}"
    #f_o=$f_kml'.kmz'
    #echo zip -r $f_o $f_kml
    # Upload
    #echo ../ftp_upload.sh $f_o KML 
done

