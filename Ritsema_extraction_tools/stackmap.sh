#!/bin/bash
#plot a map of the normed file stack in GMT

#Mapping constraints from Ritsema code

MAP=0/180/3480/6370
MAP1=0/180/0/6370
PROJ=P7i
bord=-Bf30/f500
icont=1
#The data has been normalized to 1 and added - usually we don't see anomalies larger than 1
mxc=1

LCL=$PWD
BIN=$LCL/bin

GMT=gmt4

if [ -z "$1" ]; then
	echo 'Usage: ./stackmap.sh dataset [grd file]'
	exit 1
fi 

dataset=$1

$GMT gmtset LABEL_FONT_SIZE 24
$GMT gmtset ANNOT_FONT_SIZE_PRIMARY 8
$GMT gmtset PAPER_MEDIA letter

#Create the color scheme
echo $LCL/colourscales/gu12.chj      >  in_mkrb
echo $icont                            >> in_mkrb
echo $mxc                              >> in_mkrb
$BIN/mkrb_jr < in_mkrb > col.cpt

$GMT grdimage $dataset -Ccol.cpt -X2 -Y1.8 $bord -J$PROJ -R$MAP -P -K > slice.ps

#-- Horizons
$GMT psxy $LCL/horizons/core.xy  -R$MAP1 -J$PROJ -W2ta            -A -O -K -P >> slice.ps
$GMT psxy $LCL/horizons/670km.xy -R$MAP  -J$PROJ -W2t12_12_12:12  -A -O -K -P >> slice.ps
$GMT psxy $LCL/horizons/1000km.xy -R$MAP  -J$PROJ -W2t12_12_12:12  -A -O -K -P >> slice.ps

#-- colourscale
$GMT psscale -D3.0/0.0/4.0/0.8h -Ccol.cpt -Y12.0 -X6 -O -P -Ba0.5f0.5g0.5/:"Stacked velocity perturbation": >> slice.ps

$GMT ps2raster -Tf slice.ps
open slice.pdf

