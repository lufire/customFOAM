#!/bin/bash

#- define variables
Folder=$1
echo 'Utilized folder : ' $Folder
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
#cd $casePath/Animation
rm -r temp 
mkdir temp
cp $Folder/*.jpg temp/. 
#cp Animation/*.jpg Animation/temp/.
# mogrify -resize 800x600 temp/*.jpg
 mogrify -crop 320x300+180+30 temp/*.jpg
# convert $Folder/*.jpg -crop 320x300+180+30 temp/temp.jpg
ffmpeg -r 25 -pattern_type glob -i "temp/*.jpg" -b:v 2M -vcodec msmpeg4 $Folder.wmv
rm -r temp 

