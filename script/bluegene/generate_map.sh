#!/bin/bash -e 
readonly MAPFILE="map.txt"
temp=`mktemp`
for x in {0..7}
do
    for y in {0..7} 
    do
        for z in {0..7}
        do
            t=$(($RANDOM%4))
                echo $x" "$y" "$z" "$t >> $temp
            done
        done
    done

cat $temp | sort -R > map.txt
rm $temp                
