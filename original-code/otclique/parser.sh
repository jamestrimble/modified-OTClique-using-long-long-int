#!/bin/bash
if [ $# -ne 1 ]
then
	printf "graph,limit,Pre[sec],BB[sec],BB(iterations),Total[sec],weight\n"
	exit 0
fi

graph_file=$1

n=`grep 'p edge' $graph_file | sed 's/  */ /g' | cut -d ' ' -f3`
if [ $n -le 1500 ]; then
    block_size_limit=25
else
    block_size_limit=20
fi

ifs_origin=$IFS 
IFS=$'\n'
result=(`./otclique $graph_file $block_size_limit`)
IFS=$ifs_origin

# Graph file name
printf "%s," "${graph_file##*/}"
# Subset size limit
set -- ${result[0]}; printf "%s," $5
# Precomputation phase time
set -- ${result[2]}; printf "%s," $4
# Branch-and-bound time
set -- ${result[3]}; printf "%s," $4
# Branch-and-bound iterations
set -- ${result[4]}; printf "%s," $4
# Total time
set -- ${result[5]}; printf "%s," $4
# Maximum weight
set -- ${result[6]}; printf "%s\n" $4
