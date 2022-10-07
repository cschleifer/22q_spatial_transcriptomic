#!/bin/bash
#DIRECTORY=$(cd `dirname $0` && pwd)
#echo $DIRECTORY

# set path to workbench command binary
wbcommand "/Applications/workbench/bin_macosx64/wb_command"
checkwb=$(which wb_commandd)
if [ -z $checkwb ]; then 
    echo "wb_command path not set, adding to PATH:"
    export PATH=$PATH:${wbcommand}
    echo $PATH
else 
    echo "wb_command path already set: ${check_wb}" 
fi

wb_command

project=$(git rev-parse --show-toplevel)
echo $project

