#!/bin/bash

# Move test data to the correct location

helpFunction()
{
   echo ""
   echo "Usage: $0 -a parameterA"
   echo -e "\t-a choose to setup to run either a single sample analysis (-a single) or a cohort analysis (-a cohort)"
   exit 1 # Exit script after printing help
}

while getopts "a:" opt
do
   case "$opt" in
      a ) parameterA="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterA" ]
then
   echo "Value for parameter not provided";
   helpFunction
fi

# Begin script in case all parameters are correct
cp -r ./test/"$parameterA"/* ../