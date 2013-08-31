#! /bin/bash

while read dataset output type
do
  sed -e "s|-dataset-|$dataset|g" -e "s|-output-|$output|g" pattuple_${type}WithoutTriggerInfo_Template.cfg > pattuple.cfg
  crab -create -cfg pattuple.cfg
  ./submittingJobsCrab.sh ../resultsCrab/pattuple_$output
  rm pattuple.cfg 
done < $1

exit 0

