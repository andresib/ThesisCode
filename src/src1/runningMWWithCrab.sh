#! /bin/bash

while read dataset output type
do
  mv ../resultsLocal/LHCO_${output}/Splitted data/LHCO_${output}
  sed -e "s|-input-|$output|g" mw_with_grid_template.py > mw_with_grid.py
  sed -e "s|-dataset-|$dataset|g" -e "s|-output-|$output|g" MadWeight_${type}_Template.cfg > madweight.cfg
  crab -create -cfg madweight.cfg
  ./submittingJobsCrab.sh ../resultsCrab/mw_$output 
  mv data/LHCO_${output} ../resultsLocal/LHCO_${output}/Splitted
done < $1

exit 0

