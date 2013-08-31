
#! /bin/bash                                                                                                                                                                                   
BASEDIR=/afs/cern.ch/work/j/jgomezca/public/ForUniandes/CMSSW_5_3_9/src
while read dataset output type
  do
  echo $output
  LIST=$(ls ../resultsLocal/LHCO_${output}|grep LHCO_${type}) 
  mkdir ../resultsLocal/LHCO_${output}/Splitted
  COUNTER=0
  for FILE in $LIST; do
      echo $FILE   
      sed -e "s|-input-|../resultsLocal/LHCO_${output}/$FILE|g" -e "s|-output-|../resultsLocal/LHCO_${output}/Splitted/|g" -e "s|-name-|${type}${COUNTER}|g" LHCO_Template.C > LHCO${type}${COUNTER}.C
      bsub -q 8nh runningLHCO_Template.sh $BASEDIR $type $COUNTER
      
      
      let COUNTER=COUNTER+1
# if [  $COUNTER -gt 11 ] 
#  then
#    echo "Well done"
#    break
# fi
  done
done < $1
exit 0

