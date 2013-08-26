#! /bin/bash                                                                                                                                                                                   
BASEDIR=/afs/cern.ch/work/j/jgomezca/public/ForUniandes/CMSSW_5_3_9
while read dataset output type
  do
  mkdir ../resultsLocal/mw_$output
  echo $output
  LIST=$(ls ../resultsLocal/LHCO_${output}/Splitted)
  COUNTER=0
  for FILE in $LIST; do
          if [  $COUNTER -gt 20 ]  
	      then
	      echo $FILE
	      cp -r Analysis Analysis$FILE
	      echo "MW copied"
	      cp ../resultsLocal/LHCO_${output}/Splitted/$FILE Analysis$FILE/Events/input.lhco
	      echo "LHCO copied"
	      bsub -q 8nm runningMW.sh $BASEDIR $FILE $output 
	      echo "done"
	      if [  $COUNTER -gt 1000 ]                                                                                                                         
		  then                          
		  echo "Well done"    
		  break    
	      fi  
	  fi
	  let COUNTER=COUNTER+1
  done
done < $1
exit 0

