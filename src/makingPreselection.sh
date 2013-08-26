#! /bin/bash                                                                                                                                                                                   
while read dataset input type
do
  cd ../resultsCrab/pattuple_$input/res/
  LIST=$(ls *.root)
  cd ../../../src
  for FILE in $LIST; do
      bsub -q 8nh Preselection.sh $dataset $input $type $FILE
  done
done < $1
exit 0

