#!/bin/bash
cd $1/src/Analysis$2/
eval `scramv1 runtime -sh`
./bin/madweight.py 
sed -e "s|-input-|$2|g" -e "s|-output-|$3|g" ../storingBestPerm_Template.C > storingBestPerm.C
root -b storingBestPerm.C
#rm -r storingBestPerm.C
cd ..
#rm -r Analysis$2
exit 0

