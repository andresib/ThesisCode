#! /bin/bash

ls ../resultsCrab/pattuple_$1/res/*.root                                                                                                                             
LIST=$(ls ../resultsCrab/pattuple_$1/res/*.root|tee ../resultsCrab/pattuple_$1/res/list.txt)
sed -e "s|-input-|../resultsCrab/pattuple_$1/res|g" -e "s|-name-|$1|g" AddingFilePrefix_Template.C > AddingFilePrefix_$1.C
root -b AddingFilePrefix_$1.C
rm ../resultsCrab/pattuple_$1/res/list.txt
rm AddingFilePrefix_$1.C
cd /afs/cern.ch/work/j/jgomezca/public/ForUniandes/CMSSW_5_3_9/src
eval `scramv1 runtime -sh`
edmCopyPickMerge inputFiles_load=../resultsCrab/pattuple_$1/res/listWithPrefix.txt outputFile=../resultsCrab/pattuple_$1/output.root
exit 0

