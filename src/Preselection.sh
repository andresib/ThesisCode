#! /bin/bash                                                                                                                                                                                   
pwd
cd /afs/cern.ch/work/j/jgomezca/public/ForUniandes/CMSSW_5_3_9/src
eval `scramv1 runtime -sh`
echo $2
if [ "$3" == 'MC' ];
    then
    label="preselectionStudy"
fi
if [ "$3" == 'data' ];
    then
    label="preselection"
fi
echo ${label}_$2
pwd
mkdir ../resultsLocal/${label}_$2/
mkdir ../resultsLocal/LHCO_$2/
sed -e "s|-input-|../resultsCrab/pattuple_$2/res/$4|g" -e "s|-outputLHCO-|../resultsLocal/LHCO_$2|g"  -e "s|-output-|../resultsLocal/${label}_$2|g" -e "s|-postfix-|$4|g"  ${label}_template.py > ${label}_$2_$4.py
pwd 
ls -l ${label}_$2_$4.py
#cd /afs/cern.ch/work/j/jgomezca/public/ForUniandes/CMSSW_5_3_9/src
#eval `scramv1 runtime -sh`
pwd
ls -l ${label}_$2_$4.py
cmsRun ${label}_$2_$4.py
pwd
ls -l ${label}_$2_$4.py
rm ${label}_$2_$4.py
exit 0

