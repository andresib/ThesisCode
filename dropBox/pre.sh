DOMAINNAME=`dnsdomainname`
if [[ ${DOMAINNAME} == "cern.ch" ]] 
then
    export SCRAM_ARCH=slc5_amd64_gcc472
    CMSSW_RELEASE=CMSSW_6_2_0_pre4
    source /afs/cern.ch/cms/cmsset_default.sh
    CMSSW_REL=/afs/cern.ch/cms/$SCRAM_ARCH/cms/cmssw/$CMSSW_RELEASE/
elif [[ ${DOMAINNAME} == "cms" ]]
then
    export SCRAM_ARCH=slc5_amd64_gcc462
    CMSSW_RELEASE=CMSSW_5_3_6
    source /data/cmssw/cmsset_default.sh
    CMSSW_REL=/data/cmssw/$SCRAM_ARCH/cms/cmssw/$CMSSW_RELEASE/
fi
pushd $CMSSW_REL >/dev/null
eval `scramv1 runtime -sh`
popd >/dev/null

# Enable the new key-based authentication
export COND_AUTH_SYS=1

