#!/bin/bash
cd $1
eval `scramv1 runtime -sh`
root -b LHCO$2$3.C
rm LHCO$2$3.C
exit 0

