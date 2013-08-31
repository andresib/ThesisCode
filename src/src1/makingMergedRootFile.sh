#! /bin/bash                                                                                                                                                                                   
while read dataset input type
do
 bsub -q 8nh MergingRootFiles.sh $input
 echo $input
done < $1
exit 0

