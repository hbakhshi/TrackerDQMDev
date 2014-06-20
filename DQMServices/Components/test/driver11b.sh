#!/bin/bash

eval `scramv1 r -sh`

#dqmHarvesting

tnum=11
DQMSEQUENCE=HARVESTING:DQMHarvestCommonSiStripZeroBias
step=b

cmsDriver.py test_${tnum}_${step}_1 -s ${DQMSEQUENCE} --conditions auto:com10 --data --filein file:test_${tnum}_a_1_RAW2DIGI_RECO_DQM.root --scenario pp --customise DQMServices/Components/test/customHarvesting.py --no_exec --python_filename=test_${tnum}_${step}_1.py
cmsRun -e test_${tnum}_${step}_1.py &> q${tnum}.1.log 

if [ $? -ne 0 ]; then
  return 1
fi

mv FrameworkJobReport{,_${tnum}_${step}_1}.xml

