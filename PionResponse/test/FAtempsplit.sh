#!/bin/bash

DIR=`echo $1`
ENERGY=`echo $2`
SPLIT=`echo $3`

mkdir ${DIR} -p
cat ./fullsimpionanalyzer_split${SPLIT}_temp_cfg.py \
| sed -e s/ENERGYIN/${ENERGY}/ \
> ./${DIR}/fullsimanalyzer_${ENERGY}_split${SPLIT}_cfg.py

#cmsRun ${DIR}/fullsimanalyzer_${ENERGY}_split${SPLIT}_cfg.py
