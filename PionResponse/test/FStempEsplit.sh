#!/bin/bash

DIR=`echo $1`
ENERGY=`echo $2`
SPLIT=`echo $3`
PART=`echo $4`

mkdir ${DIR} -p
cat ./SinglePi_SIM_split_temp.py \
| sed -e s/ENERGYIN/${ENERGY}/ \
| sed -e s/SPLITNUMBER/${SPLIT}/ \
| sed -e s/PNUMBER/${PART}/ \
> ./${DIR}/SinglePigun_FULLSIM_${ENERGY}_split${SPLIT}_part${PART}_cfg.py
