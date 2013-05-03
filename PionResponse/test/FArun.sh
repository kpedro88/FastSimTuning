#!/bin/bash

DIR=pion

for ENERGY in 1 2 3 5 9 10 11 15 20 30
  do
    ./FAtempsplit.sh ${DIR} ${ENERGY} 10000
  done

for ENERGY in 50 100 150 225 300
  do
    ./FAtempsplit.sh ${DIR} ${ENERGY} 5000
  done

for ENERGY in 1000 3000
  do
    ./FAtempsplit.sh ${DIR} ${ENERGY} 2000
  done
