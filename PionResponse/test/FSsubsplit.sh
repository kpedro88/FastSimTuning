#!/bin/bash

DIR=pion

for ENERGY in 1 2 3 5 9 10 11 15 20 30
  do
    for PART in {1..50}
	  do
	    ./FStempEsplit.sh ${DIR} ${ENERGY} 10000 ${PART}
	  done
  done

for ENERGY in 50 100 150 225 300
  do
    for PART in {1..100}
	  do
	    ./FStempEsplit.sh ${DIR} ${ENERGY} 5000 ${PART}
	  done
  done

for ENERGY in 1000 3000
  do
    for PART in {1..250}
	  do
	    ./FStempEsplit.sh ${DIR} ${ENERGY} 2000 ${PART}
	  done
  done
