#!/bin/bash

rm *0
rm *.nc
./microhh init dycoms
./microhh run dycoms >& log
#python 3d_to_nc.py
