#!/bin/bash

rm -rf build

mkdir build
cd build

ENZYME=/u/jpmedina/Enzyme/enzyme/build

cmake -G Ninja -DENABLE_STATIC=1 -DEnzyme_DIR=$ENZYME -DSO=0 -DGRAD=0 -DSCF=1 ..
ninja
