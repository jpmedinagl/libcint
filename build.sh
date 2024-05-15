#!/bin/bash

rm -rf build

mkdir build
cd build

ENZYME=/u/jpmedina/Enzyme/enzyme/build

cmake -G Ninja -DENABLE_EXAMPLE=1 -DENABLE_STATIC=1 -DENABLE_TEST=1 -DQUICK_TEST=1 -DEnzyme_DIR=$ENZYME ..
ninja
