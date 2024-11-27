#!/bin/bash

set -v
llvm-extract -S --func=cint_wrap --recursive --rfunc="enzyme_opt_helper_*" out.ll -o mwe.ll