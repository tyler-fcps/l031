#!/bin/bash
# Compiles and runs all .cpp files
g++ -Wall -g --std c++11 -fdiagnostics-color=always ./**.cpp -O3 -o a.out
bash run.sh