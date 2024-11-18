#!/bin/bash

diagnostics=(
    Fields
    Particles
    Currents
)

for d in ${diagnostics[*]}; do
    mpirun -np 8 $d
done