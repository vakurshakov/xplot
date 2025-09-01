#!/bin/bash

diagnostics=(
    fields
    currents
    particles
)

for d in ${diagnostics[*]}; do
    mpiexec -np 8 ./plot/$d.py
done