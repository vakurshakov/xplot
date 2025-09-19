#!/bin/bash

if [[ $# == 0 ]]; then
    echo "Directory is needed"
    exit
fi

params_path=$1
shift 1

if [[ $# == 0 ]]; then
    diagnostics=(
        Fields
        Currents
        Info_Ions
        Info_Electrons
    )
else
    diagnostics=$@
fi

for d in ${diagnostics[*]}; do
    pushd ./$params_path/$d

    # i=0
    # for f in *; do
    #     mv $f $(printf %04d $i)".png"
    #     i=$((i + 1))
    # done

    # TODO: Investigate size and quality optimizations for a large image-plots!
    ffmpeg -y -i ./%04d.png -r 15 ../Video/$d.mp4

    popd
done
