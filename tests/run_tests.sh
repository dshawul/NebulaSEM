#!/bin/bash

set -aex

allTests=( \
    backface/backface \
    cavity/cavity \
    cavity-amr/cavity \
    convec/convec \
    harg/harg \
    laplace/laplace \
    laplace-amr/laplace \
    mix/mix \
    pitzdaily/ke/pitzDaily \
    pitzdaily/les/pitzDaily \
    potential/potential \
    transport/wave2d/simple \
    transport/wave2d-amr/simple \
    transport/wave2d-amr-dg/simple \
    atmo/acoustic-sphere/sphere \
    atmo/acoustic-sphere-amr/sphere \
    atmo/acoustic-sphere-amr-dg/sphere \
    atmo/advection-sphere/sphere \
    atmo/ctbs/bubble \
    atmo/dc/bubble \
    atmo/filament/bubble \
    atmo/hydro-sphere/sphere \
    atmo/lrtb/bubble \
    atmo/srtb/bubble \
    atmo/srtb-amr/bubble \
    atmo/srtb-amr-zaxis/bubble \
    atmo/srtb-curved/bubble \
    atmo/srtb-inclined/bubble \
    atmo/srtb-piso/bubble \
    atmo/srtb-piso-amr/bubble \
)

cd ..

for i in "${allTests[@]}"; do
    echo "Running test case: $i"
    ./test.sh -s 10 -n 1 -c examples/$i
done
