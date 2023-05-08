#!/bin/sh
ver=5.1.0
link=http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-$ver.tar.gz
curl $link --output metis-$ver.tar.gz
gunzip metis-$ver.tar.gz
tar -xvf metis-$ver.tar
cd metis-$ver
make config
make
sudo make install
cd ..
rm -rf metis-$ver metis-$ver.tar.gz metis-$ver.tar
