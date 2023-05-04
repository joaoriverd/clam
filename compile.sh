#!/bin/bash
set -e # To exit script if a command fails

# Temporal script to compile all in amd-server machine

# Compile elina
cd /local/home/hectorr/repo/clam/build/crab/elina-prefix/src/elina
make ELINA_PREFIX=/local/home/hectorr/repo/clam/build/run/elina
make install ELINA_PREFIX=/local/home/hectorr/repo/clam/build/run/elina

# Compile Clam and crab
cd /local/home/hectorr/repo/clam/build
make -j20
make install

# Add elina lib binaries to correct directory (not sure why this is not done by make install)
cp -r /local/home/hectorr/repo/clam/install/elina/lib/. /local/home/hectorr/repo/clam/install/lib