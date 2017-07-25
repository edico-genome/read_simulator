#!/bin/bash

python simulator.py --dataset vlrd_chr1_normalNoise_varRate0.002 --varrate 0.002
python simulator.py --dataset vlrd_chr1_normalNoise_varRate0.005 --varrate 0.005

pushd /home/theoh/p4/sw/trunk/test
./test_runner -f suite_def_tmp/mrjd_pirs
popd

python simulator.py --dataset vlrd_normalnoise_errorrate0.002 --varrate 0.002
python simulator.py --dataset vlrd_normalnoise_errorrate0.005 --varrate 0.005


