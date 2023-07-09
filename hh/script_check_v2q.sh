#!/bin/sh
./main_v2q__time.exe -m 19 -n 13 -method hh_v2q_rl
./main_v2q__time.exe -m 19 -n 13 -method hh_v2q_rl_blas
./main_v2q__time.exe -m 19 -n 13 -method hh_v2q_rec_blas
./main_v2q__time.exe -m 19 -n 13 -method hh_v2q_ll
./main_v2q__time.exe -m 19 -n 13 -method hh_v2q_ll_blas
./main_v2q__time.exe -m 19 -n 13 -method hh_v2q_ll_tiled
./main_v2q__time.exe -m 19 -n 13 -method hh_v2q_ll_tiled_blas
./main_v2q__time.exe -m 19 -n 13 -method hh_v2q_rl_tiled
./main_v2q__time.exe -m 19 -n 13 -method hh_v2q_rl_tiled_blas
./main_v2q__time.exe -m 19 -n 13 -method org2r
./main_v2q__time.exe -m 19 -n 13 -method orgqr








