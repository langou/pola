#!/bin/sh
./main_a2v__time.exe -m 19 -n 13 -method hh_a2v_rec_blas
./main_a2v__time.exe -m 19 -n 13 -method hh_a2v_rl
./main_a2v__time.exe -m 19 -n 13 -method hh_a2v_rec_blas
./main_a2v__time.exe -m 19 -n 13 -method hh_a2v_rl_blas
./main_a2v__time.exe -m 19 -n 13 -method hh_a2v_ll
./main_a2v__time.exe -m 19 -n 13 -method hh_a2v_ll_blas
./main_a2v__time.exe -m 19 -n 13 -method hh_a2v_ll_tiled
./main_a2v__time.exe -m 19 -n 13 -method hh_a2v_ll_tiled_blas
./main_a2v__time.exe -m 19 -n 13 -method hh_a2v_rl_tiled
./main_a2v__time.exe -m 19 -n 13 -method hh_a2v_rl_tiled_blas
./main_a2v__time.exe -m 19 -n 13 -method geqr2
./main_a2v__time.exe -m 19 -n 13 -method geqrf
