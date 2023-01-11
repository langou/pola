#!/bin/sh
./main_qr.exe
./main_qr.exe -method cgs_rl       -m 10 -n  8
./main_qr.exe -method cgs_ll       -m 10 -n  8
./main_qr.exe -method mgs_pb       -m 10 -n  8
./main_qr.exe -method mgs_rl       -m 10 -n  8
./main_qr.exe -method mgs_ll       -m 10 -n  8
./main_qr.exe -method cgs2_ll      -m 10 -n  8
./main_qr.exe -method hh_a2v_v2q   -m 10 -n  8
./main_qr.exe -method hh_a2q       -m 10 -n  8
./main_qr.exe -method cgs_rl       -m 10 -n 10
./main_qr.exe -method cgs_ll       -m 10 -n 10
./main_qr.exe -method mgs_pb       -m 10 -n 10
./main_qr.exe -method mgs_rl       -m 10 -n 10
./main_qr.exe -method mgs_ll       -m 10 -n 10
./main_qr.exe -method cgs2_ll      -m 10 -n 10
./main_qr.exe -method hh_a2v_v2q   -m 10 -n 10
./main_qr.exe -method hh_a2q       -m 10 -n 10
