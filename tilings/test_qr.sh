#!/bin/sh
./main_qr__check_qr_mgs_rl__tiled.exe
./main_qr__check_qr_mgs_rl__tiled.exe -method cgs_rl       -m 10 -n  8
./main_qr__check_qr_mgs_rl__tiled.exe -method cgs_ll       -m 10 -n  8
./main_qr__check_qr_mgs_rl__tiled.exe -method mgs_pb       -m 10 -n  8
./main_qr__check_qr_mgs_rl__tiled.exe -method mgs_rl       -m 10 -n  8
./main_qr__check_qr_mgs_rl__tiled.exe -method mgs_ll       -m 10 -n  8
./main_qr__check_qr_mgs_rl__tiled.exe -method cgs2_ll      -m 10 -n  8
./main_qr__check_qr_mgs_rl__tiled.exe -method hh_a2v_v2q   -m 10 -n  8
./main_qr__check_qr_mgs_rl__tiled.exe -method hh_a2q       -m 10 -n  8
./main_qr__check_qr_mgs_rl__tiled.exe -method cgs_rl       -m 10 -n 10
./main_qr__check_qr_mgs_rl__tiled.exe -method cgs_ll       -m 10 -n 10
./main_qr__check_qr_mgs_rl__tiled.exe -method mgs_pb       -m 10 -n 10
./main_qr__check_qr_mgs_rl__tiled.exe -method mgs_rl       -m 10 -n 10
./main_qr__check_qr_mgs_rl__tiled.exe -method mgs_ll       -m 10 -n 10
./main_qr__check_qr_mgs_rl__tiled.exe -method cgs2_ll      -m 10 -n 10
./main_qr__check_qr_mgs_rl__tiled.exe -method hh_a2v_v2q   -m 10 -n 10
./main_qr__check_qr_mgs_rl__tiled.exe -method hh_a2q       -m 10 -n 10
