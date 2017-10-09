# compile-time definitions created by config.py: 
# do not change manually 
DEF MAX_DEPTH  = 15 
DEF MIN_DEPTH  = 0 
DEF MAX_CONST_DEPTH  = 3 
DEF MIN_PRUNE_DISTANCE  = 5 
DEF PRUNE_STEP  = 10000000 
DEF NUM_DNT  = 2 
DEF NUM_SNT  = 5 
DEF NUM_SYM  = 48 
DEF NUM_CONST_SYM  = 10 
DEF NUM_MAGNITUDES  = 6 
DEF NUM_ITER  = 1000000 
DEF P_MIN  = 1.001000 
DEF P_MAX  = 4.000000 
DEF DELTA  = 1.000750 
DEF NUM_NT  = NUM_SNT + NUM_DNT 
DEF SNTERMINAL_IDX = NUM_DNT 
DEF TERMINAL_IDX = NUM_NT 
DEF SYM_P = 0.020833 
DEF TERMINAL_P = 0.024390 
DEF NUM_TRAIN_CASES = 100 
DEF NUM_TEST_CASES = 100 
DEF NUM_CASES = 100 
DEF NUM_VARS = 40 
DEF SYM_PLS = 0 
DEF SYM_MIN = -1 
DEF SYM_DIV = -1 
DEF SYM_MUL = 1 
DEF SYM_POW = -1 
DEF SYM_SIN = 2 
DEF SYM_COS = 3 
DEF SYM_LOG = 4 
DEF SYM_SRT = 5 
DEF SYM_EXP = -1 
DEF SYM_QUD = -1 
DEF SYM_TAN = -1 
DEF SYM_TAH = -1 
DEF SYM_CUB = -1 
DEF SYM_ONE = -1 
DEF SYM_INV = 6 
DEF SYM_UU0 = 7 
DEF SYM_VV0 = 8 
DEF SYM_WW0 = 9 
DEF SYM_XX0 = 10 
DEF SYM_YY0 = 11 
DEF SYM_ZZ0 = 12 
DEF SYM_UU1 = 13 
DEF SYM_VV1 = 14 
DEF SYM_WW1 = 15 
DEF SYM_XX1 = 16 
DEF SYM_YY1 = 17 
DEF SYM_ZZ1 = 18 
DEF SYM_UU2 = 19 
DEF SYM_VV2 = 20 
DEF SYM_WW2 = 21 
DEF SYM_XX2 = 22 
DEF SYM_YY2 = 23 
DEF SYM_ZZ2 = 24 
DEF SYM_UU3 = 25 
DEF SYM_VV3 = 26 
DEF SYM_WW3 = 27 
DEF SYM_XX3 = 28 
DEF SYM_YY3 = 29 
DEF SYM_ZZ3 = 30 
DEF SYM_UU4 = 31 
DEF SYM_VV4 = 32 
DEF SYM_WW4 = 33 
DEF SYM_XX4 = 34 
DEF SYM_YY4 = 35 
DEF SYM_ZZ4 = 36 
DEF SYM_UU5 = 37 
DEF SYM_VV5 = 38 
DEF SYM_WW5 = 39 
DEF SYM_XX5 = 40 
DEF SYM_YY5 = 41 
DEF SYM_ZZ5 = 42 
DEF SYM_UU6 = 43 
DEF SYM_VV6 = 44 
DEF SYM_WW6 = 45 
DEF SYM_XX6 = 46 
DEF SYM_YY6 = -1 
DEF SYM_ZZ6 = -1 
DEF SYM_UU7 = -1 
DEF SYM_VV7 = -1 
DEF SYM_WW7 = -1 
DEF SYM_XX7 = -1 
DEF SYM_YY7 = -1 
DEF SYM_ZZ7 = -1 
DEF SYM_UU8 = -1 
DEF SYM_VV8 = -1 
DEF SYM_WW8 = -1 
DEF SYM_XX8 = -1 
DEF SYM_YY8 = -1 
DEF SYM_ZZ8 = -1 
DEF SYM_UU9 = -1 
DEF SYM_VV9 = -1 
DEF SYM_WW9 = -1 
DEF SYM_XX9 = -1 
DEF SYM_YY9 = -1 
DEF SYM_ZZ9 = -1 
DEF SYM_CON = 47 
DEF CON_0 = 0 
DEF CON_1 = 1 
DEF CON_2 = 2 
DEF CON_3 = 3 
DEF CON_4 = 4 
DEF CON_5 = 5 
DEF CON_6 = 6 
DEF CON_7 = 7 
DEF CON_8 = 8 
DEF CON_9 = 9 
DEF CON_P = -1 
symbols = ['+', '*', 's', 'c', 'l', 'r', 'i', 'u0', 'v0', 'w0', 'x0', 'y0', 'z0', 'u1', 'v1', 'w1', 'x1', 'y1', 'z1', 'u2', 'v2', 'w2', 'x2', 'y2', 'z2', 'u3', 'v3', 'w3', 'x3', 'y3', 'z3', 'u4', 'v4', 'w4', 'x4', 'y4', 'z4', 'u5', 'v5', 'w5', 'x5', 'y5', 'z5', 'u6', 'v6', 'w6', 'x6', 'k'] 
cdef char magnitudes[6] 
magnitudes[:] = [-6, -5, -4, -3, -2, -1] 
cdef double base_probs[48] 
base_probs[:] = [0.9239408966459561, 0.05774630604037226, 0.011406677736369828, 0.003609144127523266, 0.00147830543463353, 0.0007129173585231143, 0.00038481503400497965, 0.00022557150797020413, 0.0001408231819304917, 9.239408966459562e-05, 6.310640643712562e-05, 4.455733490769464e-05, 3.23497390373571e-05, 2.405093962531123e-05, 1.8250684378191728e-05, 1.4098219248137758e-05, 1.106237828385623e-05, 8.801448870655732e-06, 7.089731483383002e-06, 5.7746306040372265e-06, 4.750802888950366e-06, 3.944150402320351e-06, 3.3016637899591415e-06, 2.784833431730915e-06, 2.365288695413648e-06, 2.0218586898348188e-06, 1.738557801611009e-06, 1.5031837265819518e-06, 1.306327890394279e-06, 1.140667773636983e-06, 1.000454669299297e-06, 8.811387030086099e-07, 7.79091437495378e-07, 6.913986427410144e-07, 6.157040544079675e-07, 5.500905544159832e-07, 4.929890743889965e-07, 4.4310821771143763e-07, 3.993794942883593e-07, 3.6091441275232665e-07, 3.269706449504952e-07, 2.9692518055939785e-07, 2.702529034728714e-07, 2.4650940014502194e-07, 2.2531709108878675e-07, 2.0635398687244634e-07, 1.8934452818656715e-07, 1.740520894831822e-07]
cdef double double_probs[6] 
double_probs[:] = [0.9249636776898403, 0.05781022985561502, 0.011419304662837535, 0.003613139365975939, 0.0014799418843037446, 0.0007137065414273459]
cdef double const_probs[10] 
const_probs[:] = [0.9241831701947903, 0.0577614481371744, 0.011409668767836917, 0.0036100905085734, 0.0014786930723116646, 0.0007131042979898073, 0.0003849159392731322, 0.0002256306567858375, 0.00014086010824490022, 9.241831701947904e-05]
