# compile-time definitions created by config.py: 
# do not change manually 
DEF MAX_DEPTH  = 15 
DEF MIN_DEPTH  = 0 
DEF MAX_CONST_DEPTH  = 3 
DEF MIN_PRUNE_DISTANCE  = 5 
DEF PRUNE_STEP  = 10000000 
DEF NUM_DNT  = 2 
DEF NUM_SNT  = 5 
DEF NUM_SYM  = 9 
DEF NUM_CONST_SYM  = 10 
DEF NUM_MAGNITUDES  = 6 
DEF NUM_ITER  = 1000000 
DEF P_MIN  = 1.001000 
DEF P_MAX  = 4.000000 
DEF DELTA  = 1.000750 
DEF NUM_NT  = NUM_SNT + NUM_DNT 
DEF SNTERMINAL_IDX = NUM_DNT 
DEF TERMINAL_IDX = NUM_NT 
DEF SYM_P = 0.111111 
DEF TERMINAL_P = 0.500000 
DEF NUM_TRAIN_CASES = 49 
DEF NUM_TEST_CASES = 119 
DEF NUM_CASES = 119 
DEF NUM_VARS = 1 
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
DEF SYM_UUU = -1 
DEF SYM_VVV = -1 
DEF SYM_WWW = -1 
DEF SYM_XXX = 7 
DEF SYM_YYY = -1 
DEF SYM_ZZZ = -1 
DEF SYM_CON = 8 
DEF VAR_UUU = -1 
DEF VAR_VVV = -1 
DEF VAR_WWW = -1 
DEF VAR_XXX = 0 
DEF VAR_YYY = -1 
DEF VAR_ZZZ = -1 
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
cdef char symbols[9] 
symbols[:] = ['+', '*', 's', 'c', 'l', 'r', 'i', 'x', 'k'] 
cdef char magnitudes[6] 
magnitudes[:] = [-6, -5, -4, -3, -2, -1] 
cdef double base_probs[9] 
base_probs[:] = [0.9242685895423098, 0.057766786846394365, 0.011410723327682837, 0.003610424177899648, 0.0014788297432676957, 0.0007131702079801773, 0.00038495151584436063, 0.000225651511118728, 0.00014087312750225725]
cdef double double_probs[6] 
double_probs[:] = [0.9249636776898403, 0.05781022985561502, 0.011419304662837535, 0.003613139365975939, 0.0014799418843037446, 0.0007137065414273459]
cdef double const_probs[10] 
const_probs[:] = [0.9241831701947903, 0.0577614481371744, 0.011409668767836917, 0.0036100905085734, 0.0014786930723116646, 0.0007131042979898073, 0.0003849159392731322, 0.0002256306567858375, 0.00014086010824490022, 9.241831701947904e-05]