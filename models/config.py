import copy

sym_def = {'+': 'SYM_PLS',
           '-': 'SYM_MIN',
           '/': 'SYM_DIV',
           '*': 'SYM_MUL',
           '^': 'SYM_POW',
           's': 'SYM_SIN',
           'c': 'SYM_COS',
           'l': 'SYM_LOG',
           'r': 'SYM_SRT',
           'e': 'SYM_EXP',
           'q': 'SYM_QUD',
           't': 'SYM_TAN',
           'h': 'SYM_TAH',
           'b': 'SYM_CUB',
           '1': 'SYM_ONE',
           'i': 'SYM_INV',
           'u0': 'SYM_UU0',
           'v0': 'SYM_VV0',
           'w0': 'SYM_WW0',
           'x0': 'SYM_XX0',
           'y0': 'SYM_YY0',
           'z0': 'SYM_ZZ0',
           'u1': 'SYM_UU1',
           'v1': 'SYM_VV1',
           'w1': 'SYM_WW1',
           'x1': 'SYM_XX1',
           'y1': 'SYM_YY1',
           'z1': 'SYM_ZZ1',
           'u2': 'SYM_UU2',
           'v2': 'SYM_VV2',
           'w2': 'SYM_WW2',
           'x2': 'SYM_XX2',
           'y2': 'SYM_YY2',
           'z2': 'SYM_ZZ2',
           'u3': 'SYM_UU3',
           'v3': 'SYM_VV3',
           'w3': 'SYM_WW3',
           'x3': 'SYM_XX3',
           'y3': 'SYM_YY3',
           'z3': 'SYM_ZZ3',
           'u4': 'SYM_UU4',
           'v4': 'SYM_VV4',
           'w4': 'SYM_WW4',
           'x4': 'SYM_XX4',
           'y4': 'SYM_YY4',
           'z4': 'SYM_ZZ4',
           'u5': 'SYM_UU5',
           'v5': 'SYM_VV5',
           'w5': 'SYM_WW5',
           'x5': 'SYM_XX5',
           'y5': 'SYM_YY5',
           'z5': 'SYM_ZZ5',
           'u6': 'SYM_UU6',
           'v6': 'SYM_VV6',
           'w6': 'SYM_WW6',
           'x6': 'SYM_XX6',
           'y6': 'SYM_YY6',
           'z6': 'SYM_ZZ6',
           'u7': 'SYM_UU7',
           'v7': 'SYM_VV7',
           'w7': 'SYM_WW7',
           'x7': 'SYM_XX7',
           'y7': 'SYM_YY7',
           'z7': 'SYM_ZZ7',
           'u8': 'SYM_UU8',
           'v8': 'SYM_VV8',
           'w8': 'SYM_WW8',
           'x8': 'SYM_XX8',
           'y8': 'SYM_YY8',
           'z8': 'SYM_ZZ8',
           'u9': 'SYM_UU9',
           'v9': 'SYM_VV9',
           'w9': 'SYM_WW9',
           'x9': 'SYM_XX9',
           'y9': 'SYM_YY9',
           'z9': 'SYM_ZZ9',
           'k': 'SYM_CON'}

con_def = {'0': 'CON_0',
           '1': 'CON_1',
           '2': 'CON_2',
           '3': 'CON_3',
           '4': 'CON_4',
           '5': 'CON_5',
           '6': 'CON_6',
           '7': 'CON_7',
           '8': 'CON_8',
           '9': 'CON_9',
           '.': 'CON_P'}

def set(settings):
    with open('models/config.pxi', 'w') as fd:

        symbols = copy.deepcopy(settings['symbols']['dnt'])
        symbols.extend(settings['symbols']['snt'])
        symbols.extend(settings['symbols']['t'])
        num_symbols = len(symbols)
        fd.write('# compile-time definitions created by config.py: \n')
        fd.write('# do not change manually \n')
        fd.write('DEF MAX_DEPTH  = %i \n' % settings['max_depth'])
        fd.write('DEF MIN_DEPTH  = %i \n' % settings['min_depth'])
        fd.write('DEF MAX_CONST_DEPTH  = %i \n' % settings['max_const_depth'])
        fd.write('DEF MIN_PRUNE_DISTANCE  = %i \n' % settings['min_prune_distance'])
        fd.write('DEF PRUNE_STEP  = %i \n' % settings['prune_step'])
        fd.write('DEF NUM_DNT  = %i \n' % len(settings['symbols']['dnt']))
        fd.write('DEF NUM_SNT  = %i \n' % len(settings['symbols']['snt']))
        fd.write('DEF NUM_SYM  = %i \n' % num_symbols)
        fd.write('DEF NUM_CONST_SYM  = %i \n' % len(settings['constant_symbols']))
        fd.write('DEF NUM_MAGNITUDES  = %i \n' % len(settings['magnitudes']))
        fd.write('DEF NUM_ITER  = %i \n' % settings['num_iterations'])
        fd.write('DEF P_MIN  = %f \n' % settings['p_min'])
        fd.write('DEF P_MAX  = %f \n' % settings['p_max'])
        fd.write('DEF DELTA  = %f \n' % settings['delta'])
        fd.write('DEF NUM_NT  = NUM_SNT + NUM_DNT \n')
        fd.write('DEF SNTERMINAL_IDX = NUM_DNT \n')
        fd.write('DEF TERMINAL_IDX = NUM_NT \n')
        fd.write('DEF SYM_P = %f \n' % (1.0 / num_symbols))
        fd.write('DEF TERMINAL_P = %f \n' % (1.0/len(settings['symbols']['t'])))
        fd.write('DEF NUM_TRAIN_CASES = %i \n' % (settings['num_train_cases']))
        fd.write('DEF NUM_TEST_CASES = %i \n' % (settings['num_test_cases']))
        fd.write('DEF NUM_CASES = %i \n' % (max(settings['num_train_cases'], settings['num_test_cases'])))
        fd.write('DEF NUM_VARS = %i \n' % (len(settings['variables'])))
        fd.write('DEF NUM_BATCH_CASES = %i \n' % (int(settings['batch_size'])*settings['num_train_cases']))

        idx = 0
        ordered_symbols = []
        for s in sym_def:
            if s in symbols:
                fd.write('DEF %s = %i \n' % (sym_def[s], idx))
                idx += 1
                ordered_symbols.append(s)
            else:
                fd.write('DEF %s = %i \n' % (sym_def[s], -1))
        idx = 0
        for s in con_def:
            if s in settings['constant_symbols']:
                fd.write('DEF %s = %i \n' % (con_def[s], idx))
                idx += 1
            else:
                fd.write('DEF %s = %i \n' % (con_def[s], -1))
        #fd.write('cdef bytes symbols[%i] \n' % num_symbols)
        fd.write('symbols = %s \n' % str(ordered_symbols))
        fd.write('cdef char magnitudes[%i] \n' % len(settings['magnitudes']))
        fd.write('magnitudes[:] = %s \n' % str(settings['magnitudes']))
        fd.write('cdef double base_probs[%i] \n' % num_symbols)

        # build base probs
        base_probs = []
        p_sum = 0
        for i in range(num_symbols):
            base_probs.append(pow(i+1, -settings['p_max']))
            p_sum += base_probs[i]
        for i in range(num_symbols):
            base_probs[i] = base_probs[i]/p_sum
        fd.write('base_probs[:] = %s\n' % str(base_probs))
        fd.write('cdef double double_probs[%i] \n' % len(settings['magnitudes']))

        # build double probs
        double_probs = []
        p_sum = 0
        for i in range(len(settings['magnitudes'])):
            double_probs.append(pow(i+1, -settings['p_max']))
            p_sum += double_probs[i]
        for i in range(len(settings['magnitudes'])):
            double_probs[i] = double_probs[i]/p_sum
        fd.write('double_probs[:] = %s\n' % str(double_probs))

        fd.write('cdef double const_probs[%i] \n' % len(settings['constant_symbols']))

        # build const probs
        const_probs = []
        p_sum = 0
        for i in range(len(settings['constant_symbols'])):
            const_probs.append(pow(i + 1, -settings['p_max']))
            p_sum += const_probs[i]
        for i in range(len(settings['constant_symbols'])):
            const_probs[i] = const_probs[i] / p_sum
        fd.write('const_probs[:] = %s\n' % str(const_probs))
