"""
    Implementation for P-Tree Programming
    by Christian Oesch, University of Basel, Switzerland
    
    This implementation relies on Cython.
    
    To customize the function set, modify the symbols dictionary. 
        'dnt' contains functions with arity 2, available choices are +, -, *, /, ^ (power)
        'snt' contains functions with arity 1, available choices are s (sin), c (cos), l (log), r (sqrt), e (exp), 
            q (^2), b (^3), t (tan), h (tanh), i (inverse: ^⁻1)
        't' contains terminals, available choices are u, v, w, x, y, z (variables), 1 (constant 1), k (constant) 
"""
from models.config import set
from joblib import Parallel, delayed
import multiprocessing
import benchmarks.symreg as sr
import matplotlib.pyplot as plt
import pyximport
import time
import numpy as np

# Settings for benchmark Multivar:
num_experiments = 100                           # numbers of runs to perform
symbols = {'dnt': ['+', '*'],                   # functions with arity 2
           'snt': ['s', 'c', 'l', 'r', 'i'],    # functions with arity 1
           't': ['k']}                     # terminals (that are not variables)
variables = ['u0', 'v0', 'w0', 'x0', 'y0', 'z0',
             'u1', 'v1', 'w1', 'x1', 'y1', 'z1',
             'u2', 'v2', 'w2', 'x2', 'y2', 'z2',
             'u3', 'v3', 'w3', 'x3', 'y3', 'z3',
             'u4', 'v4', 'w4', 'x4', 'y4', 'z4',
             'u5', 'v5', 'w5', 'x5', 'y5', 'z5',
             'u6', 'v6', 'w6', 'x6']
magnitudes = [-6, -5, -4, -3, -2, -1]   # positions of floating point in a constant
constant_symbols = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'}  # symbols to concatenate an integer
symbols['t'].extend(variables)
settings = {'problem': sr.get_benchmark_multivar,
            'problem_degree': 40,       # difficulty degree of problem, if applicable
            'p_max': 4,                 # exponent for the power law
            'p_min': 1.001,             # branch depth discounting factor
            'delta': 1.00075,           # general discounting factor
            'prune_step': 10000000,     # how many evaluations before the prototype tree is pruned,
                                        # ineffective if > num_evals
            'min_prune_distance': 5,    # how far away from the best solution the branch is trimmed
            'max_depth': 15,            # maximum depth of the prototype tree
            'min_depth': 0,
            'max_const_depth': 3,       # maximum number of digits to concatenate in a constant
            'symbols': symbols,
            'constant_symbols': constant_symbols,
            'variables': variables,
            'magnitudes': magnitudes,
            'num_train_cases': 100,      # number of training cases
            'num_test_cases': 100,      # number of testing cases
            'num_iterations': 1000000}  # number of evaluations to perform

# create compile time definitions
set(settings)

# compile and import ptree
pyximport.install(reload_support=1)
from models.cProbtree import start

num_cores = multiprocessing.cpu_count() - 1  # how many cores to use for the runs
t = time.time()
# results contains the end-of-training, the test errors as well as the error evolution for every 5000th evaluation
results = Parallel(n_jobs=num_cores)(delayed(start)(k, settings['problem'], settings['problem_degree'])
                                         for k in range(num_experiments))
elapsed = time.time()-t
print('Took %f seconds' % (elapsed))
