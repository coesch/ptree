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

# Settings for benchmark Keijzer-6:
num_experiments = 1                           # numbers of runs to perform
symbols = {'dnt': ['+', '*'],                   # functions with arity 2
           'snt': ['s', 'c', 'l', 'r', 'i'],    # functions with arity 1
           't': ['x0', 'k']}                     # terminals
variables = ['x']
magnitudes = [-6, -5, -4, -3, -2, -1]   # positions of floating point in a constant
constant_symbols = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'}  # symbols to concatenate an integer
settings = {'p_max': 4,                 # exponent for the power law
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
            'num_train_cases': 49,      # number of training cases
            'num_test_cases': 119,      # number of testing cases
            'num_iterations': 10000}  # number of evaluations to perform

# create compile time definitions
set(settings)

# compile and import ptree
pyximport.install(reload_support=1)
from models.cProbtree import start

num_cores = multiprocessing.cpu_count() - 1  # how many cores to use for the runs
t = time.time()
# results contains the end-of-training, the test errors as well as the error evolution for every 5000th evaluation
results = Parallel(n_jobs=num_cores)(delayed(start)(k, sr.get_benchmark_keijzer, 6)
                                         for k in range(num_experiments))
elapsed = time.time()-t
print('Took %f seconds' % (elapsed))
