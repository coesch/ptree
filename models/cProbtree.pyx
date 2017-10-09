# cython: profile=False
import random
import numpy as np

from libc.stdlib cimport malloc, free#,
from libc.stdio cimport sscanf
from libc.math cimport sin, cos, log, exp, fabs, pow, sqrt, tan, tanh, isnan
cdef extern from "stdlib.h":
    ctypedef void const_void "const void"
    void qsort(void *base, int nmemb, int size,
                int(*compar)(const_void *, const_void *)) nogil
    double drand48()
    void srand48(long int seedval)

# compile-time definitions:
include "config.pxi"

# base node type
cdef struct Node:
    unsigned int depth
    unsigned int num_evals
    unsigned int current_choice
    unsigned int best
    Node* nodes_left[NUM_NT]
    Node* nodes_right[NUM_DNT]
    FloatNode* cnode
    double probs[NUM_SYM]
    double fitness[NUM_SYM]


cdef struct ConstNode:
    unsigned int depth
    unsigned int num_evals
    unsigned int current_choice
    unsigned int best
    ConstNode* nodes[NUM_CONST_SYM]
    double probs[NUM_CONST_SYM]
    double fitness[NUM_CONST_SYM]


cdef struct FloatNode:
    unsigned int depth
    unsigned int num_evals
    unsigned int current_choice
    unsigned int best
    ConstNode* nodes[NUM_MAGNITUDES]
    double probs[NUM_MAGNITUDES]
    double fitness[NUM_MAGNITUDES]


cdef create_node(Node *self, unsigned int depth):
    self.depth = depth
    self.num_evals = 0
    self.current_choice = SYM_XXX

    # first evaluation will choose a terminal
    cdef int k
    for k in range(NUM_SYM):
        self.probs[k] = 0.0
        self.fitness[k] = 1E120

    for k in range(TERMINAL_IDX,NUM_SYM):
        self.probs[k] = <double> TERMINAL_P
    for k in range(NUM_NT):
        self.nodes_left[k] = NULL
    for k in range(NUM_DNT):
        self.nodes_right[k] = NULL
    self.cnode = NULL


cdef create_const_node(ConstNode *self, unsigned int depth):
    self.depth = depth
    self.current_choice = 0
    self.num_evals = 0
    cdef int k
    for k in range(NUM_CONST_SYM):
        self.probs[k] = 1.0/<double>NUM_CONST_SYM
        self.fitness[k] = 1E120
        self.nodes[k] = NULL


cdef create_double_node(FloatNode *self, unsigned int depth):
    self.depth = depth
    self.current_choice = 0
    self.num_evals = 0
    cdef int k
    for k in range(NUM_MAGNITUDES):
        self.probs[k] = 1.0/<double>NUM_MAGNITUDES
        self.fitness[k] = 1E120
        self.nodes[k] = NULL


cdef create_tree(Node* self):
    self.num_evals += 1
    cdef double r = random.random()
    cdef double p_sum = 0
    cdef unsigned int choice
    # choose any symbol
    for choice in range(NUM_SYM):
        p_sum += self.probs[choice]
        if r < p_sum:
            break
    self.current_choice = choice
    if choice < SNTERMINAL_IDX: # its a double non terminal
        if self.nodes_left[choice] == NULL:
            self.nodes_left[choice] = <Node *>malloc(sizeof(Node))
            self.nodes_right[choice] = <Node *>malloc(sizeof(Node))
            create_node(self.nodes_left[choice], self.depth+1)
            create_node(self.nodes_right[choice], self.depth+1)
        create_tree(self.nodes_left[choice])
        create_tree(self.nodes_right[choice])
    elif choice < TERMINAL_IDX: # its a single non terminal
        if self.nodes_left[choice] == NULL:
            self.nodes_left[choice] = <Node *>malloc(sizeof(Node))
            create_node(self.nodes_left[choice], self.depth+1)
        create_tree(self.nodes_left[choice])
    elif choice == SYM_CON:
        if self.cnode == NULL:
            self.cnode = <FloatNode *>malloc(sizeof(FloatNode))
            create_double_node(self.cnode, 0)
        create_double_tree(self.cnode)


cdef create_const_tree(ConstNode* self):
    self.num_evals += 1
    cdef double r = random.random()
    cdef double p_sum = 0
    cdef unsigned int choice
    # choose any symbol
    for choice in range(NUM_CONST_SYM):
        p_sum += self.probs[choice]
        if r < p_sum:
            #self.current_choice = choice
            break
    self.current_choice = choice
    if self.depth < MAX_CONST_DEPTH-1:
        if self.nodes[choice] == NULL:
            self.nodes[choice] = <ConstNode *>malloc(sizeof(ConstNode))
            create_const_node(self.nodes[choice], self.depth+1)

        create_const_tree(self.nodes[choice]) # do not allow any further periods


cdef create_double_tree(FloatNode* self):
    self.num_evals += 1
    cdef double r = random.random()
    cdef double p_sum = 0
    cdef unsigned int choice
    # choose any symbol
    for choice in range(NUM_MAGNITUDES):
        p_sum += self.probs[choice]
        if r < p_sum:
            break
    self.current_choice = choice
    if self.nodes[choice] == NULL:
        self.nodes[choice] = <ConstNode *>malloc(sizeof(ConstNode))
        create_const_node(self.nodes[choice], 0)
    create_const_tree(self.nodes[choice]) # do not allow any further periods

cdef get_result(Node* self, const int num_cases, const double[][NUM_VARS] model_input, double[] result):
    cdef double temp1, temp2
    cdef unsigned int choice = self.current_choice
    if choice == SYM_ONE:
        for i in range(num_cases):
            result[i] = 1.0
        return
    elif choice == SYM_UUU:
        for i in range(num_cases):
            result[i] = model_input[i][VAR_UUU]
        return
    elif choice == SYM_VVV:
        for i in range(num_cases):
            result[i] = model_input[i][VAR_VVV]
        return
    elif choice == SYM_WWW:
        for i in range(num_cases):
            result[i] = model_input[i][VAR_WWW]
        return
    elif choice == SYM_XXX:
        for i in range(num_cases):
            result[i] = model_input[i][VAR_XXX]
        return
    elif choice == SYM_YYY:
        for i in range(num_cases):
            result[i] = model_input[i][VAR_YYY]
        return
    elif choice == SYM_ZZZ:
        for i in range(num_cases):
            result[i] = model_input[i][VAR_ZZZ]
        return

    cdef double res1[NUM_CASES]
    cdef double res2[NUM_CASES]
    if choice == SYM_PLS:
        get_result(self.nodes_left[choice], num_cases, model_input, res1)
        get_result(self.nodes_right[choice], num_cases, model_input, res2)
        for i in range(num_cases):
            result[i] = res1[i]+res2[i]
            #print('%f + %f = %f' %(res1[i], res2[i], result[i]))
        return
    if choice == SYM_MIN:
        get_result(self.nodes_left[choice], num_cases, model_input, res1)
        get_result(self.nodes_right[choice], num_cases, model_input, res2)
        for i in range(num_cases):
            result[i] = res1[i]-res2[i]
        return
    elif choice == SYM_DIV:
        get_result(self.nodes_left[choice], num_cases, model_input, res1)
        get_result(self.nodes_right[choice], num_cases, model_input, res2)
        for i in range(num_cases):
            if fabs(res2[i]) > 1E-32:
                result[i] = res1[i]/res2[i]
            else:
                result[i] = 1.0
        return
    elif choice == SYM_MUL:
        get_result(self.nodes_left[choice], num_cases, model_input, res1)
        get_result(self.nodes_right[choice], num_cases, model_input, res2)
        for i in range(num_cases):
            result[i] = res1[i]*res2[i]
        return
    elif choice == SYM_POW:
        get_result(self.nodes_left[choice], num_cases, model_input, res1)
        get_result(self.nodes_right[choice], num_cases, model_input, res2)
        for i in range(num_cases):
            result[i] = pow(res1[i], res2[i])
        return
    elif choice == SYM_SIN:
        get_result(self.nodes_left[choice], num_cases, model_input, res1)
        for i in range(num_cases):
            result[i] = sin(res1[i])
        return
    elif choice == SYM_COS:
        get_result(self.nodes_left[choice], num_cases, model_input, res1)
        for i in range(num_cases):
            result[i] = cos(res1[i])
        return
    elif choice == SYM_SRT:
        get_result(self.nodes_left[choice], num_cases, model_input, res1)
        for i in range(num_cases):
            result[i] = sqrt(fabs(res1[i]))
        return
    elif choice == SYM_LOG:
        get_result(self.nodes_left[choice], num_cases, model_input, res1)
        for i in range(num_cases):
            res1[i] = fabs(res1[i])
            if res1[i] > 1E-32:
                result[i] = log(res1[i])
            else:
                result[i] = 1.0
        return
    elif choice == SYM_EXP:
        get_result(self.nodes_left[choice], num_cases, model_input, res1)
        for i in range(num_cases):
            if res1[i] < 500:
                result[i] = exp(res1[i])
            else:
                result[i] = 1.0
        return
    elif choice == SYM_QUD:
        get_result(self.nodes_left[choice], num_cases, model_input, res1)
        for i in range(num_cases):
            result[i] = pow(res1[i],2)
        return
    elif choice == SYM_CUB:
        get_result(self.nodes_left[choice], num_cases, model_input, res1)
        for i in range(num_cases):
            result[i] = pow(res1[i],3)
        return

    elif choice == SYM_TAN:
        get_result(self.nodes_left[choice], num_cases, model_input, res1)
        for i in range(num_cases):
            result[i] = tan(res1[i])
        return
    elif choice == SYM_TAH:
        get_result(self.nodes_left[choice], num_cases, model_input, res1)
        for i in range(num_cases):
            result[i] = tanh(res1[i])
        return
    elif choice == SYM_CON:
        res1[0] = get_double_result(self.cnode)
        for i in range(num_cases):
            result[i] = res1[0]
        return
    elif choice == SYM_INV:
        get_result(self.nodes_left[choice], num_cases, model_input, res1)
        for i in range(num_cases):
            if fabs(res1[i]) > 1E-32:
                result[i] = 1.0/res1[i]
            else:
                result[i] = 1.0
        return
    else:
        print('Symbol not found: %i' % choice)
        return 0


cdef double get_const_result(ConstNode* self, char[] result):
    cdef unsigned int choice = self.current_choice
    cdef double temp1
    cdef char magnitude
    if self.nodes[choice] != NULL:
        get_const_result(self.nodes[choice], result)
    if choice == CON_0:
        result[self.depth] = '0'
    elif choice == CON_1:
        result[self.depth] = '1'
    elif choice == CON_2:
        result[self.depth] = '2'
    elif choice == CON_3:
        result[self.depth] = '3'
    elif choice == CON_4:
        result[self.depth] = '4'
    elif choice == CON_5:
        result[self.depth] = '5'
    elif choice == CON_6:
        result[self.depth] = '6'
    elif choice == CON_7:
        result[self.depth] = '7'
    elif choice == CON_8:
        result[self.depth] = '8'
    elif choice == CON_9:
        result[self.depth] = '9'
    else:
        print('Const symbol not found: %i' % choice)


cdef double get_double_result(FloatNode* self):
    cdef char const_result[MAX_CONST_DEPTH]
    cdef unsigned int choice = self.current_choice
    cdef float temp1
    cdef char magnitude
    get_const_result(self.nodes[choice], const_result)
    sscanf(const_result, '%f', &temp1);
    temp1 = temp1*pow(10, magnitudes[choice])
    return temp1

cdef unsigned int propagate_fitness(Node* self, const double fitness):
    cdef unsigned int depth1 = self.depth
    cdef unsigned int depth2 = self.depth
    if self.current_choice < SNTERMINAL_IDX:
        depth1 = propagate_fitness(self.nodes_left[self.current_choice], fitness)
        depth2 = propagate_fitness(self.nodes_right[self.current_choice], fitness)
        depth1 = max(depth1, depth2)
    elif self.current_choice < TERMINAL_IDX:
        depth1 = propagate_fitness(self.nodes_left[self.current_choice], fitness)
    elif self.current_choice == SYM_CON:
        # do not take depth of constant into account
        propagate_double_fitness(self.cnode, fitness)
    if fitness < self.fitness[self.current_choice]:
        self.fitness[self.current_choice] = fitness*pow(P_MIN, depth1 - self.depth)
    else:
        self.fitness[self.current_choice] *= DELTA
    # calculate probabilities
    calc_probs(self.num_evals, self.fitness, self.probs, self.depth)
    return depth1

cdef unsigned int propagate_const_fitness(ConstNode* self, const double fitness):
    cdef unsigned int depth = self.depth
    if self.nodes[self.current_choice] != NULL:
        depth = propagate_const_fitness(self.nodes[self.current_choice], fitness)
    if fitness < self.fitness[self.current_choice]:
        self.fitness[self.current_choice] = fitness*pow(P_MIN, depth - self.depth)
    else:
        self.fitness[self.current_choice] *= DELTA
    # calculate probabilities
    calc_const_probs(self.num_evals, self.fitness, self.probs, self.depth)
    return depth

cdef unsigned int propagate_double_fitness(FloatNode* self, const double fitness):
    cdef unsigned int depth = self.depth
    if self.nodes[self.current_choice] != NULL:
        depth = propagate_const_fitness(self.nodes[self.current_choice], fitness)
    if fitness < self.fitness[self.current_choice]:
        self.fitness[self.current_choice] = fitness
    else:
        self.fitness[self.current_choice] *= DELTA
    # calculate probabilities
    calc_double_probs(self.num_evals, self.fitness, self.probs, self.depth)
    return depth


cdef unsigned int set_best(Node* self):
    if self.current_choice < SNTERMINAL_IDX:
        set_best(self.nodes_left[self.current_choice])
        set_best(self.nodes_right[self.current_choice])
    elif self.current_choice < TERMINAL_IDX:
        set_best(self.nodes_left[self.current_choice])
    elif self.current_choice == SYM_CON:
        set_double_best(self.cnode)
    self.best = self.current_choice


cdef unsigned int set_const_best(ConstNode* self):
    if self.nodes[self.current_choice] != NULL:
        set_const_best(self.nodes[self.current_choice])
    self.best = self.current_choice


cdef unsigned int set_double_best(FloatNode* self):
    if self.nodes[self.current_choice] != NULL:
        set_const_best(self.nodes[self.current_choice])
    self.best = self.current_choice


cdef unsigned int set_current_to_best(Node* self):
    if self.best < SNTERMINAL_IDX:
        set_current_to_best(self.nodes_left[self.best])
        set_current_to_best(self.nodes_right[self.best])
    elif self.best < TERMINAL_IDX:
        set_current_to_best(self.nodes_left[self.best])
    elif self.best == SYM_CON:
        set_double_current_to_best(self.cnode)
    self.current_choice = self.best


cdef unsigned int set_const_current_to_best(ConstNode* self):
    if self.nodes[self.best] != NULL:
        set_const_current_to_best(self.nodes[self.best])
    self.current_choice = self.best


cdef unsigned int set_double_current_to_best(FloatNode* self):
    if self.nodes[self.best] != NULL:
        set_const_current_to_best(self.nodes[self.best])
    self.current_choice = self.best


cdef unsigned int get_depth(Node* self):
    cdef unsigned int depth1 = self.depth
    cdef unsigned int depth2 = self.depth
    if self.current_choice < SNTERMINAL_IDX:
        depth1 = get_depth(self.nodes_left[self.current_choice])
        depth2 = get_depth(self.nodes_right[self.current_choice])
        depth1 = max(depth1, depth2)
    elif self.current_choice < TERMINAL_IDX:
        depth1 = get_depth(self.nodes_left[self.current_choice])
    return depth1


cdef calc_probs(int num_evals, const double[] fitness, double[] probs, int depth):
    cdef unsigned int i, j
    if depth > MAX_DEPTH-1:
        for i in range(TERMINAL_IDX, NUM_SYM):
            if fitness[i]>1E119:
                for j in range(TERMINAL_IDX, NUM_SYM):
                    probs[j] = <double> TERMINAL_P
                return
    else:
        for i in range(NUM_SYM):
            if fitness[i]>1E119:
                for j in range(NUM_SYM):
                    probs[j] = SYM_P
                return
    cdef double sorted_fit[NUM_SYM][2]
    cdef double p_sum = 0
    for i in range(NUM_SYM):
        sorted_fit[i][0] = fitness[i]
        sorted_fit[i][1] = <double>i
    qsort(sorted_fit, NUM_SYM, sizeof(sorted_fit[0]), compare)
    # unevaluated will have fitness 1E120
    for i in range(NUM_SYM):
        if depth < MAX_DEPTH:
            probs[<int>sorted_fit[i][1]] = base_probs[i]
        elif <int>sorted_fit[i][1] >= TERMINAL_IDX:
            probs[<int>sorted_fit[i][1]] = <double> TERMINAL_P
        else:
            probs[<int>sorted_fit[i][1]] = 0.0  # if we are at max depth set all non-terminals probs to 0


cdef calc_const_probs(int num_evals, const double[] fitness, double[] probs, int depth):
    cdef unsigned int i, j
    for i in range(NUM_CONST_SYM):
        if fitness[i]>1E119:
            # use starting probs
            return
    cdef double sorted_fit[NUM_CONST_SYM][2]
    cdef double p_sum = 0
    for i in range(NUM_CONST_SYM):
        sorted_fit[i][0] = fitness[i]
        sorted_fit[i][1] = <double>i
    qsort(sorted_fit, NUM_CONST_SYM, sizeof(sorted_fit[0]), compare)
    for i in range(NUM_CONST_SYM):
        probs[<int>sorted_fit[i][1]] = const_probs[i] # pow(i+1, -P_MAX)


cdef calc_double_probs(int num_evals, const double[] fitness, double[] probs, int depth):
    cdef unsigned int i, j
    for i in range(NUM_MAGNITUDES):
        if fitness[i]>1E119:
            # use starting probs
            return
    cdef double sorted_fit[NUM_MAGNITUDES][2]
    cdef double p_sum = 0
    for i in range(NUM_MAGNITUDES):
        sorted_fit[i][0] = fitness[i]
        sorted_fit[i][1] = <double>i
    qsort(sorted_fit, NUM_MAGNITUDES, sizeof(sorted_fit[0]), compare)
    # unevaluated will have fitness 1E120
    for i in range(NUM_MAGNITUDES):
        probs[<int>sorted_fit[i][1]] = double_probs[i] # pow(i+1, -P_MAX)


cdef int compare(const_void * pa, const_void * pb):
    cdef double a = (<double *>pa)[0]
    cdef double b = (<double *>pb)[0]
    if a < b:
        return -1
    elif a > b:
        return 1
    else:
        return 0


cdef prune(Node* self, unsigned int distance):
    cdef double min_fit = 1E121
    cdef unsigned int i
    for i in range(NUM_SYM):
        if self.fitness[i] < min_fit:
            min_fit = self.fitness[i]
    for i in range(NUM_NT):
        if self.nodes_left[i] != NULL:
            if self.fitness[i] == min_fit and distance == 0:
                prune(self.nodes_left[i], distance)
            else:
                prune(self.nodes_left[i], distance+1)
    for i in range(NUM_DNT):
        if self.nodes_right[i] != NULL:
            if self.fitness[i] == min_fit and distance == 0:
                prune(self.nodes_right[i], distance)
            else:
                prune(self.nodes_right[i], distance+1)
    if distance > MIN_PRUNE_DISTANCE:
        for i in range(NUM_NT):
            if i!=SYM_CON:
                if self.nodes_left[i] != NULL:
                    destroy_tree(self.nodes_left[i], 0)
                    self.nodes_left[i] = NULL
            else:
                if self.cnode != NULL:
                    destroy_double_tree(self.cnode)
                    self.cnode = NULL
        for i in range(NUM_DNT):
            if self.nodes_right[i] != NULL:
                destroy_tree(self.nodes_right[i], 0)
                self.nodes_right[i] = NULL


cdef destroy_tree(Node* self, const unsigned int root):
    cdef unsigned int k
    for k in range(NUM_NT):
        if self.nodes_left[k] != NULL:
            destroy_tree(self.nodes_left[k], 0)
    for k in range(NUM_DNT):
        if self.nodes_right[k] != NULL:
            destroy_tree(self.nodes_right[k], 0)
    if self.cnode != NULL:
        destroy_double_tree(self.cnode)
    if root != 1:
        free(self)


cdef destroy_const_tree(ConstNode* self):
    for k in range(NUM_CONST_SYM):
        if self.nodes[k] != NULL:
            destroy_const_tree(self.nodes[k])
    free(self)


cdef destroy_double_tree(FloatNode* self):
    for k in range(NUM_MAGNITUDES):
        if self.nodes[k] != NULL:
            destroy_const_tree(self.nodes[k])

    free(self)


cdef best_to_string(Node* self):
    cdef unsigned int i_min = 0
    i_min = self.best
    if i_min < SNTERMINAL_IDX:
        return '(%s %c %s)' %   (best_to_string(self.nodes_left[i_min]), symbols[i_min],
                                best_to_string(self.nodes_right[i_min]))
    elif i_min < TERMINAL_IDX:
        return '%c(%s)' % (symbols[i_min], best_to_string(self.nodes_left[i_min]))
    elif i_min == SYM_CON:
        return '%f' % get_double_result(self.cnode)
    else:
        return '%c' % symbols[i_min]

def start(seed, target_function, target_degree=0):
    train_x_np, train_y_np, test_x_np, test_y_np = target_function(random.Random(), target_degree)
    cdef Node n
    cdef double train_x[NUM_TRAIN_CASES][NUM_VARS]
    cdef double train_result[NUM_TRAIN_CASES]
    cdef double train_y[NUM_TRAIN_CASES]
    cdef double test_x[NUM_TEST_CASES][NUM_VARS]
    cdef double test_result[NUM_TEST_CASES]
    cdef double test_y[NUM_TEST_CASES]
    cdef int i, j
    fitness_evolution = []

    for i in range(NUM_TRAIN_CASES):
        for j in range(NUM_VARS):
            train_x[i][j] = train_x_np[i,j]
        train_y[i] = train_y_np[i, 0]
    for i in range(NUM_TEST_CASES):
        for j in range(NUM_VARS):
            test_x[i][j] = test_x_np[i,j]
        test_y[i] = test_y_np[i, 0]

    cdef double min_fitness = 1E120
    cdef double min_squared_error = 1E120
    cdef double error, squared_error, fitness
    cdef int hit
    cdef int depth
    cdef double last_pre_fitness = 1E120

    create_node(&n, 0)
    for i in range(NUM_ITER):
        error = 0
        squared_error = 0
        fitness = 0
        hit = 0
        create_tree(&n)

        get_result(&n, NUM_TRAIN_CASES, train_x, train_result)
        for j in range(NUM_TRAIN_CASES):
            error = fabs(train_result[j]-train_y[j])
            squared_error += (train_result[j]-train_y[j])**2
            fitness += error
            if error < 0.01:
                hit += 1
        fitness = squared_error
        propagate_fitness(&n, fitness)
        if i % PRUNE_STEP == 0 and i > 0:
            prune(&n, 0)
        if fitness < min_fitness:
            min_fitness = fitness
            set_best(&n)
        if squared_error < min_squared_error:
            min_squared_error = squared_error
        if i%5000 == 0:
            fitness_evolution.append(min_squared_error/<double>NUM_TRAIN_CASES)
        if i % 50000 == 0 and i>0:
            print('%i Evaluation %i: MSE: %E' %(seed, i, min_fitness/<double>NUM_TRAIN_CASES))

    set_current_to_best(&n)
    get_result(&n, NUM_TEST_CASES, test_x, test_result)
    squared_error = 0
    for j in range(NUM_TEST_CASES):
        squared_error += (test_result[j]-test_y[j])**2

    print('%i Done. Evaluation %i: Train: %f Test: %f' %(seed, i, min_fitness/<double>NUM_TRAIN_CASES,
        squared_error/<double>NUM_TEST_CASES))
    print(best_to_string(&n))
    destroy_tree(&n, 1)
    return min_fitness/<double>NUM_TRAIN_CASES, squared_error/<double>NUM_TEST_CASES, fitness_evolution
