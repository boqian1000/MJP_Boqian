import time
import sys
import string
from path_class import *
from gibbs_MJPs import *
from hmc_new_model import *
import json
import pickle
from numpy import random

def my_save(list, filename):
    import json
    f = open(filename, 'w')
    json.dump(list, f)
    f.close()

def parseArgs(args):
    #"""Parses arguments vector, looking for switches of the form -key {optional value}.
    #For example:
	#parseArgs([ 'template.py', '-v', '10', '-n' , '100'}"""
    args_map = {}
    curkey = None
    for i in xrange(1, len(args)):
        if args[i][0] == '-':
            args_map[args[i]] = True
            curkey = args[i]
        else:
            assert curkey
            args_map[curkey] = args[i]
            curkey = None
    return args_map

def validateInput(args):
    args_map = parseArgs(args)
    if '-n' in args_map:
        number = int(args_map['-n'])
    return number

def simu_one(number):
    #set a seed
    random.seed(number)
    L = [1, 3, 5, 10]
    time_hmc_list = []
    import json
    with open('_data_EXPTime_interval_d3', 'r') as f:
        T_interval = json.load(f)
        f.close()
    with open('_data_EXPprior_d3', 'r') as f:
        mu, lamb, omega, theta = json.load(f)
        f.close()
    with open('_data_observation_d3', 'r') as f:
        observation = pickle.load(f)
        f.close()
    pi_0 = [1./3., 1./3., 1./3.]
    sample_n = 20000
    k = 2
    stepsize = 0.05
    m1 = 1.
    m2 = 1.
    for l in L:
        print number, "variance for current L is : ", l
        if len(time_hmc_list) > 0:
            print "hmc_time :", time_hmc_list
        observation1 = copy.deepcopy(observation)
        t0 = time.clock()
        samples_HMC, hmcalpha_list, hmcbeta_list = HMC(observation1, pi_0, sample_n, T_interval, k, l, stepsize, m1, m2, mu, lamb, omega, theta)
        print "done"
        dt = time.clock() - t0
        time_hmc_list.append(dt)
        import json
        filename1 = "__" + str(number) + "_hmc_alpha" + str(l)
        filename2 = "__" + str(number) + "_hmc_beta" + str(l)
        my_save(hmcalpha_list, filename1)
        my_save(hmcbeta_list, filename2)
        filename = "__" + str(number) + "_timelist_hmc"
        filename4 = "__" + str(number) + "L_list"
        my_save(time_hmc_list, filename)
        my_save(L, filename4)

def main():
    arguments = validateInput(sys.argv)
    i = arguments
    simu_one(i)

if __name__ == '__main__':
    main()


