import time
import sys
import string
from path_class import *
from gibbs_MJPs import *
from gibbs_new_model import *
from PMCMC_new_model import *
import json
import pickle
from numpy import random

# t0 = time.time()
#print time.time() - t0

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
    print "seed:", number
    random.seed(number)
    var = [0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 1.7]
    time_p_list = []
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
    sample_n = 3000
    N = 20
    for cur_var in var:
        print number, "variance for current iteration is : ", cur_var
        if len(time_p_list) > 0:
            print "p_time :", time_p_list
        observation2 = copy.deepcopy(observation)
        t0 = time.clock()
        samples_p, palpha_list, pbeta_list = PMMHsampler(observation2, pi_0, N, sample_n, T_interval, mu, lamb, omega, theta, cur_var)
        dt = time.clock() - t0
        time_p_list.append(dt)
        import json
        filename1 = "__" + str(number) + "_p_alpha" + str(cur_var)
        filename2 = "__" + str(number) + "_p_beta" + str(cur_var)
        my_save(palpha_list, filename1)
        my_save(pbeta_list, filename2)
    filename2 = "__" + str(number) + "_timelist_p"
    filename3 = "__" + str(number) + "p_var_list"
    my_save(time_p_list, filename2)
    my_save(var, filename3)


def main():
    arguments = validateInput(sys.argv)
    i = arguments
    simu_one(i)

if __name__ == '__main__':
    main()