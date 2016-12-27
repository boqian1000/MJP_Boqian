import time
import sys
import string
from path_class import *
from MH_new_model_iv import *
import json
import pickle
from numpy import random
import numpy as np
from numpy import array

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
    numpy.random.seed(number)
    with open('_data_EXPTime_interval_d3', 'r') as f:
        T_interval = json.load(f)
        f.close()
    with open('_data_EXPprior_d3', 'r') as f:
        mu, lamb, omega, theta = json.load(f)
        f.close()
    with open('_data_observation_d3', 'r') as f:
        observation = pickle.load(f)
        f.close()
#    observation.info()
    pi_0 = [1./3.] * 3
    sample_n = 49
    k = 2
    cur_var = 0.7
    ivlist = [[1,1], [5,5], [10,10], [15,15], [20, 20], [25, 25], [30,30], [35,35], [40,40], [45,45],[50,50],[60,60],[70,70],[80,80],[90,90], [100,100], [150, 150], [200, 200], [300, 300], [400, 400], [500, 500], [1000, 1000]]
    result = []
    for iv in ivlist:
        trans = []
        for i in range(20):
            observation2 = copy.deepcopy(observation)
            samples_MH, mhalpha_list, mhbeta_list = MHusampler(observation2, pi_0, sample_n, T_interval, k, mu, lamb, omega, theta, cur_var, iv)
            trans_time = []
            for path in samples_MH:
                trans_time.append(len(path.T))
            trans.append(trans_time)
        trans = array(trans)
        ret = np.mean(trans,axis=0)
        ret = list(ret)
        result.append(ret)
    my_save(result, "transit_time_RESULT")
def main():
    arguments = validateInput(sys.argv)
    i = arguments
    simu_one(1)

if __name__ == '__main__':
    main()
#set seed as 1