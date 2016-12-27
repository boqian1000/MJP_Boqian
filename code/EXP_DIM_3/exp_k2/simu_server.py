import time
import sys
import string
from path_class import *
from gibbs_MJPs import *
from gibbs_new_model import *
from MH_new_model3 import *
from MH_new_model_old import *
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
    print "seed:", number
    random.seed(number)
    var = [0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 1.7]
    time_gbs_list = []
    time_mh_list = []
    time_omh_list = []
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
    pi_0 = [1./3.] * 3
    sample_n = 10000
    k = 2
    for cur_var in var:
        print number, "variance for current iteration is : ", cur_var
        if len(time_mh_list) > 0:
            print "mh_time :", time_mh_list
            print "gbs_time :", time_gbs_list
            print "omh_time", time_omh_list
        observation1 = copy.deepcopy(observation)
        observation2 = copy.deepcopy(observation)
        observation3 = copy.deepcopy(observation)
        t0 = time.clock()
        samples_GBS, gbsalpha_list, gbsbeta_list =BGsampler(observation1, pi_0, sample_n, T_interval, k, mu, lamb, omega, theta, cur_var)
        dt = time.clock() - t0
        time_gbs_list.append(dt)
        import json
        filename1 = "__" + str(number) + "_GBS_alpha" + str(cur_var)
        filename2 = "__" + str(number) + "_GBS_beta" + str(cur_var)
        my_save(gbsalpha_list, filename1)
        my_save(gbsbeta_list, filename2)
        t0 = time.clock()
        samples_MH, mhalpha_list, mhbeta_list = MHusampler(observation2, pi_0, sample_n, T_interval, k, mu, lamb, omega, theta, cur_var)
        dt = time.clock() - t0
        time_mh_list.append(dt)
        import json
        filename1 = "__" + str(number) + "_MH_alpha" + str(cur_var)
        filename2 = "__" + str(number) + "_MH_beta" + str(cur_var)
        my_save(mhalpha_list, filename1)
        my_save(mhbeta_list, filename2)
        t0 = time.clock()
        samples_oMH, omhalpha_list, omhbeta_list = MHusampler_old(observation3, pi_0, sample_n, T_interval, k, mu, lamb, omega, theta, cur_var)
        dt = time.clock() - t0
        time_omh_list.append(dt)
        import json
        filename1 = "__" + str(number) + "_oMH_alpha" + str(cur_var)
        filename2 = "__" + str(number) + "_oMH_beta" + str(cur_var)
        my_save(omhalpha_list, filename1)
        my_save(omhbeta_list, filename2)
    filename1 = "__" + str(number) + "_timelist_gbs"
    filename2 = "__" + str(number) + "_timelist_mh"
    filename3 = "__" + str(number) + "_timelist_omh"
    filename4 = "__" + str(number) + "var_list"
    my_save(time_gbs_list, filename1)
    my_save(time_mh_list, filename2)
    my_save(time_omh_list, filename3)
    my_save(var, filename4)


def main():
    arguments = validateInput(sys.argv)
    i = arguments
    simu_one(i)
if __name__ == '__main__':
    main()