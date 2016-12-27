import time
import sys
import string
from path_class import *
from gibbs_MJPs import *
from dist_gbs import *
from mh_immi import *
from mh_immi_old import *
from model_parameters import *
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
    import json
    with open('_data_EXPTime_interval_immi_d3', 'r') as f:
        T_interval = json.load(f)
        f.close()
    with open('_data_EXPprior_immi_d3', 'r') as f:
        mu, lamb, omega, theta = json.load(f)
        f.close()
    with open('_data_observation_immi_d3', 'r') as f:
        observation = pickle.load(f)
        f.close()
    observation.info()
    pi_0 = [1. / 3] * 3
    sample_n = 5000
    k = 2
    observation1 = copy.deepcopy(observation)
    observation2 = copy.deepcopy(observation)
    sample_list, alpha_list, beta_list, loc_alpha_list, loc_beta_list, rate_alpha_list, rate_beta_list = BGsampler_immi(observation1, pi_0, 2 * sample_n, T_interval, k, mu, lamb, omega, theta)
    sample_mid = copy.deepcopy(sample_list[sample_n - 1])
    alpha_list2, beta_list2, loc_alpha_list2, loc_beta_list2, rate_alpha_list2, rate_beta_list2 = BGsampler_ST_immi(observation2, sample_mid, sample_n, T_interval, k, mu, lamb, omega, theta)
    my_save(alpha_list, "alpha1")
    my_save(beta_list, "beta1")
    my_save(alpha_list2, "alpha2")
    my_save(beta_list2, "beta2")
    my_save(loc_alpha_list, "loc_alpha1")
    my_save(loc_alpha_list2, "loc_alpha2")
    my_save(loc_beta_list, "loc_beta1")
    my_save(loc_beta_list2, "loc_beta2")
    my_save(rate_alpha_list, "rate_alpha1")
    my_save(rate_alpha_list2, "rate_alpha2")
    my_save(rate_beta_list, "rate_beta1")
    my_save(rate_beta_list2, "rate_beta2")


def main():
    arguments = validateInput(sys.argv)
    i = arguments
    simu_one(1)

if __name__ == '__main__':
    main()