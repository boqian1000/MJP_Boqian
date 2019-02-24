import time
import sys
import string
import numpy as np
import pickle
import time
from MH import *
from GBS import *
from ExperimentConfig import *
import json
import matplotlib.pyplot as plt

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
    if '-v' in args_map:
        varindex = float(args_map['-v'])
    return number, varindex

def experiment(seedIndex, kappa, pi_0, sample_n, T_interval, k, priors, burnin_size):
    """
    :param seedIndex:
    :param kappa:
    :param pi_0:
    :param sample_n:
    :param T_interval:
    :param k:
    :param priors:
    :param kappa:
    :param burnin_size:
    :return:
    """
    print "seed:", seedIndex
    np.random.seed(seedIndex)
#    multi = [.1, .3, .5, 1., 3., 5., 10.]

    with open('lag_inner', 'r') as f:
        observation = pickle.load(f)
        f.close()

    burnin_alpha_list, burnin_beta_list, burnin_lamb1_list, burnin_lamb2_list = GBSsampler(observation, pi_0, burnin_size, T_interval, k, priors)
    SimulationIndex = seedIndex

    filename1 = "__" + str(SimulationIndex) + "_burnin_MH_alpha_" + str(kappa)
    filename2 = "__" + str(SimulationIndex) + "_burnin_MH_beta_" + str(kappa)
    filename4 = "__" + str(SimulationIndex) + "_burnin_MH_lamb1_" + str(kappa)
    filename5 = "__" + str(SimulationIndex) + "_burnin_MH_lamb2_" + str(kappa)
    my_save(burnin_alpha_list, filename1)
    my_save(burnin_beta_list, filename2)
    my_save(burnin_lamb1_list, filename4)
    my_save(burnin_lamb2_list, filename5)
    X = np.stack((burnin_alpha_list, burnin_beta_list, burnin_lamb1_list, burnin_lamb2_list), axis=0)

    covariance = np.cov(X)
    cov = covariance * kappa


    alpha_list, beta_list, lamb1_list, lamb2_list = MHsampler(observation, pi_0, sample_n, T_interval, k, priors, cov)
    filename1 = "__" + str(SimulationIndex) + "_MH_alpha_" + str(kappa)
    filename2 = "__" + str(SimulationIndex) + "_MH_beta_" + str(kappa)
    filename4 = "__" + str(SimulationIndex) + "_MH_lamb1_" + str(kappa)
    filename5 = "__" + str(SimulationIndex) + "_MH_lamb2_" + str(kappa)
    my_save(alpha_list, filename1)
    my_save(beta_list, filename2)
    my_save(lamb1_list, filename4)
    my_save(lamb2_list, filename5)
    filename3 = "__" + str(SimulationIndex) + "_timelist_mh_" + str(kappa)

    alphagbs_list, betagbs_list, lamb1gbs_list, lamb2gbs_list = GBSsampler(observation, pi_0, sample_n + burnin_size, T_interval, k, priors)

    filename1 = "__" + str(SimulationIndex) + "_GBS_alpha_"
    filename2 = "__" + str(SimulationIndex) + "_GBS_beta_"
    filename4 = "__" + str(SimulationIndex) + "_GBS_lamb1_"
    filename5 = "__" + str(SimulationIndex) + "_GBS_lamb2_"
    alphagbs_list = alphagbs_list[burnin_size: len(alphagbs_list)]
    betagbs_list = betagbs_list[burnin_size: len(betagbs_list)]
    lamb1gbs_list = lamb1gbs_list[burnin_size: len(lamb1gbs_list)]
    lamb2gbs_list = lamb2gbs_list[burnin_size: len(lamb2gbs_list)]
    my_save(alphagbs_list, filename1)
    my_save(betagbs_list, filename2)
    my_save(lamb1gbs_list, filename4)
    my_save(lamb2gbs_list, filename5)

    num_bins = 20
    plt.hist(alpha_list, num_bins, histtype='bar', color='blue', label='MH')
    plt.hist(alphagbs_list, num_bins, histtype='bar', color='red', label='GBS')
    plt.legend(prop={'size': 20})
    plt.title("posterior alpha")
    plt.show()
    plt.savefig('alpha.pdf')



def main():
    arguments = validateInput(sys.argv)
    seedNumber, kappa = arguments
#    seedNumber, kappa = 1, 1
    experiment(seedNumber, kappa, pi_0, sample_n, T_interval, k, priors, burnin_size)
    ## -v from covariance multiplicative number
    ## -n from seed

if __name__ == '__main__':
    main()
