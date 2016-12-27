from path_class import *
from gibbs_MJPs import sampleUI
from gibbs_MJPs import get_likelihood
from gibbs_MJPs import FF
from gibbs_MJPs import BS
from math import log
from numpy import random
from math import pi
from math import sqrt
from model_parameters import *

def MHu_sampler_one(observation, ST_old, k, alpha, mu, lamb, var):
    # basic information
    initial_pi = copy.deepcopy(ST_old.initial_pi)
    N = len(initial_pi)
    t_start = ST_old.t_start
    t_end = ST_old.t_end
    # copy old parameters
    matrix_old = copy.deepcopy(ST_old.rate_matrix)
    alpha_old = alpha
    # Step 1 Propose a new alpha based on log normal distribution with variance var
    alpha_new = propose_alpha(alpha_old, var)
    matrix_new = constructor_rate_matrix(alpha_new)
    OMEGA_old = get_omega(alpha_old, k / 2.)
    OMEGA_new = get_omega(alpha_new, k / 2.)
    OMEGA = OMEGA_old + OMEGA_new
    uipath_old = sampleUI(ST_old, OMEGA)
    # calculate likelihood
    likelihood = get_likelihood(observation, uipath_old.T, N, [t_start, t_end])
    # Forward calculate P(Y | W, alpha_old, beta_old)
    logp_old,  ALPHA_old = FF(likelihood, initial_pi, OMEGA, uipath_old)
    # Temperarily change the rate matrix in order to cal marginal probability
    uipath_old.rate_matrix = matrix_new
    # Forward calculate P(Y | W, alpha_old, beta_new)
    logp_new,  ALPHA_new = FF(likelihood, initial_pi, OMEGA, uipath_old)
    # Step 2 Sample S,T
    accept_rate = logp_new - logp_old  + mu * (log(alpha_new) - log(alpha_old)) - lamb * (alpha_new - alpha_old)
    accept_rate = min(0, accept_rate)
    if log(random.uniform()) < accept_rate: # proposed alpha is accepted
        alpha_old = alpha_new
        ST_new = BS(initial_pi, ALPHA_new, likelihood, OMEGA, uipath_old)
    else: # rejected
        uipath_old.rate_matrix = matrix_old
        ST_new = BS(initial_pi, ALPHA_old, likelihood, OMEGA, uipath_old)
    ST_new.delete_virtual()
    return ST_new, alpha_old

def MHusampler(observation, pi_0, sample_n, T_interval, k, mu, lamb, var):
    alpha_list = []
    ST_list = []
    alpha_old = 3.0
    rate_matrix = constructor_rate_matrix(alpha_old)
    sample_old = MJPpath(t_start=T_interval[0], t_end=T_interval[1], rate_matrix=rate_matrix, initial_pi=pi_0)
    sample_old.generate_newpath()
    ST_old = sample_old
    ST_list.append(ST_old)
    for i in range(sample_n):
        ST_new, alpha_new = MHu_sampler_one(observation, ST_old, k, alpha_old, mu, lamb, var)
        ST_list.append(copy.deepcopy(ST_new))
        alpha_list.append(alpha_new)
        ST_old, alpha_old = ST_new, alpha_new
    return ST_list, alpha_list


def Full_GBS_u(timebreaks, pi_0, sample_n, T_interval, k, mu, lamb, var):
    observation = Observation()
    alpha_list = []
    uipath_list = []
    alpha_old = 1.0
    rate_matrix = constructor_rate_matrix(alpha_old)
    sample_old = MJPpath(t_start=T_interval[0], t_end=T_interval[1], rate_matrix=rate_matrix, initial_pi=pi_0)
    sample_old.generate_newpath()
    observation.sample_observation(timebreaks, sample_old)
    uipath_list.append(sample_old)
    for i in range(sample_n):
        my_print(i, 1000)
        observation.sample_observation(timebreaks, sample_old)
        sample_new, alpha_new = MHu_sampler_one(observation, sample_old, k, alpha_old, mu, lamb, var)
        uipath_list.append(copy.deepcopy(sample_new))
        alpha_list.append(alpha_new)
        sample_old, alpha_old = sample_new, alpha_new
    return uipath_list, alpha_list
