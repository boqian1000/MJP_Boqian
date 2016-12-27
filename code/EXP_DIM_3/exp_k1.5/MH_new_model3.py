## This one is the latest version AUG 21th
from path_class import *
from gibbs_MJPs import sampleUI
from gibbs_MJPs import get_likelihood
from math import log
from numpy import random
from scipy.stats import gamma
from scipy.stats import poisson
from math import pi
from math import sqrt

def get_omega(mat, k):
    return(max(-mat.diagonal()) * k)

def my_f(alpha, beta, i, j):
    return alpha * exp(- beta * 1.0 / (i + j))

def my_F(alpha, beta, i, d):
    ret = 0
    for j in range(d):
        if i != j:
            ret += my_f(alpha, beta, i, j)
    return ret

def constructor_rate_matrix(alpha, beta, d):
    mat = numpy.identity(d, 'float')
    for i in range(d):
        mat[i][i] = 0
        for j in range(d):
            if i != j:
                mat[i][j] = my_f(alpha, beta, i, j)
                mat[i][i] = mat[i][i] - mat[i][j]
    return mat

def propose(alpha_old, beta_old, var=0.01):
    alpha_new = exp(random.normal(log(alpha_old), sqrt(var)))
    beta_new = exp(random.normal(log(beta_old), sqrt(var)))
    return alpha_new, beta_new

def propose_p(beta_old, beta_new, var): # ignoring the const term
    prob = (log(beta_old) - log(beta_new)) ** 2 / var
    prob = exp(-prob / 2.) / beta_new #/ (sqrt(2 * pi)) / sqrt(var)
    return prob

def prior(beta, omega, theta):
    prob = gamma.pdf(beta, omega, scale=1./theta)
    return prob

def FF(likelihood, initial_pi, OMEGA, path): #
    # X is the corresponding observations X1, ... ,XT
    N = len(initial_pi)
    rate_matrix = copy.deepcopy(path.rate_matrix)
    B = numpy.identity(N) + rate_matrix / float(OMEGA)
    alpha = [initial_pi]
    M = []
    M.append(1)
    # forward filtering:
    # (len(path.T) + 1) * N dimensional likelihood matrix
    for t in range(1, len(path.T) + 1):
        temp = [0.0] * N
        # j-th coordinate
        for j in range(N):
            for k in range(N):
                temp[j] += alpha[t - 1][k] * likelihood[t - 1][k] * B[k][j]
        maxt = max(temp)
        M.append(maxt)
        temp2 = [x / maxt for x in temp]
        alpha.append(temp2)
    # backward sampling:
    newS = []
    beta = array(likelihood[-1]) * array(alpha[-1])
    p_marginal = sum(beta)
    logm = 0.0
    for m in M:
        logm += log(m)
    log_p = log(p_marginal) + logm
    return log_p, alpha

def BS(initial_pi, alpha, likelihood, OMEGA, path): # Here path means path containing virtual jumps
    # X is the corresponding observations X1, ... ,XT
    N = len(initial_pi)
    rate_matrix = copy.deepcopy(path.rate_matrix)
    path_times = copy.deepcopy(path.T)
    B = numpy.identity(N) + rate_matrix / float(OMEGA)
    t_start = path.t_start
    t_end = path.t_end
    newS = []
    beta = array(likelihood[-1]) * array(alpha[-1])
    temp = sample_from_Multi(beta)
    newS.append(temp)
    iter = range(len(path_times))
    iter.reverse()
    for t in iter:
        beta = [0.0] * N
        for i in range(N):
            beta[i] = alpha[t][i] * likelihood[t][i] * B[i][newS[-1]]
        temp = sample_from_Multi(beta)
        newS.append(temp)
    newS.reverse()
    MJPpath_new= MJPpath(newS, path_times, t_start, t_end, rate_matrix, initial_pi)
    return MJPpath_new

# parameters = [alpha, beta]
def MHu_sampler_one(observation, ST_old, k, parameters, mu, lamb, omega, theta, var):
    # basic information
    initial_pi = copy.deepcopy(ST_old.initial_pi)
    N = len(initial_pi)
    t_start = ST_old.t_start
    t_end = ST_old.t_end
    w = len(ST_old.T)
    # copy old parameters
    matrix_old = copy.deepcopy(ST_old.rate_matrix)
    alpha_old = parameters[0]
    beta_old = parameters[1]
    # Step 1 Propose a theta* based on log normal distribution with variance var
    alpha_new, beta_new = propose(alpha_old, beta_old, var)
    matrix_new = constructor_rate_matrix(alpha_new, beta_new, N)
    OMEGA_old = get_omega(matrix_old, k)
    OMEGA_new = get_omega(matrix_new, k)
    OMEGA = max(OMEGA_new, OMEGA_old)
    #    OMEGA = OMEGA_new + OMEGA_old
    # Step 2 Sample W*:
    uipath_old = sampleUI(ST_old, OMEGA)
    # calculate likelihood
    likelihood = get_likelihood(observation, uipath_old.T, N, [t_start, t_end])
    # Forward calculate P(Y | W, alpha_old, beta_old)
    logp_old,  ALPHA_old = FF(likelihood, initial_pi, OMEGA, uipath_old)
    # Temperarily change the rate matrix in order to cal marginal probability
    uipath_old.rate_matrix = matrix_new
    # Forward calculate P(Y | W, alpha_old, beta_new)
    logp_new,  ALPHA_new = FF(likelihood, initial_pi, OMEGA, uipath_old)
    # Step 3 decide whether exchange theta* and theta
    accept_rate = logp_new - logp_old + mu * (log(alpha_new) - log(alpha_old)) - lamb * (alpha_new - alpha_old) + omega * (log(beta_new) - log(beta_old)) - theta * (beta_new - beta_old)
#    accept_rate *= propose_p(beta_new, beta_old, var) * propose_p(alpha_new, alpha_old, var) / propose_p(beta_old, beta_new, var) / propose_p(alpha_old, alpha_new, var)
    accept_rate = min(0, accept_rate)
    if log(random.uniform()) < accept_rate: # proposed beta is accepted
        beta_old = beta_new
        alpha_old = alpha_new
        ST_new = BS(initial_pi, ALPHA_new, likelihood, OMEGA, uipath_old)
    else: # rejected
        uipath_old.rate_matrix = matrix_old
        ST_new = BS(initial_pi, ALPHA_old, likelihood, OMEGA, uipath_old)
    # Step 4 Delete virtual jumps
    ST_new.delete_virtual()
    return ST_new, alpha_old, beta_old

def MHusampler(observation, pi_0, sample_n, T_interval, k, mu, lamb, omega, theta, var):
    alpha_list = []
    beta_list = []
    ST_list = []
    alpha_old = 2.0
    beta_old = 1.5
    rate_matrix = constructor_rate_matrix(alpha_old, beta_old, len(pi_0))
    ST_old = MJPpath(t_start=T_interval[0], t_end=T_interval[1], rate_matrix=rate_matrix, initial_pi=pi_0)
    ST_old.generate_newpath()
    ST_list.append(copy.deepcopy(ST_old))
    for i in range(sample_n):
        #        my_print(i, 1000)
        ST_new, alpha_new, beta_new= MHu_sampler_one(observation, ST_old, k, [alpha_old, beta_old], mu, lamb, omega, theta, var)
        ST_list.append(copy.deepcopy(ST_new))
        alpha_list.append(alpha_new)
        beta_list.append(beta_new)
        ST_old, alpha_old, beta_old = ST_new, alpha_new, beta_new
    return ST_list, alpha_list, beta_list
