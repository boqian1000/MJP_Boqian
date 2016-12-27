from path_class import *
from gibbs_MJPs import sampleUI
from gibbs_MJPs import get_likelihood
from gibbs_MJPs import FFBS
from math import log
from numpy import random
from scipy.stats import gamma
from scipy.stats import poisson
from math import pi
from math import sqrt

def my_f(alpha0, beta0, i, j):
    alpha = abs(alpha0)
    beta = abs(beta0)
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

def dA(alpha, beta, d):
    dalpha = numpy.identity(d, 'float')
    for i in range(d):
        dalpha[i][i] = 0
        for j in range(d):
            if i != j:
                dalpha[i][j] = my_f(1., beta, i, j)
                dalpha[i][i] = dalpha[i][i] - dalpha[i][j]
    dalpha *= alpha / abs(alpha)
    dbeta = numpy.identity(d, 'float')
    for i in range(d):
        dbeta[i][i] = 0
        for j in range(d):
            if i != j:
                dbeta[i][j] = my_f(alpha, beta, i, j) * (beta / abs(beta)) / (-i - j)
                dbeta[i][i] = dbeta[i][i] - dbeta[i][j]
    return dalpha, dbeta

def get_OMEGA(alpha, beta, k, mat):
    N = len(mat)
    mxd = max(-mat.diagonal())
    dig = list(-mat.diagonal())
    I = dig.index(mxd)
    OMEGA = mxd * k
    dOMEGAdalpha = OMEGA / alpha
    dOMEGAdbeta = 0.0
    for j in range(1, N):
        dOMEGAdbeta += mat[I][j] * (-1.0 / (j + I)) * (beta / abs(beta))
    return OMEGA, dOMEGAdalpha, dOMEGAdbeta

def DB(alpha, beta, k, mat):
    d = len(mat)
    dadalpha, dadbeta = dA(alpha, beta, d)
    OMEGA, dOMEGAdalpha, dOMEGAdbeta = get_OMEGA(alpha, beta, k, mat)
    dBdalpha = dadalpha / OMEGA - mat * dOMEGAdalpha / OMEGA / OMEGA
    dBdbeta = dadbeta / OMEGA - mat * dOMEGAdbeta / OMEGA / OMEGA
    return dBdalpha, dBdbeta

def prior(par, mu, lamb, omega, theta):
    return gamma.pdf(abs(par[0]), mu, scale=1./lamb) * gamma.pdf(abs(par[1]), omega, scale=1./theta) / 4.

def HFF(initial_pi, observation, k, path, alpha, beta): #
    # X is the corresponding observations X1, ... ,XT
    rate_matrix = copy.deepcopy(path.rate_matrix)
    OMEGA = get_OMEGA(alpha, beta, k, rate_matrix)[0]
    path_times = copy.deepcopy(path.T)
    B = numpy.identity(rate_matrix.shape[0]) + rate_matrix / float(OMEGA)
    N = len(initial_pi)
    t_start = observation.t_start
    t_end = observation.t_end
    ALPHA = [initial_pi]
    M = [1.]
    # forward filtering:
    likelihood = get_likelihood(observation, path_times, N, [t_start, t_end])
    # (len(path.T) + 1) * N dimensional likelihood matrix
    DBDalpha, DBDbeta = DB(alpha, beta, k, rate_matrix)
    dalpha = []
    dbeta = []
    dalpha.append([0] * N)
    dbeta.append([0] * N)
    for t in range(1, len(path_times) + 1):
        temp = [0.0] * N
        tempdalpha = [0.0] * N
        tempdbeta = [0.0] * N
        # j-th coordinate
        for s in range(N):
            for v in range(N):
                temp[s] += ALPHA[-1][v] * likelihood[t - 1][v] * B[v][s]
                tempdalpha[s] += (dalpha[-1][v] * B[v][s] + ALPHA[-1][v] * DBDalpha[v][s])* likelihood[t - 1][v]
                tempdbeta[s] += (dbeta[-1][v] * B[v][s] + ALPHA[-1][v] * DBDbeta[v][s])* likelihood[t - 1][v]
        maxt = max(max(temp), max(tempdalpha), max(tempdbeta))
        M.append(maxt)
        temp2 = [x / maxt for x in temp]
        tempdalpha2 = [x / maxt for x in tempdalpha]
        tempdbeta2 = [x / maxt for x in tempdbeta]
        ALPHA.append(temp2)
        dalpha.append(tempdalpha2)
        dbeta.append(tempdbeta2)
    newS = []
    BETA = array(likelihood[-1]) * array(ALPHA[-1])
    Dalpha = array(likelihood[-1]) * array(dalpha[-1])
    Dbeta = array(likelihood[-1]) * array(dbeta[-1])
    logm = 0.0
    for m in M:
        logm += log(m)
    logp = log(sum(BETA)) + logm
    dalpha = sum(Dalpha) * exp(logm)
    dbeta = sum(Dbeta) * exp(logm)
    # print p_marginal, sum(Dalpha), sum(Dbeta)
    return logp, dalpha, dbeta, ALPHA, likelihood

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

def leapfrog(alpha_old, beta_old, k, observation, uipath, L=1, stepsize=0.1, m1=1, m2=1, mu=2, lamb=2, omega0=2, theta=2):
    path = copy.deepcopy(uipath)
    T = path.t_end - path.t_start
    N_W = len(path.T)
    initial_pi = copy.deepcopy(path.initial_pi)
    d = len(path.initial_pi)
    alpha = alpha_old
    beta = beta_old
    mat = constructor_rate_matrix(alpha, beta, d)
    OMEGA_old, dOMEGAdalpha, dOMEGAdbeta = get_OMEGA(alpha, beta, k, mat)
    p1 = numpy.random.normal(0, sqrt(m1))
    p2 = numpy.random.normal(0, sqrt(m2))
    p1_old = p1
    p2_old = p2
    logp, dp_dalpha, dp_dbeta, ALPHA_OLD, likelihood = HFF(initial_pi, observation, k, path, alpha, beta)
    OMEGA = OMEGA_old
    logp_old = logp
    dU_dalpha = - dp_dalpha * exp(-logp) + T * dOMEGAdalpha - N_W * dOMEGAdalpha / OMEGA - (mu - 1) / alpha + lamb * (alpha / abs(alpha))
    p1 -= dU_dalpha * stepsize / 2.
    dU_dbeta = - dp_dbeta * exp(-logp) + T * dOMEGAdbeta - N_W * dOMEGAdbeta / OMEGA - (omega0 - 1) / beta + theta * (beta / abs(beta))
    p2 -= dU_dbeta * stepsize / 2.
    for i in range(L):
        alpha += stepsize * p1 / m1
        beta += stepsize * p2 / m2
        path.rate_matrix = constructor_rate_matrix(alpha, beta, d)
        logp, dp_dalpha, dp_dbeta, ALPHA_NEW, l = HFF(initial_pi, observation, k, path, alpha, beta)
        OMEGA, dOMEGAdalpha, dOMEGAdbeta = get_OMEGA(alpha, beta, k, path.rate_matrix)
        dU_dalpha = - dp_dalpha * exp(-logp) + T * dOMEGAdalpha - N_W * dOMEGAdalpha / OMEGA - (mu - 1) / alpha + lamb * (alpha / abs(alpha))
        dU_dbeta = - dp_dbeta * exp(-logp) + T * dOMEGAdbeta - N_W * dOMEGAdbeta / OMEGA - (omega0 - 1) / beta + theta * (beta / abs(beta))
        if i < L - 1:
            p1 -= dU_dalpha * stepsize
            p2 -= dU_dbeta * stepsize
        else:
            p1 -= dU_dalpha * stepsize / 2.
            p2 -= dU_dbeta * stepsize / 2.
    p1 = - p1
    p2 = - p2
#    Uproposed = -logp + OMEGA * T - N_W * log(OMEGA) - (mu - 1) * log(abs(alpha)) - (omega0 - 1) * log(abs(beta)) + lamb * abs(alpha) + theta * abs(beta)
    Kproposed = p1 * p1 / 2. / m1 + p2 * p2 / 2. / m2
#    Uold = -logp_old + OMEGA_old * T - N_W * log(OMEGA_old) - (mu - 1) * log(abs(alpha_old)) - (omega0 - 1) * log(abs(beta_old)) + lamb * abs(alpha_old) + theta * abs(beta_old)
    Kold = p1_old * p1_old / 2. / m1 + p2_old * p2_old / 2. / m2
    acc_rate = logp - logp_old + (OMEGA_old - OMEGA) * T + (mu - 1) * (log(abs(alpha)) - log(abs(alpha_old))) + (omega0 - 1) * (log(abs(beta)) - log(abs(beta_old))) + lamb * (abs(alpha_old) - abs(alpha)) + theta * (abs(beta_old) - abs(beta))
    acc_rate += Kold - Kproposed
    if OMEGA / OMEGA_old == 0:
        alpha, beta, ALPHA_NEW = alpha_old, beta_old, ALPHA_OLD
        return alpha, beta, ALPHA_NEW, likelihood
    acc_rate += N_W * log(OMEGA / OMEGA_old)
    acc_rate = min(0, acc_rate)
#    acc_rate = exp(-Uproposed - Kproposed + Uold + Kold)
#    acc_rate = min(1, acc_rate)
    if log(random.uniform()) > acc_rate:
        alpha, beta, ALPHA_NEW = alpha_old, beta_old, ALPHA_OLD
    return alpha, beta, ALPHA_NEW, likelihood

# in MHsampler, all the trajectories contain virtual jumps
# parameters = [alpha, beta]
def HMC_sampler_one(alpha_old, beta_old, observation, uipath_old, k, parameters, L=10, stepsize=0.01, m1=1, m2=1, mu=2, lamb=2, omega=2, theta=2):
    matrix_old = copy.deepcopy(uipath_old.rate_matrix)
    initial_pi = copy.deepcopy(uipath_old.initial_pi)
    alpha, beta, ALPHA, likelihood = leapfrog(alpha_old, beta_old, k, observation, uipath_old, L, stepsize, m1, m2, mu, lamb, omega, theta)
    uipath_old.rate_matrix = constructor_rate_matrix(alpha, beta, len(matrix_old))
    OMEGA = get_OMEGA(alpha, beta, k, uipath_old.rate_matrix)[0]
#    print uipath_old.rate_matrix, OMEGA
    #ST_new = BS(initial_pi, ALPHA, likelihood, OMEGA, uipath_old)
    ST_new, p_old = FFBS(initial_pi,observation, OMEGA, uipath_old)
    ST_new.delete_virtual()
    uipath_new = sampleUI(ST_new, OMEGA)
    return uipath_new, alpha, beta

def HMC(observation, pi_0, sample_n, T_interval, k, L, stepsize, m1, m2, mu, lamb, omega, theta):
    uipath_list = []
    alpha_old = 2.
    beta_old = 2.
    alpha_list = [alpha_old]
    beta_list = [beta_old]
    rate_matrix = constructor_rate_matrix(alpha_old, beta_old, len(pi_0))
    OMEGA = get_OMEGA(alpha_old, beta_old, k, rate_matrix)[0]
    sample_old = MJPpath(t_start=T_interval[0], t_end=T_interval[1], rate_matrix=rate_matrix, initial_pi=pi_0)
    sample_old.generate_newpath()
    uipath_old = sampleUI(sample_old, OMEGA)
    uipath_list.append(copy.deepcopy(uipath_old))
    for i in range(sample_n):
        #if i % 1000 == 0 and i > 0:
        #            print i
        uipath_new, alpha_new, beta_new = HMC_sampler_one(alpha_old, beta_old, observation, uipath_old, k, [alpha_old, beta_old],L, stepsize, m1, m2, mu, lamb, omega, theta)
        uipath_list.append(copy.deepcopy(uipath_new))
        alpha_list.append(alpha_new)
        beta_list.append(beta_new)
        uipath_old, alpha_old, beta_old = uipath_new, alpha_new, beta_new
    return uipath_list, alpha_list, beta_list
'''

def Full_GBS_HMC(timebreaks, pi_0, sample_n, T_interval, k, L, stepsize, m1, m2, mu, lamb, omega, theta):
    observation = Observation()
    uipath_list = []
    alpha_old = 2.
    beta_old = 2.
    alpha_list = [alpha_old]
    beta_list = [beta_old]
    rate_matrix = constructor_rate_matrix(alpha_old, beta_old, len(pi_0))
    OMEGA = get_OMEGA(alpha_old, beta_old, k, rate_matrix)[0]
    sample_old = MJPpath(t_start=T_interval[0], t_end=T_interval[1], rate_matrix=rate_matrix, initial_pi=pi_0)
    sample_old.generate_newpath()
    observation.sample_observation(timebreaks, sample_old)
    uipath_old = sampleUI(sample_old, OMEGA)
    uipath_list.append(copy.deepcopy(uipath_old))
    for i in range(sample_n):
        if i % 1000 == 0 and i > 0:
            print i
        sample_old = copy.deepcopy(uipath_old)
        sample_old.delete_virtual()
        observation.sample_observation(timebreaks, sample_old)
        uipath_new, alpha_new, beta_new = HMC_sampler_one(alpha_old, beta_old, observation, uipath_old, k, [alpha_old, beta_old],L, stepsize, m1, m2, mu, lamb, omega, theta)
        uipath_list.append(copy.deepcopy(uipath_new))
        alpha_list.append(alpha_new)
        beta_list.append(beta_new)
        uipath_old, alpha_old, beta_old = uipath_new, alpha_new, beta_new
    return uipath_list, alpha_list, beta_list

 ## main() __
pi_0 = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]
timebreaks = [2, 4, 6, 8]
sample_n = 20000
T_interval = [0, 10]
# shape of alpha
mu = 3.
# scale of alpha
lamb = 2.
# shape of beta
omega = 5.
# scale of beta
theta = 2.
L = 1
stepsize = 0.2
m1 = 1.
m2 = 1.
samples, alpha_list, beta_list = Full_GBS_HMC(timebreaks, pi_0, sample_n, T_interval, 2, L, stepsize, m1, m2, mu, lamb, omega, theta)

print numpy.mean(alpha_list)
print numpy.mean(beta_list)

shape = mu
rate = lamb
scale = 1. / rate

import matplotlib.pyplot as plt
import scipy.special as sps
count, bins, ignored = plt.hist(alpha_list, 100, normed=True)
y = bins**(shape-1)*(numpy.exp(-bins* rate) /
                      (sps.gamma(shape)*scale**shape))
plt.plot(bins, y, linewidth=2, color='r')
plt.show()

shape = omega
rate = theta
scale = 1. / rate

count, bins, ignored = plt.hist(beta_list, 100, normed=True)
y = bins**(shape-1)*(numpy.exp(-bins* rate) /
                      (sps.gamma(shape)*scale**shape))
plt.plot(bins, y, linewidth=2, color='r')
plt.show()

import json

f = open('hmca', 'w')
json.dump(alpha_list, f)
f.close()

g = open('hmcb', 'w')
json.dump(beta_list, g)
g.close()
'''