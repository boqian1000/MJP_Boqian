from path_class import *
from gibbs_MJPs import sampleUI
from gibbs_MJPs import get_likelihood
from gibbs_MJPs import FFBS
from gibbs_MJPs import BGsampler_one
from math import log

#[t_start, t_1, t_2,...,t_N, t_end]
# #[s_0, s_1, s_2,...,s_N]
#input: alpha(scale par), beta , dimension (d)
#output: immigration model transition matrix

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

def pos_beta(path, omega, lamb, mu, theta, alpha, beta):
    S = path.S
    rate_new = lamb
    T = copy.deepcopy(path.T)
    T.insert(0, path.t_start)
    T.append(path.t_end)
    d = len(path.initial_pi)
    n = len(S)
    for i in range(n):
        rate_new += my_F(alpha, beta, S[i], d) * (T[i + 1] - T[i]) / alpha
    left = 1.
    for i in range(n - 1):
        left *= (my_f(alpha, beta, S[i], S[i + 1]) / alpha)
    return left * (beta ** (omega - 1)) * exp( - theta * beta) / (rate_new ** (mu + n - 1))
# NOTICE:!! here, n is not exactly the same as the "n" in the papar, n = "n" + 1 because of the s0 added.
# prob_new / prob_old (for efficient to do so)
def pos_beta_ratio(path, omega, lamb, mu, theta, alpha, beta, beta_new):
    S = path.S
    rate = lamb
    rate_new = lamb
    T = copy.deepcopy(path.T)
    T.insert(0, path.t_start)
    T.append(path.t_end)
    d = len(path.initial_pi)
    n = len(S)
    for i in range(n):
        rate_new += my_F(alpha, beta_new, S[i], d) * (T[i + 1] - T[i]) / alpha
        rate += my_F(alpha, beta, S[i], d) * (T[i + 1] - T[i]) / alpha
    left = 1.
    for i in range(n - 1):
        left *= my_f(alpha, beta_new, S[i], S[i + 1]) / my_f(alpha, beta, S[i], S[i + 1])
    return left * ((beta_new / beta) ** (omega - 1)) * exp( - theta * (beta_new - beta)) / ((rate_new / rate) ** (mu + n - 1))

def propose_beta(beta_old, var=0.01):
    beta_new = exp(random.normal(log(beta_old), sqrt(var)))
    return beta_new

def propose_p_beta(beta_old, beta_new, var): # ignoring the const term
    prob = (log(beta_old) - log(beta_new)) ** 2 / var
    prob = exp(-prob / 2.) / beta_new #/ (sqrt(2 * pi)) / sqrt(var)
    return prob

def BGsampler_one_A(observation, MJPpath0, k, mu, lamb, omega, theta, alpha, beta, var=0.01):
    #initialization
    d = len(MJPpath0.rate_matrix)
#    OMEGA = get_omega(MJPpath0.rate_matrix, k)
    OMEGA = k * alpha * (d - 1)
    new_path = BGsampler_one(observation, MJPpath0, OMEGA)
    # Sample alpha & beta| S, T, y:
    # Apply Metropolis Hasting to sample beta | S, T, y
    n = len(new_path.T)
    # propose new beta
    beta_new = propose_beta(beta, var)
    #print pos_beta(new_path, omega, lamb, mu, theta, f, alpha, beta_new), propose_p_beta(beta_new, beta, var), pos_beta(new_path, omega, lamb, mu, theta, f, alpha, beta), propose_p_beta(beta, beta_new, var)
    acc_rate = pos_beta_ratio(new_path, omega, lamb, mu, theta, alpha, beta, beta_new) * propose_p_beta(beta_new, beta, var) / propose_p_beta(beta, beta_new, var)
    acc_rate = min(1, acc_rate)
    if random.uniform() < acc_rate: # proposed beta is accepted
        beta = beta_new
    # Apply Gibbs to sample alpha | beta , S, T, y
    rate_new = lamb
    T = copy.deepcopy(new_path.T)
    T.insert(0, new_path.t_start)
    T.append(new_path.t_end)
    for i in range(len(T) - 1):
        rate_new += my_F(alpha, beta, new_path.S[i], d) * (T[i + 1] - T[i]) / alpha
    alpha = numpy.random.gamma((mu + n), 1./ rate_new)
    new_path.rate_matrix = constructor_rate_matrix(alpha, beta, d)
    return new_path, alpha, beta

def BGsampler(observation, pi_0, sample_n, T_interval, k, mu, lamb, omega, theta, var):
    # Initialization:
    alpha = 2.0
    beta = 3.0
    rate_matrix = constructor_rate_matrix(alpha, beta, len(pi_0))
    sample_old = MJPpath(t_start=T_interval[0], t_end=T_interval[1], rate_matrix=rate_matrix, initial_pi=pi_0)
    sample_old.generate_newpath()
    sample_list = []
    sample_list.append(copy.deepcopy(sample_old))
    alpha_list = []
    beta_list = []
    # Gibbs Iteration:
    for i in range(sample_n):
        my_print(i, 1000)
        sample_new, alpha, beta = BGsampler_one_A(observation, sample_old, k, mu, lamb, omega, theta, alpha, beta, var)
        sample_list.append(copy.deepcopy(sample_new))
        sample_old = sample_new
        alpha_list.append(alpha)
        beta_list.append(beta)
    return sample_list, alpha_list, beta_list


def Full_GBS(timebreaks, pi_0, sample_n, T_interval, k, mu, lamb, omega, theta, var):
    observation = Observation()
    alpha = 1.0
    beta = 1.0
    rate_matrix = constructor_rate_matrix(alpha, beta, len(pi_0))
    sample_old = MJPpath(t_start=T_interval[0], t_end=T_interval[1], S=[0, 1, 2, 1, 0, 1, 0], T=[2, 3, 4, 5, 6, 9], rate_matrix=rate_matrix, initial_pi=pi_0)
    sample_list = []
    sample_list.append(copy.deepcopy(sample_old))
    alpha_list = []
    beta_list = []
    for i in range(sample_n):
#        observation.sample_observation(timebreaks, sample_list[-1])
        observation.sample_observation(timebreaks, copy.deepcopy(sample_old))
        sample_new, alpha, beta = BGsampler_one_A(observation, sample_old, k, mu, lamb, omega, theta, alpha, beta, var)
        sample_list.append(copy.deepcopy(sample_new))
        sample_old = sample_new
        alpha_list.append(alpha)
        beta_list.append(beta)
    return sample_list, alpha_list, beta_list
''' ## main() __
pi_0 = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]
timebreaks = range(1, 50)
sample_n = 20000
T_interval = [0, 50]
# shape of alpha
mu = 3.
# scale of alpha
lamb = 2.
# shape of beta
omega = 3.
# scale of beta
theta = 2.
var = 0.25
print "started!"
samples, alpha_list, beta_list = Full_GBS(timebreaks, pi_0, sample_n, T_interval, 2, mu, lamb, omega, theta, var)


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

import matplotlib.pyplot as plt
import scipy.special as sps
count, bins, ignored = plt.hist(beta_list, 100, normed=True)
y = bins**(shape-1)*(numpy.exp(-bins* rate) /
                      (sps.gamma(shape)*scale**shape))
plt.plot(bins, y, linewidth=2, color='r')
plt.show()
import json

f = open('alphalist_GBS3252_apr15', 'w')
json.dump(alpha_list, f)
f.close()

g = open('betalist_GBS3252_apr15', 'w')
json.dump(beta_list, g)
g.close()


#__main()__:
pi_0 = [1./3., 1./3., 1./3.]
timebreaks = [2, 4, 6, 8]
sample_n = 3000
T_interval = [0, 10]

# shape of alpha
mu = 2.
# rate of alpha
lamb = 3.
# shape of beta
omega = 5.
# rate of beta
theta = 2.

samples, alpha_list, beta_list = Full_GBS_2(timebreaks, pi_0, sample_n, T_interval, 2, mu, lamb, omega, theta)
print numpy.mean(alpha_list), mu / lamb
print numpy.mean(beta_list), omega / theta
print numpy.var(alpha_list), (mu / lamb / lamb)
print numpy.var(beta_list), (omega / theta / theta)

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

import matplotlib.pyplot as plt
import scipy.special as sps
count, bins, ignored = plt.hist(beta_list, 100, normed=True)
y = bins**(shape-1)*(numpy.exp(-bins* rate) /
                      (sps.gamma(shape)*scale**shape))
plt.plot(bins, y, linewidth=2, color='r')
plt.show()

import json
f = open('alphalist_IMMI_', 'w')
json.dump(alpha_list, f)
f.close()

g = open('betalist_IMMI_', 'w')
json.dump(beta_list, g)
g.close()


h = open('N', 'w')
json.dump(N, h)
h.close()

import plotly.plotly as py
import plotly.graph_objs as go


trace0 = go.Histogram(
    x=alpha_list,
    opacity=0.75
)

trace1 = go.Histogram(
    x=beta_list,
    opacity=0.75
)

trace2 = go.Histogram(
    x=numpy.random.gamma(mu, 1./lamb, sample_n),
    opacity=0.75
)

trace3 = go.Histogram(
    x=numpy.random.gamma(omega, 1./theta, sample_n),
    opacity=0.75
)

data = [trace0, trace1, trace2, trace3]
layout = go.Layout(
    barmode='overlay'
)
fig = go.Figure(data=data, layout=layout)
plot_url = py.plot(fig, filename='alpha_beta_GIBS')


'''