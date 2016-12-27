from path_class import *
from gibbs_MJPs import *
from gibbs_new_model import *
from MH_new_model3 import *
import json
import pickle
#from numpy import random
#from scipy.stats import gamma

pi_0 = [0.1] * 10
timebreaks = range(1, 20)
sample_n = 100000
T_interval = [0, 20]
# shape of alpha
mu = 3.
# scale of alpha
lamb = 2.
# shape of beta
omega = 5.
# scale of beta
theta = 2.
var = 1
SEED = 3
print "seed:", SEED
random.seed(SEED)
alpha = numpy.random.gamma(mu, 1./ lamb)
beta = numpy.random.gamma(omega, 1./ theta)
print alpha, beta
f = open('_data_EXPtruealpha_beta_d3', 'w')
json.dump([alpha, beta], f)
f.close()

f = open('_data_EXPprior_d3', 'w')
json.dump([mu, lamb, omega, theta], f)
f.close()

f = open('_data_EXPTime_interval_d3', 'w')
json.dump(T_interval, f)
f.close()

def f(alpha, beta, i, j):
    return alpha * exp(- beta / (i + j))

rate_matrix = constructor_rate_matrix(alpha, beta, 10)
base_path = MJPpath(t_start= T_interval[0], t_end= T_interval[1], rate_matrix= rate_matrix, initial_pi= pi_0)
base_path.generate_newpath()

file = open('_data_based_path_d3','w')
pickle.dump(base_path,file,0)
file.close()

observation = Observation(t_start= T_interval[0], t_end= T_interval[1])
observation.sample_observation(timebreaks, base_path)
observation.info()
file = open('_data_observation_d3','w')
pickle.dump(observation,file,0)
file.close()
