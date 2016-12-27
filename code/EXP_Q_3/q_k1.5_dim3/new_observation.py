from path_class import *
from gibbs_MJPs import *
from gbs_immi import *
from model_parameters import *
from mh_immi import *
import json
import pickle
from numpy import random
#from numpy import random
#from scipy.stats import gamma
DIM = 3
pi_0 = [1.0 / DIM] * DIM
timebreaks = range(1, 20)
T_interval = [0, 20]
# shape of alpha
mu = 3.
# scale of alpha
lamb = 2.
# shape of beta
omega = 5.
# scale of beta
theta = 2.

random.seed(DIM)
alpha = numpy.random.gamma(mu, 1./ lamb)
beta = numpy.random.gamma(omega, 1./ theta)
print alpha, beta
f = open('_data_EXPtruealpha_beta_immi_d3', 'w')
json.dump([alpha, beta], f)
f.close()

f = open('_data_EXPprior_immi_d3', 'w')
json.dump([mu, lamb, omega, theta], f)
f.close()

f = open('_data_EXPTime_interval_immi_d3', 'w')
json.dump(T_interval, f)
f.close()

rate_matrix = constructor_rate_matrix(alpha, beta, DIM)
base_path = MJPpath(t_start= T_interval[0], t_end= T_interval[1],rate_matrix= rate_matrix, initial_pi= pi_0)
base_path.generate_newpath()

file = open('_data_based_path_immi_d3','w')
pickle.dump(base_path,file,0)
file.close()
base_path.info()

observation = Observation(t_start= T_interval[0], t_end= T_interval[1])
observation.sample_observation(timebreaks, base_path)
observation.info()

file = open('_data_observation_immi_d3','w')
pickle.dump(observation,file,0)
file.close()

