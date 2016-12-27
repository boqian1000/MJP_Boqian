from path_class import *
from gibbs_MJPs import *
from gibbs_JC import *
from MH_JC3 import *
from model_parameters import *
import json
import pickle
from numpy import random

random.seed(4)
pi_0 = [1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1. / 4.]
timebreaks = range(1, 20)
T_interval = [0, 20]
# shape of alpha
mu = 3.
# scale of alpha
lamb = 2.

alpha = numpy.random.gamma(mu, 1./ lamb)
print alpha
f = open('_data_EXPtruealpha', 'w')
json.dump(alpha, f)
f.close()

f = open('_data_EXPprior', 'w')
json.dump([mu, lamb], f)
f.close()

f = open('_data_EXPTime_interval', 'w')
json.dump(T_interval, f)
f.close()

rate_matrix = constructor_rate_matrix(alpha)
print rate_matrix
base_path = MJPpath(t_start= T_interval[0], t_end= T_interval[1],rate_matrix= rate_matrix, initial_pi= pi_0)
base_path.generate_newpath()
base_path.info()

file = open('_data_based_path','w')
pickle.dump(base_path,file,0)
file.close()

observation = Observation(t_start= T_interval[0], t_end= T_interval[1])
observation.sample_observation(timebreaks, base_path)
observation.info()
file = open('_data_observation','w')
pickle.dump(observation,file,0)
file.close()
