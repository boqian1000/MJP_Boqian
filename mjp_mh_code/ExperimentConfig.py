#from path_class import *
import numpy as np
import json
def my_save(list, filename):
    f = open(filename, 'w')
    json.dump(list, f)
    f.close()

expInd = '7'
pi_0 = [.5, .5]
sample_n = 2000
T_interval = [0.0, 2319.838]
#T_interval = [0.0, 11]
k = 2
priors = np.array([[2, 2], [2, 3], [3, 2], [1, 2]])
#vars = [0.01, 0.1, 0.0001, 0.02]
#6
vars = [0.03, 0.3, 0.0003, 0.05]
burnin_size = 20000
#cov = np.matrix([[0.0008091786, 0.0025761028, -1.394539e-04, -0.0018867815],
#       [0.0025761028, 0.0451922495, -1.171265e-04, 0.0090474670],
#       [0.0001394539, -0.0001171265,  5.322672e-05 , 0.0004702824],
#       [-0.0018867815, 0.0090474670, 4.702824e-04, 0.0198454516]])