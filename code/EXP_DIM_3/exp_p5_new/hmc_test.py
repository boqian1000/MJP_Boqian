from path_class import *
from hmc_new_model import *
from gibbs_new_model import *
from PMCMC_new_model import *
from MH_new_model2 import *
from MH_new_model_old import *

import json
import pickle
import time

print "exp started!"
with open('_data_EXPTime_interval', 'r') as f:
    T_interval = json.load(f)
    f.close()


with open('_data_EXPprior', 'r') as f:
    mu, lamb, omega, theta = json.load(f)
    f.close()
with open('_data_observation', 'r') as f:
    observation = pickle.load(f)
    f.close()

with open('_data_based_path', 'r') as f:
    base = pickle.load(f)
    f.close()


pi_0 = [1./3., 1./3., 1./3.]
sample_n = 10000
var = 0.5
k = 2
N = 30
time_list = []
def my_save(list, filename):
    import json
    f = open(filename, 'w')
    json.dump(list, f)
    f.close()

def compare_hist(list1, list2):
    import plotly.plotly as py
    import plotly.graph_objs as go
    trace0 = go.Histogram(
        x=list1,
        opacity=0.75
    )
    trace1 = go.Histogram(
        x=list2,
        opacity=0.75
    )
    data = [trace0, trace1]
#    data = [trace0, trace1, trace2]
    layout = go.Layout(
        barmode='overlay'
    )
    fig = go.Figure(data=data, layout=layout)
    plot_url = py.plot(fig, filename='list1 vs list2')

def histo_one(path):
    S = copy.deepcopy(path.S)
    N = path.total_state
    his = [float(S.count(x))/float(len(S)) for x in range(N)]
    return his


def histo(path_list):
    histo_list = map(histo_one, path_list)
    return numpy.sum(histo_list, axis=0) / len(path_list)


def getnum_one(path):
    return len(path.S)


def getnum(path_list):
    num = map(getnum_one, path_list)
    return numpy.mean(num)


observation1 = copy.deepcopy(observation)
observation2 = copy.deepcopy(observation)
L = 10
stepsize = 0.2
m1 = 1.
m2 = 1.

t0 = time.clock()
samples_HMC, hmcalpha_list, hmcbeta_list = HMC(observation1, pi_0, sample_n, T_interval, k, L, stepsize, m1, m2, mu, lamb, omega, theta)
print "done"
dt = time.clock() - t0
time_list.append(dt)

t0 = time.clock()
samples_MH, mhalpha_list, mhbeta_list = MHusampler(observation2, pi_0, sample_n, T_interval, k, mu, lamb, omega, theta, var)
print "done"
dt = time.clock() - t0
time_list.append(dt)



my_save(mhalpha_list, "pos_mh_alpha3")
my_save(hmcalpha_list, "pos_hmc_alpha3")


my_save(mhbeta_list, "pos_mh_beta3")
my_save(hmcbeta_list, "pos_hmc_beta3")


my_save(time_list, "pos_time3")
