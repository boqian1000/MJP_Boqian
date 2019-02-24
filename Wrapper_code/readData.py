from path_class import *
import pickle
file = open('Ecoli_lead_inner.txt', 'r')
tempdata = file.readlines()
data = [float(z) for z in tempdata[0].split(' ')]
data = data[2: len(data)]
obs1 = Observation(O=[1] * len(data), T = data, t_start= 0.0, t_end = 2319.838)
f = open('lead_inner', 'w')
pickle.dump(obs1, f, 0)
f.close()

