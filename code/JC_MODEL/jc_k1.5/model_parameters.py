from math import log
from math import exp
from math import sqrt
from numpy import random
from numpy import array
import numpy

def my_print(i, j):
    if(i % j == 0):
        print i

def get_omega(alpha, k):
    ret = alpha * k * 3.
    return ret

def constructor_rate_matrix(alpha):
    mat = numpy.identity(4, 'float')
    for i in range(4):
        mat[i][i] = -3. * alpha
        for j in range(4):
            if i != j:
                mat[i][j] = alpha
    return mat

def prior(alpha, mu, lamb):
    prob = gamma.pdf(alpha, mu, scale=1./lamb)
    return prob

def propose_alpha(alpha_old, var=0.01):
    alpha_new = exp(random.normal(log(alpha_old), sqrt(var)))
    return alpha_new

def propose_p_alpha(alpha_old, alpha_new, var): # ignoring the const term
    prob = (log(alpha_old) - log(alpha_new)) ** 2 / var
    prob = exp(-prob / 2.) / alpha_new #/ (sqrt(2 * pi)) / sqrt(var)
    return prob

