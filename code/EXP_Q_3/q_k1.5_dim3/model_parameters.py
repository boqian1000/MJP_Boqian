from math import log
from math import exp
from math import sqrt
from numpy import random
from numpy import array
import numpy

def my_print(i, j):
    if(i % j == 0):
        print i

def get_omega(matrix, k):
    OMEGA = k * (max(-matrix[-1][-1], -matrix[-2][-2]))
    return OMEGA

def constructor_rate_matrix(alpha, beta, d):
    mat = numpy.identity(d, 'float')
    for i in range(1, d - 1):
        mat[i][i + 1] = alpha
        mat[i][i - 1] = i * beta
        mat[i][i] = -alpha - i * beta
    mat[0][1] = alpha
    mat[0][0] = -alpha
    mat[-1][-2] = (d - 1) * beta
    mat[-1][-1] = -(d - 1) * beta
    return mat

def prior(beta, omega, theta):
    prob = gamma.pdf(beta, omega, scale=1./theta)
    return prob

def propose(alpha_old, beta_old, var=0.01):
    alpha_new = exp(random.normal(log(alpha_old), sqrt(var)))
    beta_new = exp(random.normal(log(beta_old), sqrt(var)))
    return alpha_new, beta_new

def propose_p(beta_old, beta_new, var): # ignoring the const term
    prob = (log(beta_old) - log(beta_new)) ** 2 / var
    prob = exp(-prob / 2.) / beta_new #/ (sqrt(2 * pi)) / sqrt(var)
    return prob
