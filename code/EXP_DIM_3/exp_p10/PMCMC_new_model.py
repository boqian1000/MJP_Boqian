from path_class import *
from PMCMC import *
from math import log
from scipy.stats import gamma


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

def prior(beta, omega, theta):
    prob = gamma.pdf(beta, omega, scale=1./theta)
    return prob


#par_new = [alpha, beta]
def propose(alpha_old, beta_old, var=0.01):
    alpha_new = exp(random.normal(log(alpha_old), sqrt(var)))
    beta_new = exp(random.normal(log(beta_old), sqrt(var)))
    return alpha_new, beta_new

def propose_p(beta_old, beta_new, var): # ignoring the const term
    prob = (log(beta_old) - log(beta_new)) ** 2 / var
    prob = exp(-prob / 2.) / beta_new #/ (sqrt(2 * pi)) / sqrt(var)
    return prob


def PMMHsampler(observation, pi_0, N, sample_n, T_interval, mu, lamb, omega, theta, var):
    y = observation.O
    T = observation.T
    MJP_samples = []
    # Step 1
    alpha = 2.
    beta = 2.
    d = len(pi_0)
    rate_matrix = constructor_rate_matrix(alpha, beta, d)
    alpha_list = [alpha]
    beta_list = [beta]
    particles, weights, old_p = SMC_MJPs(y, T, rate_matrix, pi_0, T_interval, N)
    old_sample = sample_from_particles(particles, weights)
    MJP_samples.append(copy.deepcopy(old_sample))
    for iter in range(sample_n):
        alpha_new, beta_new = propose(alpha, beta, var)
        new_rate_matrix = constructor_rate_matrix(alpha_new, beta_new, d)
        particles, weights, new_p = SMC_MJPs(y, T, new_rate_matrix, pi_0, T_interval, N)
        new_sample = sample_from_particles(particles, weights)
        accept_rate = new_p * alpha_new * beta_new / (old_p * alpha * beta)
        accept_rate *= prior(beta_new, omega, theta) * prior(alpha_new, mu, lamb) / prior(beta, omega, theta) / prior(alpha, mu, lamb)
        accept_rate = min(1, accept_rate)
        if random.uniform() < accept_rate:
            old_sample = new_sample
            old_p = new_p
            alpha, beta = alpha_new, beta_new
        MJP_samples.append(copy.deepcopy(old_sample))
        alpha_list.append(alpha)
        beta_list.append(beta)
    return MJP_samples, alpha_list, beta_list

