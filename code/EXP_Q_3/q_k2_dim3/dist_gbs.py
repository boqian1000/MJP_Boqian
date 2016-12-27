from path_class import *
from gibbs_MJPs import sampleUI
from gibbs_MJPs import get_likelihood
from gibbs_MJPs import FFBS
from gibbs_MJPs import BGsampler_one
from model_parameters import constructor_rate_matrix
from model_parameters import my_print

def BGsampler_one_A_immi(observation, MJPpath0, k, mu, lamb, omega, theta):
    #initialization
    OMEGA = k * max(-numpy.diag(MJPpath0.rate_matrix))
    new_path = BGsampler_one(observation, MJPpath0, OMEGA)
    # Sample alpha & beta| S, T, y:
    move = numpy.diff(new_path.S)
    Tau = new_path.stay_time_list()
    N = new_path.total_state
    rate_alpha = lamb + MJPpath0.t_end - MJPpath0.t_start - Tau[-1]
    rate_beta = theta + numpy.inner(Tau ,range(N))
    loc_alpha = mu + sum(move > 0)
    loc_beta = omega + sum(move < 0)
    alpha = numpy.random.gamma(loc_alpha, 1./ rate_alpha)
    beta = numpy.random.gamma(loc_beta, 1./ rate_beta)
    A_new = constructor_rate_matrix(alpha, beta, N)
    new_path.rate_matrix = copy.deepcopy(A_new)
    return new_path, alpha, beta, loc_alpha, loc_beta, rate_alpha, rate_beta

def BGsampler_one_ST_immi(observation, MJPpath0, mu, lamb, omega, theta):
    #initialization
    # Sample alpha & beta| S, T, y:
    new_path = MJPpath0
    move = numpy.diff(new_path.S)
    Tau = new_path.stay_time_list()
    N = new_path.total_state
    rate_alpha = lamb + MJPpath0.t_end - MJPpath0.t_start - Tau[-1]
    rate_beta = theta + numpy.inner(Tau ,range(N))
    loc_alpha = mu + sum(move > 0)
    loc_beta = omega + sum(move < 0)
    alpha = numpy.random.gamma(loc_alpha, 1./ rate_alpha)
    beta = numpy.random.gamma(loc_beta, 1./ rate_beta)
    return alpha, beta, loc_alpha, loc_beta, rate_alpha, rate_beta

def BGsampler_ST_immi(observation, MJPpath0, sample_n, T_interval, k, mu, lamb, omega, theta):
    # Initialization:
    alpha_list = []
    beta_list = []
    loc_alpha_list = []
    loc_beta_list = []
    rate_alpha_list = []
    rate_beta_list = []
    # Gibbs Iteration:
    for i in range(sample_n):
        alpha, beta, loc_alpha, loc_beta, rate_alpha, rate_beta = BGsampler_one_ST_immi(observation, MJPpath0, mu, lamb, omega, theta)
        alpha_list.append(alpha)
        beta_list.append(beta)
        loc_alpha_list.append(loc_alpha)
        loc_beta_list.append(loc_beta)
        rate_alpha_list.append(rate_alpha)
        rate_beta_list.append(rate_beta)
    return alpha_list, beta_list, loc_alpha_list, loc_beta_list, rate_alpha_list, rate_beta_list


def BGsampler_immi(observation, pi_0, sample_n, T_interval, k, mu, lamb, omega, theta):
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
    loc_alpha_list = []
    loc_beta_list = []
    rate_alpha_list = []
    rate_beta_list = []
    # Gibbs Iteration:
    for i in range(sample_n):
        sample_new, alpha, beta, loc_alpha, loc_beta, rate_alpha, rate_beta = BGsampler_one_A_immi(observation, copy.deepcopy(sample_old), k, mu, lamb, omega, theta)
        sample_list.append(copy.deepcopy(sample_new))
        sample_old = sample_new
        alpha_list.append(alpha)
        beta_list.append(beta)
        loc_alpha_list.append(loc_alpha)
        loc_beta_list.append(loc_beta)
        rate_alpha_list.append(rate_alpha)
        rate_beta_list.append(rate_beta)
    return sample_list, alpha_list, beta_list, loc_alpha_list, loc_beta_list, rate_alpha_list, rate_beta_list


