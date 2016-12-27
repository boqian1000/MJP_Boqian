from path_class import *
from gibbs_MJPs import BGsampler_one
from math import log
from model_parameters import constructor_rate_matrix
from model_parameters import my_print
from model_parameters import get_omega
#[t_start, t_1, t_2,...,t_N, t_end]
# #[s_0, s_1, s_2,...,s_N]
#input: alpha(scale par), beta , dimension (d)
#output: immigration model transition matrix


def BGsampler_one_A(observation, MJPpath0, k, mu, lamb, alpha):
    #initialization
    OMEGA = get_omega(alpha, k)
    new_path = BGsampler_one(observation, MJPpath0, OMEGA)
    # Sample alpha | S, T, y:
    n = len(new_path.T)
    rate_new = lamb + (new_path.t_end - new_path.t_start) * 3.
    alpha_new = numpy.random.gamma(mu + n, 1./ rate_new)
    new_path.rate_matrix = constructor_rate_matrix(alpha_new)
    return new_path, alpha_new

def BGsampler_jc(observation, pi_0, sample_n, T_interval, k, mu, lamb):
    # Initialization:
    alpha_old = 2.0
    rate_matrix = constructor_rate_matrix(alpha_old)
    sample_old = MJPpath(t_start=T_interval[0], t_end=T_interval[1], rate_matrix=rate_matrix, initial_pi=pi_0)
    sample_old.generate_newpath()
    sample_list = []
    sample_list.append(copy.deepcopy(sample_old))
    alpha_list = []
    # Gibbs Iteration:
    for i in range(sample_n):
        sample_new, alpha_new = BGsampler_one_A(observation, sample_old, k, mu, lamb, alpha_old)
        sample_list.append(copy.deepcopy(sample_new))
        alpha_list.append(alpha_new)
        sample_old, alpha_old = sample_new, alpha_new
    return sample_list, alpha_list

def Full_GBS(timebreaks, pi_0, sample_n, T_interval, k, mu, lamb):
    observation = Observation()
    alpha_list = []
    ST_list = []
    alpha_old = 2.0
    rate_matrix = constructor_rate_matrix(alpha_old)
    sample_old = MJPpath(t_start=T_interval[0], t_end=T_interval[1], rate_matrix=rate_matrix, initial_pi=pi_0)
    sample_old.generate_newpath()
    for i in range(sample_n):
        my_print(i, 1000)
        observation.sample_observation(timebreaks, sample_old)
        sample_new, alpha_new = BGsampler_one_A(observation, sample_old, k, mu, lamb, alpha_old)
        ST_list.append(copy.deepcopy(sample_new))
        alpha_list.append(alpha_new)
        sample_old, alpha_old = sample_new, alpha_new
    return ST_list, alpha_list
