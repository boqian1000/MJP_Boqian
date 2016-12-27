from path_class import *
from gibbs_MJPs import sampleUI
from gibbs_MJPs import get_likelihood
from gibbs_MJPs import FFBS
from gibbs_MJPs import BGsampler_one
from model_parameters import constructor_rate_matrix
from model_parameters import my_print
#[t_start, t_1, t_2,...,t_N, t_end]
# #[s_0, s_1, s_2,...,s_N]
#input: alpha, beta , dimension (d)
#output: immigration model transition matrix


def BGsampler_one_A_immi(observation, MJPpath0, k, mu, lamb, omega, theta):
    #initialization
    OMEGA = k * max(-numpy.diag(MJPpath0.rate_matrix))
    new_path = BGsampler_one(observation, MJPpath0, OMEGA)
    # Sample alpha & beta| S, T, y:
    move = numpy.diff(new_path.S)
    Tau = new_path.stay_time_list()
    N = new_path.total_state
    alpha = numpy.random.gamma(mu + sum(move > 0), 1./ (lamb + MJPpath0.t_end - MJPpath0.t_start - Tau[-1]))
    beta = numpy.random.gamma(omega + sum(move < 0), 1./ (theta + numpy.inner(Tau ,range(N))))
    A_new = constructor_rate_matrix(alpha, beta, N)
    new_path.rate_matrix = copy.deepcopy(A_new)
    return new_path, alpha, beta

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
    # Gibbs Iteration:
    for i in range(sample_n):
        sample_new, alpha, beta = BGsampler_one_A_immi(observation, copy.deepcopy(sample_old), k, mu, lamb, omega, theta)
        sample_list.append(copy.deepcopy(sample_new))
        sample_old = sample_new
        alpha_list.append(alpha)
        beta_list.append(beta)
    return sample_list, alpha_list, beta_list


