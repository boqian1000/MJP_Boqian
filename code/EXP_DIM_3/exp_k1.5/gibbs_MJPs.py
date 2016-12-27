from path_class import *
from math import log
#[t_0(t_start), t_1, t_2,...,t_N, t_N+1(t_end)]
# #[s_0, s_1, s_2,...,s_N]

# Input  : Old trajectory and bound constant OMEGA
# Output : Trajectory with true jumps and virtual jumps
def sampleUI(MJPpath0, OMEGA):
    # for more convenience calling
    A = MJPpath0.rate_matrix
#    t_start = MJPpath0.t_start
#    t_end = MJPpath0.t_end
    S = MJPpath0.S
    T = copy.deepcopy(MJPpath0.T)
    T.append(MJPpath0.t_end)
    T.insert(0, MJPpath0.t_start)
    AS = A.diagonal()
    vj_times = [] # virtual jump times
    vj_states = [] # virtual jump states
    #Here is the trick to generate non-homogeneous poisson process
    for i in range(len(T) - 1 ):
        current_state = S[i]
        rate = OMEGA + AS[current_state]
#        print rate
        t = float(T[i])
        z = 0.0 # time_step
        while t + z < T[i + 1]:
            t += z
            # add 0 as virtual jump time. delete it later
            vj_times.append(t) # in the end we will get the union of virtual jumps and effective jumps
            vj_states.append(current_state)
            scale = 1.0 / rate
            z = random.exponential(scale)
    new_MJP = copy.deepcopy(MJPpath0)
    vj_times.pop(0)
    new_MJP.S = vj_states
    new_MJP.T = vj_times
    return new_MJP


def get_likelihood(observation, path_T_, N, t_interval):
    t_start = t_interval[0]
    t_end = t_interval[1]
    path_T = copy.deepcopy(path_T_)
    path_T.append(t_end)
    path_T.insert(0, t_start)
    O = observation.O
    OT = observation.T
    likelihood_list = []
    No = len(OT)
    for i in range(len(path_T) - 1):
        likelihood = [1.0] * N
        for s in range(N):
            for j in range(No):
                if OT[j] < path_T[i + 1] and OT[j] >= path_T[i]:
                    likelihood[s] *= trans_Likelihood(O[j], s)
        likelihood_list.append(likelihood)
    # 0 , 1, 2, ..., n-1
    return likelihood_list

def FF(likelihood, initial_pi, OMEGA, path): #
    # X is the corresponding observations X1, ... ,XT
    N = len(initial_pi)
    rate_matrix = copy.deepcopy(path.rate_matrix)
    B = numpy.identity(N) + rate_matrix / float(OMEGA)
    alpha = [initial_pi]
    M = []
    M.append(1)
    # forward filtering:
    # (len(path.T) + 1) * N dimensional likelihood matrix
    for t in range(1, len(path.T) + 1):
        temp = [0.0] * N
        # j-th coordinate
        for j in range(N):
            for k in range(N):
                temp[j] += alpha[t - 1][k] * likelihood[t - 1][k] * B[k][j]
        maxt = max(temp)
        M.append(maxt)
        temp2 = [x / maxt for x in temp]
        alpha.append(temp2)
    # backward sampling:
    newS = []
    beta = array(likelihood[-1]) * array(alpha[-1])
    p_marginal = sum(beta)
    logm = 0.0
    for m in M:
        logm += log(m)
    log_p = log(p_marginal) + logm
    return log_p, alpha


def BS(initial_pi, alpha, likelihood, OMEGA, path): # Here path means path containing virtual jumps
    # X is the corresponding observations X1, ... ,XT
    N = len(initial_pi)
    rate_matrix = copy.deepcopy(path.rate_matrix)
    path_times = copy.deepcopy(path.T)
    B = numpy.identity(N) + rate_matrix / float(OMEGA)
    t_start = path.t_start
    t_end = path.t_end
    newS = []
    beta = array(likelihood[-1]) * array(alpha[-1])
    temp = sample_from_Multi(beta)
    newS.append(temp)
    iter = range(len(path_times))
    iter.reverse()
    for t in iter:
        beta = [0.0] * N
        for i in range(N):
            beta[i] = alpha[t][i] * likelihood[t][i] * B[i][newS[-1]]
        temp = sample_from_Multi(beta)
        newS.append(temp)
    newS.reverse()
    MJPpath_new= MJPpath(newS, path_times, t_start, t_end, rate_matrix, initial_pi)
    return MJPpath_new

def FFBS(initial_pi, observation, OMEGA, path): # Here path means path containing virtual jumps
    # X is the corresponding observations X1, ... ,XT
    rate_matrix = copy.deepcopy(path.rate_matrix)
#    X = observation.O
#    times = copy.deepcopy(observation.T)
    path_times = copy.deepcopy(path.T)
    B = numpy.identity(rate_matrix.shape[0]) + rate_matrix / float(OMEGA)
    N = len(initial_pi)
    t_start = path.t_start
    t_end = path.t_end
    alpha = [initial_pi]
    # forward filtering:
    likelihood = get_likelihood(observation, path.T, N, [t_start, t_end])
    # (len(path.T) + 1) * N dimensional likelihood matrix
    for t in range(1, len(path_times) + 1):
        temp = [0.0] * N
        # j-th coordinate
        for j in range(N):
            for k in range(N):
                temp[j] += alpha[t - 1][k] * likelihood[t - 1][k] * B[k][j] # error
        alpha.append(temp)
    # backward sampling:
    newS = []
    beta = array(likelihood[-1]) * array(alpha[-1])
    p_marginal = sum(beta)
    temp = sample_from_Multi(beta)
    newS.append(temp)
    iter = range(len(path_times))
    iter.reverse()
    for t in iter:
        beta = [0.0] * N
        for i in range(N):
            beta[i] = alpha[t][i] * likelihood[t][i] * B[i][newS[-1]]
        temp = sample_from_Multi(beta)
        newS.append(temp)
    newS.reverse()
    MJPpath_new = MJPpath(newS, path_times, t_start, t_end, rate_matrix, initial_pi)
    return MJPpath_new, p_marginal


# given rate matrix(i.e. all the parameters), and the trajectory with virtual jumps,
# sample a new trajectory.
def BGsampler_one(observation, MJPpath0, OMEGA):
    initial_pi = copy.deepcopy(MJPpath0.initial_pi)
    uipath = sampleUI(MJPpath0, OMEGA)
    likelihood = get_likelihood(observation, uipath.T, len(initial_pi), [uipath.t_start, uipath.t_end])
    # Forward calculate P(Y | W)
    p,  ALPHA = FF(likelihood, initial_pi, OMEGA, uipath)
    new_path = BS(initial_pi, ALPHA, likelihood, OMEGA, uipath)
    new_path.delete_virtual()
    return new_path


def BGsampler(observation, rate_matrix, pi_0, sample_n, T_interval, OMEGA):
    sample_old = MJPpath(t_start=T_interval[0], t_end=T_interval[1], rate_matrix=rate_matrix, initial_pi=pi_0)
    sample_old.generate_newpath()
    sample_list = []
    sample_list.append(copy.deepcopy(sample_old))
    for i in range(sample_n - 1):
        sample_new = BGsampler_one(observation, sample_old, OMEGA)
        sample_list.append(copy.deepcopy(sample_new))
        sample_old = sample_new
    return sample_list

