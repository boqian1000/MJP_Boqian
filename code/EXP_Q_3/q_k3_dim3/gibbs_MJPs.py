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
        t = float(T[i])
        z = 0.0 # time_step
        while t + z < T[i + 1]:
            t += z
            # add 0 as virtual jump time. delete it later
            vj_times.append(t) # in the end we will get the union of virtual jumps and effective jumps
            vj_states.append(current_state)
            z = random.exponential(1.0 / rate)
    new_MJP = copy.deepcopy(MJPpath0)
    vj_times.pop(0)
    new_MJP.S = vj_states
    new_MJP.T = vj_times
    return new_MJP

def get_likelihood(observation, path_T, N, t_interval):
    t_start = t_interval[0]
    t_end = t_interval[1]
    path_T = copy.deepcopy(path_T)
    path_T.append(t_end)
    path_T.insert(0, t_start)
    O = observation.O
    OT = observation.T
    likelihood_list = []
    for i in range(len(path_T) - 1):
        likelihood = [1] * N
        for s in range(N):
            for j in range(len(O)):
                if OT[j] < path_T[i + 1] and OT[j] >= path_T[i]:
                    likelihood[s] *= trans_Likelihood(O[j], s)
        likelihood_list.append(likelihood)
    # 0 , 1, 2, ..., n-1
    return likelihood_list

'''
    def llh_row(o, N):
        return map(trans_Likelihood, [o] * N, range(N))
    uncom_ll_list = map(llh_row, O, [N] * len(O))
    index = map(getindex, [path_T], OT)
    for i in range(len(index)):
    likelihood_list = [[1] * N] * (len(path_T) - 1)
    for i in range(len(likelihood_list)):

        for s in range(N):
            for j in range(len(O)):
                if OT[j] < path_T[i + 1] and OT[j] >= path_T[i]:
                    likelihood[s] *= trans_Likelihood(O[j], s)'''

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
            for k in range(max(j - 1, 0), min(j + 2, N)):
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
    N = len(initial_pi)
    t_start = path.t_start
    t_end = path.t_end
    likelihood = get_likelihood(observation, path.T, N, [t_start, t_end])
    log_p, alpha = FF(likelihood, initial_pi, OMEGA, path)
    MJPpath_new = BS(initial_pi, alpha, likelihood, OMEGA, path)
    return MJPpath_new


# given rate matrix(i.e. all the parameters), and the trajectory with virtual jumps,
# sample a new trajectory.
def BGsampler_one(observation, MJPpath0, OMEGA):
    #rate_matrix = copy.deepcopy(MJPpath0.rate_matrix)
    initial_pi = copy.deepcopy(MJPpath0.initial_pi)
    uipath = sampleUI(MJPpath0, OMEGA)
    #new_path = FFBS(rate_matrix,initial_pi, observation, OMEGA, current_path)
    new_path = FFBS(initial_pi, observation, OMEGA, uipath)
    new_path.delete_virtual()
    return new_path


def BGsampler(observation, rate_matrix, pi_0, sample_n, T_interval, OMEGA):
    sample_old = MJPpath(t_start=T_interval[0], t_end=T_interval[1], rate_matrix=rate_matrix, initial_pi=pi_0)
    sample_old.generate_newpath()
    sample_list = []
    sample_list.append(copy.deepcopy(sample_old))
    for i in range(sample_n - 1):
        sample_new = BGsampler_one(observation, copy.deepcopy(sample_old), OMEGA)
        sample_list.append(copy.deepcopy(sample_new))
        sample_old = sample_new
    return sample_list


'''
rate_matrix = numpy.array([[-1.0, 0.5, 0.5], [0.5, -1.0, 0.5], [0.5, 0.5, -1.0]], dtype='float')
pi_0 = [1.0/3.0, 1./3., 1./3.]
y = [2, 2, 2, 1, 2, 0]
T = [1, 2, 3, 4, 5, 6]
T_interval = [0, 10]
observation = Observation(y, T, 0, 10)
sample_n = 20000

test = BGsampler(observation, rate_matrix, pi_0, sample_n, T_interval, 2)
'''