from path_class import *
from config import *

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


# model dependent
def get_likelihood(observation, path_T_, t_interval, lambs):
    t_start = t_interval[0]
    t_end = t_interval[1]
    grid = copy.deepcopy(path_T_)
    grid.append(t_end)
    grid.insert(0, t_start)
    Obs = observation.O
    ObsTime = observation.T
    likelihood_list = []
    NumofObs = len(ObsTime)
    CurItvIndex = 0
    for i in xrange(len(grid) - 1):
        likelihood = np.ones(2)
        obscount = 0
        deltaT = grid[i + 1] - grid[i]
        while ObsTime[CurItvIndex] < grid[i + 1]:
            obscount += 1
            CurItvIndex += 1
        likelihood *= np.array([(lambs[0] ** obscount) * np.exp(-lambs[0] * deltaT), (lambs[1] ** obscount) * np.exp(-lambs[1] * deltaT)]
        likelihood_list.append(likelihood)
    return likelihood_list

# given rate matrix(i.e. all the parameters), and the trajectory with virtual jumps,
# sample a new trajectory.
def BGsampler_one(observation, MJPpath0, OMEGA, lambs):
    initial_pi = copy.deepcopy(MJPpath0.initial_pi)
    uipath = sampleUI(MJPpath0, OMEGA)
    likelihood = get_likelihood(observation, uipath.T, [uipath.t_start, uipath.t_end], lambs)
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
