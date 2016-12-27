from path_class import *


def IS_density(path, yn ,path_old):
    return path.con_likelihood


def cal_weight(yn, path):
    return trans_Likelihood(yn, path.S[-1])


# For proposal step in SCM
def generate_particles(N, rate_matrix, pi_0, T_interval, s0_list, known_first=False):
    if not known_first:
        s0_list = numpy.random.choice(len(pi_0), N, p=pi_0)
        #s0_list = map(sample_from_Multi,[pi_0] * N)
    particles = []
    for i in range(N):
        temp_particle = MJPpath(t_start=T_interval[0], t_end=T_interval[1], S = [s0_list[i]], rate_matrix=rate_matrix, initial_pi=pi_0)
        temp_particle.generate_newpath(known_first)
        particles.append(copy.deepcopy(temp_particle))
    return particles


def particles_combine(particles, particles_old, A_index):
    for i in range(len(A_index)):
        particles[i].backward_combine(particles_old[A_index[i]])

def SMC_MJPs(y, T, rate_matrix, pi_0, T_interval, N):
    #step 1:
    TT = copy.deepcopy(T)
    TT.append(T_interval[1])
    m_likelihood_list = []
    #particles_list = []
    particles = generate_particles(N, rate_matrix, pi_0, [T_interval[0], T[0]],[-1] * N, False)
    old_particles = particles
    weights_list = []
    weights = map(cal_weight, [y[0]] * N, particles)
    A_list = []
    A = sample_from_Multi(weights, N)
    A_list.append(A)
    s0_list = [particles[x].S[-1] for x in A]
    m_likelihood = numpy.mean(weights)
    m_likelihood_list.append(m_likelihood)
    for i in range(1, len(TT)):
        particles = generate_particles(N, rate_matrix, pi_0, [TT[i - 1], TT[i]], s0_list, True)
        if i < len(TT) - 1:
            weights = map(cal_weight, [y[i]] * N, particles)
        else:
            weights = [1.0 / N] * N
        particles_combine(particles, old_particles, A_list[-1])
        old_particles = particles
        A = sample_from_Multi(weights, N)
        A_list.append(A)
        s0_list = [particles[x].S[-1] for x in A]
        temp_p = [weights[a] for a in A]
        m_likelihood_list.append(numpy.mean(temp_p))
    # combine the last one
    return particles, weights, numpy.prod(m_likelihood_list) / m_likelihood_list[-1]



def sample_from_particles(particles, weights):
    a = sample_from_Multi(weights)
    return particles[a]

def PIMHsampler(observation, rate_matrix, pi_0, N, sample_n, T_interval):
    y = observation.O
    T = observation.T
    MJP_samples = []
    # Step 1
    particles, weights, old_p = SMC_MJPs(y, T, rate_matrix, pi_0, T_interval, N)
    old_sample = sample_from_particles(particles, weights)
    MJP_samples.append(copy.deepcopy(old_sample))
    for iter in range(sample_n):
#        if iter % 100 == 0:
#            print iter
#        particles_list, weights_list, new_p = SMC_MJPs(y, T, rate_matrix, pi_0, T_interval, N)
        particles, weights, new_p = SMC_MJPs(y, T, rate_matrix, pi_0, T_interval, N)
        new_sample = sample_from_particles(particles, weights)
        accept_rate = min(1, new_p / old_p)
        if random.uniform() < accept_rate:
            old_sample = new_sample
            old_p = new_p
        MJP_samples.append(copy.deepcopy(old_sample))
    return MJP_samples
