import numpy as np
import copy
from config import *



#[t_start, t_1, t_2,...,t_N, t_end]
#[s_0, s_1, s_2,...,s_N]

class Observation:
    def __init__(self, O=None, T=None,t_start=0, t_end=10):
        """
        Constructor
        :param O:       list    observations
        :param T:       list    times when observing
        :param t_start: float
        :param t_end:   float
        :return:        An Observation instance
        """
        self.O = O
        self.T = T
        self.t_start = t_start
        self.t_end = t_end

    def sample_observation(self, time_points, MJPpath0):
        """
        Sample observations at given times, given an MJP path
        :param time_points: times at which to sample the observations
        :param MJPpath0:    MJP path
        :return:            Updated Observation instance
        """
        self.t_start = MJPpath0.t_start
        self.t_end = MJPpath0.t_end
        self.T = time_points
        S = copy.deepcopy(MJPpath0.S)
        mjpT = copy.deepcopy(MJPpath0.T)
        mjpT.insert(0, self.t_start)
        index = [getindex(mjpT, x) for x in time_points]
        States = [S[x - 1] for x in index]
        observation = random.normal(States)
        self.O = list(observation)

    def extend(self, T_list):
        """
        Add more observation times
        :param T_list:  list
        :return:    updated Observation instance
        """
        T = copy.deepcopy(self.T)
        O = copy.deepcopy(self.O)
        new_O = [O[getindex(T, t) - 1] for t in T_list]
        T.extend(T_list)
        O.extend(new_O)
        TO = zip(T,O)
        TO.sort()
        self.T = [x[0] for x in TO]
        self.O = [x[1] for x in TO]

    def info(self):
        print "O:", self.O
        print "T:", self.T
        print "t_start:", self.t_start
        print "t_end:", self.t_end


class MJPpath:
    def __init__(self,S=[0], T=[], t_start=0, t_end=10,
                 rate_matrix=np.array([[-1, 0.5, 0.5], [0.2, -0.7, 0.5], [0.8, 0.1, -0.9]], dtype='float'),
                 initial_pi=[0.1, 0.2, 0.7]):
        """
        Constructor
        :param S:   list of states with initial state
        :param T:   list of jump time without initial time
        :param t_start: starting time
        :param t_end:   ending time
        :param rate_matrix: rate matrix     np.array([[],[]]) 2D np.array matrix
        :param initial_pi: initial distribution pi_0
        :return:
        """
        self.s0 = S[0]
        self.S = S
        self.T = T
        self.t_start = t_start
        self.t_end = t_end
        self.total_state = len(initial_pi)
        self.rate_matrix = rate_matrix
        self.initial_pi = initial_pi
        if len(S) -1 != len(T) :
            print "ERROR:The length of S and the length of T are not matching!:("
    def backward_combine(self, path1):
        if self.t_start == path1.t_end and path1.S[-1] == self.s0:
            SS = copy.deepcopy(self.S)
            SS.pop(0)
            S = copy.deepcopy(path1.S)
            T = copy.deepcopy(path1.T)
            S.extend(SS)
            T.extend(copy.deepcopy(self.T))
            self.S = S
            self.s0 = path1.s0
            self.T = T
            self.t_start = path1.t_start
        else:
            print "NOT Matching :("

    def forward_combine(self,path1):
        if self.t_end == path1.t_start and self.S[-1] == path1.s0:
            SS = copy.deepcopy(path1.S)
            SS.pop(0)
            self.S.extend(SS)
            self.T.extend(copy.deepcopy(path1.T))
            self.t_end = path1.t_end

    def info(self):
        print "s0:", self.s0
        print "S:", self.S
        print "T:", self.T
        print "t_start:", self.t_start
        print "t_end:", self.t_end
        print "N:", self.total_state
        S = copy.deepcopy(self.S)
        T = copy.deepcopy(self.T)
        T.append(self.t_end)
        T.insert(0, self.t_start)
        S.append(S[-1])
        if len(self.T) ==1 and self.T[0] == -1:
            print "This is a prototype of MJP."
    #generate a new path given the rate matrix.
    def generate_newpath(self, known_first=False):
        t_start = self.t_start
        t_end = self.t_end
        t = t_start
        z = 0
        rate_matrix = copy.deepcopy(self.rate_matrix)
        if not known_first:
            pi0 = self.initial_pi
            s0 = sample_from_Multi(pi0)
            self.s0 = s0
        rate_list = abs(rate_matrix.diagonal())
        S = []
        T = []
        temp_state = self.s0
        while (t + z) < t_end:
            t += z
            current_state = temp_state
            S.append(current_state)
            if t > t_start:
                T.append(t)
            rate = rate_list[current_state]
            z = np.random.exponential(1.0 / rate)
            beta = rate_matrix[current_state]
            beta[beta < 0] = 0
            temp_state = sample_from_Multi(beta)
        self.S = S
        self.T = T

    def append(self, MJPpath0):
        if self.s0 == MJPpath0.s0 and self.t_start == MJPpath0.t_start and self.t_end == MJPpath0.t_end:
            S = copy.deepcopy(self.S)
            S = S[1: len(S)]
            z = copy.deepcopy(MJPpath0.S)
            S.extend(z[1: len(z)])
            T = copy.deepcopy(self.T)
            T.extend(copy.deepcopy(MJPpath0.T))
            # sorting
            TS = zip(T, S)
            TS = list(set(TS))
            TS.sort()
            self.S = [x[1] for x in TS]
            self.S.insert(0, self.s0)
            self.T = [x[0] for x in TS]
        else:
            print "Error: s0s, t_starts and t_ends are not all matching :("

    def append_virtual(self, Virtual_T):
        Virtual_S = [self.s0]
        time = copy.deepcopy(self.T)
        S = copy.deepcopy(self.S)
        for v in Virtual_T:
            index = getindex(time, v)
            Virtual_S.append(S[index])
        Virtualmjp = copy.deepcopy(self)
        Virtualmjp.S = Virtual_S
        Virtualmjp.T = Virtual_T
        Virtualmjp.info()
        self.append(Virtualmjp)

    def delete_virtual(self):
        S = copy.deepcopy(self.S)
        S = S[1: len(S)]
        T = copy.deepcopy(self.T)
        s0 = self.s0
        t_0 = self.t_start
        newS = []
        newT = []
        old_state = s0
        for i in range(len(S)):
            if S[i] != old_state:
                newS.append(S[i])
                newT.append(T[i])
                old_state = S[i]
        newS.insert(0, s0)
        self.S = newS
        self.T = newT

    # likelihood
    def likelihood(self):
        S = self.S
        T = copy.deepcopy(self.T)
        T.insert(0, self.t_start)
        A = self.rate_matrix
        s0 = S[0]
        product = self.initial_pi[s0]
        for i in range(1, len(T)):
            product *= A[S[i - 1]][S[i]] * np.exp(A[S[i - 1]][S[i - 1]] * (T[i] - T[i - 1]))
        product *= np.exp(A[S[-1]][S[-1]] * (self.t_end - T[-1]))
        return product
    # likelihood conditioned on s0

    def con_likelihood(self):
        return self.likelihood() / self.initial_pi[self.s0]

    def jump_times(self,s, s_to=-1):
        if s_to == -1:
            total_n = self.S.count(s)
            if self.S[-1] == s:
                total_n -= 1
            return total_n
        else:
            if s_to >= 0:
                total_n = 0
                s_from = s
                for i in range(len(self.S) - 1):
                    if self.S[i] == s_from and self.S[i + 1] == s_to:
                        total_n += 1
                return total_n
            else:
                print "only accept 2 parameters!"

    def stay_time(self,s):
        T = copy.deepcopy(self.T)
        T.insert(0, self.t_start)
        T.append(self.t_end)
        diff_T = [T[i + 1] - T[i] for i in range(len(T) - 1)]
#        diff_T = np.diff(array(T))
        staying_time = 0
        for i in range(len(diff_T)):
            if self.S[i] == s:
                staying_time += diff_T[i]
        return staying_time
#        return float(sum(array(diff_T)[np.array(self.S) == s]))

    def stay_time_list(self):
        list = []
        for s in range(self.total_state):
            list.append(self.stay_time(s))
        return list

    def get_virtual_jumps(self):
        V_S = []
        V_T = []
        V_index = []
        True_T = []
        True_S = [self.s0]
        for i in range(1 , len(self.S)):
            if self.S[i] == self.S[i - 1]:
                V_index.append(i)
                V_S.append(self.S[i])
                V_T.append(self.T[i -1])
            else:
                True_S.append(self.S[i])
                True_T.append(self.T[i -1])
        # Numbers of jump times in each true jumping time interval.
        tt = copy.deepcopy(True_T)
        tt.insert(0, self.t_start)
        tt.append(self.t_end)
        array_vt = np.array(V_T)
        V_Num_in_T_interval = []
        for i in range(len(tt) - 1):
            V_Num_in_T_interval.append(sum((array_vt<tt[i + 1]) * (array_vt > tt[i])))
        return V_index, V_T, V_S, True_S, True_T, V_Num_in_T_interval

#def p_trans(MJPpath, yn):
#    return trans_Likelihood(yn,MJPpath.S[-1])
