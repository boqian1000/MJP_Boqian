import numpy as np
from path_class import *
import copy
def SampleMMPP(MJPpath, lambs):
    NumJumps = len(MJPpath.T)
    grid = copy.deepcopy(MJPpath.T)
    grid.insert(0, MJPpath.t_start)
    grid.append(MJPpath.t_end)
    MMPPTime = []
    CurInterval = 0
    while CurInterval <= NumJumps:
        rate = lambs[MJPpath.S[CurInterval]]
        deltaT = grid[CurInterval + 1] - grid[CurInterval]
        Npoisson = np.random.poisson(deltaT * rate)
        mmpps = np.random.uniform(grid[CurInterval], grid[CurInterval + 1], Npoisson)
        mmpps.sort()
        mmpps = list(mmpps)
        MMPPTime.extend(mmpps)
        CurInterval += 1
    obs = Observation(O = [1] * len(MMPPTime), T = MMPPTime, t_start= MJPpath.t_start, t_end= MJPpath.t_end)
    return obs

