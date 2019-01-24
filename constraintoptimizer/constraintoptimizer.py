import os
import sys
import numpy as np

class TConstraintOptimizer():
    def __init__(self,N):
        self.N = N
        self.Ncon = int(self.N * (self.N - 1) / 2)
        self.sv = self._fastsingleviolators(self.N)

    def _conftolhzconf(self, conf):
        lhzconf = []
        nconf = [1]
        nconf.extend(conf)
        N = len(nconf)

        for bl in range(1, N):
            for bi in range(N-bl):
                lhzconf.append(nconf[bi]*nconf[bi+bl])
        return lhzconf

    def _decomposejij(self,jijmatrix):
        jij = []
        jijlist = []

        hi = jijmatrix[0][:]
        for ii in range(self.N - 1):
            jij.append(jijmatrix[ii + 1][1:])

        for bl in range(1, self.N + 1):
            for bi in range(self.N + 1 - bl):
                jijlist.append(jijmatrix[bi][bi + bl - 1])

        return hi, jij, jijlist

    def getenergydirect(self,hi,jij,conf):
        e = 0
        for ii in range(self.N):
            e += conf[ii]*hi[ii]
        for bl in range(self.N-1):
            for bi in range(self.N-1-bl):
                e += conf[bi] * conf[bi+bl+1] * jij[bi][bi+bl]
        return e

    def inttoconf(self,x,NN):
        return [-1 if int(d) == 0 else 1 for d in str(bin(x))[2:].zfill(NN)]

    def getlowestdirect(self,hi,jij,NN):
        eh = []

        for x in range(2 ** NN):
            conf = self.inttoconf(x,NN)
            eh.append(self.getenergydirect(hi,jij,conf))

        sl = sorted(eh)

        return sl

    def singledown(self,Ncon,jijlist):
        emins = []
        for jo in range(Ncon):
            ar = []
            for lhzconf in self.sv[jo]:
                hrest = np.dot(jijlist,lhzconf)
                ar.append(hrest)
            emins.append(min(ar))
        return(np.array(emins))

    def cindextocoordinates(self,N):
        coordinates = []
        for bl in range(1, N):
            for bi in range(N-bl):
                coordinates.append([bi+1,bi+bl])
        return coordinates

    def listify(self,N):
        lx= []
        ly = []
        for bl in range(1, N+1):
            for bi in range(N+1-bl):
                lx.append(bi)
                ly.append(bi+bl-1)
        return np.array(lx),np.array(ly)

    def _fastsingleviolators(self, N):

        listification = self.listify(N)
        coord = self.cindextocoordinates(N)
        singleviolators = []
        alllhzconf = []
        for x in range(2**N):
            conf = self.inttoconf(x, N)
            lhzconf = self._conftolhzconf(conf)
            alllhzconf.append(np.array(lhzconf))

        for coni,item in enumerate(coord):
            singleviolators.append([])
            mt = np.ones((N, N))
            mt[:item[0],item[1]:] *= -1
            mtlist = mt[listification]
            for lhzconf in alllhzconf:
                cf = np.multiply(lhzconf,mtlist)
                singleviolators[coni].append(cf)

        return singleviolators

    def getconstraints(self,upto,jijTotal):

        hi, jij, jijlist = self._decomposejij(jijTotal)

        denergy = self.getlowestdirect(hi, jij, self.N)
        emins = self.singledown(self.Ncon,jijlist)
        fs = -1.0 * emins + denergy[upto-1]
        return fs


