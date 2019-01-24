import constraintoptimizer as co
import random



def generateJij(N):
    Jij = []
    JijTotal = []
    Jijlist = []

    for ii in range(N):
        JijTotal.append([0] * (N))

    for ii in range(N):
        for jj in range(ii+1):
            JijTotal[jj][ii] = (random.random()-0.5)

    hi = JijTotal[0][:]
    for ii in range(N-1):
        Jij.append(JijTotal[ii+1][1:])

    for bl in range(1, N+1):
        for bi in range(N+1-bl):
            Jijlist.append(JijTotal[bi][bi+bl-1])

    return hi,Jij,JijTotal,Jijlist



def generateJijGauss(N):
    Jij = []
    JijTotal = []
    Jijlist = []

    for ii in range(N):
        JijTotal.append([0] * (N))

    for ii in range(N):
        for jj in range(ii+1):
            JijTotal[jj][ii] = (random.gauss(0,1))

    hi = JijTotal[0][:]
    for ii in range(N-1):
        Jij.append(JijTotal[ii+1][1:])

    for bl in range(1, N+1):
        for bi in range(N+1-bl):
            Jijlist.append(JijTotal[bi][bi+bl-1])

    return hi,Jij,JijTotal,Jijlist



if __name__=="__main__":

    NN = 4
    fileN = 0
    upto = 1

    coo = co.TConstraintOptimizer(NN)

    hi, jij, jijtotal, jijlist = generateJij(NN)

    fs = coo.getconstraints(upto,jijtotal)
    print(jijtotal)
    print(fs)

