#! /bin/env python
import os,sys,numpy,random
sys.path.append(os.environ['FIASCO'])
import fiasco_numpy

def basicTests():
    L= 3
    M= 4

    print "##############################################"
    print "######                                ########"
    print "######          Basic Tests           ########"
    print "######                                ########"
    print "##############################################"
    print "###### Beginning tests of KalmanState ########"
    state= fiasco_numpy.createKalmanState(M)

    print "Getting L and M"
    if state.getM()!=M: print "$$$$$ Error: M value is invalid!"

    print "##### Dumping uninitialized KalmanState ######"
    state.dumpSelf(sys.stdout)

    print "Setting fields"
    x= numpy.array([1.0,1.0,1.0,1.0])
    state.setX(x)
    P= numpy.matrix([[0.25,0.0,0.0,0.0],
                     [0.0,1.0,0.5,0.0],
                     [0.0,0.5,1.0,0.0],
                     [0.0,0.0,0.0,0.75]])
    state.setP(P)

    print "##### Dumping initialized KalmanState ######"
    state.dumpSelf(sys.stdout)

    print "##############################################"
    print "###### Beginning tests of KalmanProcess ######"
    print "##############################################"
    process= fiasco_numpy.createKalmanProcess(L,M)

    print "##### Dumping uninitialized KalmanState ######"
    process.dumpSelf(sys.stdout)

    print "Setting and getting debug"
    if process.getDebug()!=0: print "$$$$$ Error: debug not initialized to zero!"
    process.setDebug(1)
    if process.getDebug()!=1: print "$$$$$ Error: debug not settable!"
    process.setDebug(0)

    print "Getting L and M"
    if process.getL()!=L: print "$$$$$ Error: L value is invalid!"
    if process.getM()!=M: print "$$$$$ Error: M value is invalid!"

    print "#### Testing set/get #####"
    print "Setting and getting A"
##     A= numpy.matrix([[1,0.2,0.3,0.4],
##                      [0.1,1.1,0.31,0.41],
##                      [0.12,0.22,0.72,0.42],
##                      [0.13,0.23,-0.33,1.43]])
    A_T= numpy.matrix([[1,0.2,0.3,0.4],
                     [0.1,1.1,0.31,0.41],
                     [0.12,0.22,0.72,0.42],
                     [0.13,0.23,-0.33,1.43]])
    A= A_T.transpose()
    process.setA(A)
    AOut= numpy.zeros([4,4])
    process.getA(AOut)
    if (A!=AOut).any():
        print "$$$$$ Error: failure setting and getting A!"
        print "A:"
        print A
        print "AOut:"
        print AOut

    print "Setting and getting H"
##     H= numpy.matrix([[0.7,0.8,0.9,1.0],
##                      [0.71,0.81,0.91,1.01],
##                      [0.72,0.82,0.92,1.02]])
    H= numpy.matrix([[0.7,1.0,0.91,0.82],
                     [0.8,0.71,1.01,0.92],
                     [0.9,0.81,0.72,1.02]])
    process.setH(H)
    HOut= numpy.zeros([3,4])
    process.getH(HOut)
    if (H!=HOut).any():
        print "$$$$$ Error: failure setting and getting H!"
        print "H:"
        print H
        print "HOut:"
        print HOut

    print "Setting and getting Q"
    Q= numpy.matrix([[0.25,0.0,0.0,0.0],
                     [0.0,1.0,0.5,0.0],
                     [0.0,0.5,1.0,0.0],
                     [0.0,0.0,0.0,0.75]])
    process.setQ(Q)
    QOut= numpy.zeros([4,4])
    process.getQ(QOut)
    if (Q!=QOut).any():
        print "$$$$$ Error: failure setting and getting Q!"
        print "Q:"
        print Q
        print "QOut:"
        print QOut

    print "Setting and getting R"
    R= numpy.matrix([[0.025,0.0,0.0],
                     [0.0,0.1,0.035],
                     [0.0,0.035,0.1]])
    process.setR(R)
    ROut= numpy.zeros([3,3])
    process.getR(ROut)
    if (R!=ROut).any():
        print "$$$$$ Error: failure setting and getting R!"
        print "R:"
        print R
        print "ROut:"
        print ROut

    print "##### Dumping initialized KalmanProcess ######"
    process.dumpSelf(sys.stdout)

    print "##### Testing one Kalman iteration ######"
    print "Initial state follows"
    state.dumpSelf(sys.stdout)
    z= numpy.array([0.7,0.8,0.9])
    #process.setDebug(1)
    logLikelihoodDelta= process.apply(state, z, 0, 1, 1)
    print "logLikelihoodDelta= %g"%logLikelihoodDelta
    print "Updated state follows"
    state.dumpSelf(sys.stdout)

###############
# This test follows the example in Welch & Bishop's "An Introduction
# to the Kalman Filter".  Convergence is a little slower than
# is demonstrated there because our initial P (initialized to match Q)
# starts out inconveniently small.
###############
def estRandomConstant():
    print "##############################################"
    print "######                                ########"
    print "######     Testing KalmanFilter by    ########"
    print "######  estimating a random constant  ########"
    print "######                                ########"
    print "##############################################"
    valToEstimate= -0.37727
    sigma= 0.1
    nsteps= 200
    A= numpy.matrix([[1.0]])
    H= numpy.matrix([[1.0]])
    Q= numpy.matrix([[1.0e-5]])
    R= numpy.matrix([[sigma*sigma]])
    proc= fiasco_numpy.createKalmanProcess(1,1)
    proc.setA(A)
    proc.setH(H)
    proc.setQ(Q)
    proc.setR(R)
##    proc.setDebug(1)
    filter= fiasco_numpy.createKalmanFilter(proc)
    print "Filter follows:"
    filter.dumpSelf(sys.stdout)
    zIn= numpy.zeros([nsteps,1])
    random.seed(123) # so every run generates the same sequence
    for i in xrange(nsteps):
        zIn[i,0]= valToEstimate+random.normalvariate(0.0,sigma)
    xOut= numpy.zeros([nsteps,1])
    POut= numpy.zeros([nsteps,1,1])
    logLikelihood= filter.run(nsteps,zIn,xOut,POut,1)
    print "Final log likelihood: %g"%logLikelihood
    print "Last 20 Z's follow"
    print zIn[-20:]
    print "Last 20 X's follow"
    print xOut[-20:]
    print "Last 20 P's follow"
    print POut[-20:]

###############
# This test is like estRandomConstant(), but with *2* constants, to
# avoid falling into the H=1 special case in the kalmanfilter.c code.
###############
def estTwoRandomConstants():
    print "##############################################"
    print "######                                ########"
    print "######     Testing KalmanFilter by    ########"
    print "######  estimating 2 random constants ########"
    print "######                                ########"
    print "##############################################"
    val1= -0.37727
    val2= 0.213752
    sigma1= 0.1
    sigma2= 0.2
    nsteps= 500
    A= numpy.matrix([[1.0, 0.0],
                     [0.0, 1.0]])
    H= numpy.matrix([[1.0, 0.0],
                     [0.0, 1.0]])
    Q= numpy.matrix([[1.0e-5, 0.0],
                     [0.0, 1.0e-5]])
    R= numpy.matrix([[sigma1*sigma1, 0.0],
                     [0.0, sigma2*sigma2]])
    proc= fiasco_numpy.createKalmanProcess(2,2)
    proc.setA(A)
    proc.setH(H)
    proc.setQ(Q)
    proc.setR(R)
##     proc.setDebug(1)
    filter= fiasco_numpy.createKalmanFilter(proc)
    print "Filter follows:"
    filter.dumpSelf(sys.stdout)
    zIn= numpy.zeros([nsteps,2])
    random.seed(123) # so every run generates the same sequence
    for i in xrange(nsteps):
        zIn[i,0]= val1+random.normalvariate(0.0,sigma1)
        zIn[i,1]= val2+random.normalvariate(0.0,sigma2)
    xOut= numpy.zeros([nsteps,2])
    POut= numpy.zeros([nsteps,2,2])
    logLikelihood= filter.run(nsteps,zIn,xOut,POut,1)
    print "Final log likelihood: %g"%logLikelihood
##     print "Z follows"
##     print zIn
    print "Final X follows"
    print xOut[-1]
    print "final P follows"
    print POut[-1]

###############
# Main
###############

basicTests()

estRandomConstant()

estTwoRandomConstants()
