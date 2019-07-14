#! /usr/bin/env python
import os,sys,numpy
sys.path.append(os.environ['FIASCO'])
import fiasco_numpy

opt= fiasco_numpy.createNelminOptimizer(0.0001,1.0,5);

print "str is <%s>"%opt
print "repr is <%s>"%repr(opt)

def valFun( array, hook ):
    #print "In valFun! <%s> <%s>"%(repr(array),repr(hook))
    return (1.27*array[0] + 3.14*array[0]*array[0] \
            + 17.3*array[1] + 12.0 * array[1]*array[1] \
            + 3.2*array[2] + 3.7*array[2]*array[2] \
            + 17.3)

def resetFun( array, hook ):
    print "In resetFun! <%s> <%s>"%(repr(array),repr(hook))
    pass

sf= fiasco_numpy.buildSimpleScalarFunction( (valFun, resetFun, 3, "hello") )

valArray= numpy.array( [101.0,102.0,103.0] )
print "value: %g"%sf.value(valArray)
print "reset follows"
sf.reset(valArray)

opt.setDebugLevel(3)
r= opt.go(sf, valArray)
print "Optimization result is <%s>"%r
