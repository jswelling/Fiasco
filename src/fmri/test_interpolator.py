#! /bin/env python
import os,sys,numpy
sys.path.append(os.environ['FIASCO'])
import fiasco_numpy

#intrp= fiasco_numpy.createInterpolator3DByType(fiasco_numpy.INTRP_LINEAR,
#                                           5,7,11,3)
intrp= fiasco_numpy.createInterpolator1DByType(fiasco_numpy.INTRP_LINEAR,5,3)
print("typename is <%s>"%intrp.getTypeName())
print("Debug is %d"%intrp.getInt(fiasco_numpy.INTRP_OPT_DEBUG))
intrp.setInt(fiasco_numpy.INTRP_OPT_DEBUG,1)
print("Now debug is %d"%intrp.getInt(fiasco_numpy.INTRP_OPT_DEBUG))
# print "Dumping now"
# intrp.dumpSelf(sys.stdout)
inData= numpy.zeros([11,7,5,3])
for v in range(3):
    for i in range(5):
        for j in range(7):
            for k in range(11):
                inData[k,j,i,v]= 1000.0*k+100.0*j+10*i+v
inData2= inData[1,1,:,:]
intrp.prep(inData2)
#intrp.prep(inData)
outData= numpy.zeros([3])
print("Output field before: %s"%outData)
intrp.calc(outData,(2.5,0.2,3.1),3,0)
print("Output field after: %s"%outData)

# Should really test warpApply as well!
