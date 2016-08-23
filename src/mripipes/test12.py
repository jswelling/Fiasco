#! /bin/env python
import sys
import os
sys.path.append(os.environ['FIASCO'])
import mripipes

def myCallback(size,offset):
#    print "My God, I made it! %d %Ld"%(size,offset)
    return (17,31)
#    return (17,"hello")

#file1= 'crand1'
#file1= 'crand1_stretch2'
#file1= 'test1_paste'
#file2= 'crand2'
file1= 'crand1_stretch'
file2= 'crand3_stretch'

a= mripipes.createArena();
r1= mripipes.createZeroSrcTool(a,"xyz","3:9:4")

t2= mripipes.createBlockMapTool(a,"y",5,0,myCallback)
#t2= mripipes.createBlockMapTool(a,"z",1,0,myCallback)
#t2= mripipes.createBlockMapTool(a,"x",2,0,myCallback)
#t2= mripipes.createPassthruTool(a)
t2.setDebug()
t2.setVerbose()

out= mripipes.createMRIFileOutputTool(a,'junk')
out.setDebug()
out.getSink(0).connect(t2.getSource(0))
t2.getSink(0).connect(r1.getSourceByName('images'))
if not a.init():
    print "Initialization failed!"
else:
    print "Initialization OK!"
    if not a.execute():
        print "Execution failed!"
    
