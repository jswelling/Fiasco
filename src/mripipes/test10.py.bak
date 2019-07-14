#! /bin/env python
import sys
import os
sys.path.append(os.environ['FIASCO'])
import mripipes

#file1= 'crand1'
#file1= 'crand1_stretch2'
#file1= 'test1_paste'
#file2= 'crand2'
file1= 'crand1_stretch'
file2= 'crand3_stretch'

a= mripipes.createArena();
r1= mripipes.createZeroSrcTool(a,"xyz","3:9:4")

t1= mripipes.createRpnMathTool(a,'$x,100,*,$y,10,*,+,$z,+,1000,+,dup,1,if_print_1')
t1.setDebug()

t2= mripipes.createSubsetTool(a,"y",5,0)
#t2= mripipes.createSubsetTool(a,"z",1,0)
#t2= mripipes.createSubsetTool(a,"x",2,0)
#t2= mripipes.createPassthruTool(a)
t2.setDebug()
t2.setVerbose()
t2.setDebug()

out= mripipes.createMRIFileOutputTool(a,'prod')
out.setDebug()
out.getSink(0).connect(t2.getSource(0))
t2.getSink(0).connect(t1.getSource(0))
t1.getSink(0).connect(r1.getSourceByName('images'))
if not a.init():
    print "Initialization failed!"
else:
    print "Initialization OK!"
    if not a.execute():
        print "Execution failed!"
    
