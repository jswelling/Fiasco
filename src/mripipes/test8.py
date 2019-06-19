from __future__ import print_function
#! /usr/bin/env python
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
r1= mripipes.createZeroSrcTool(a,"xyz","64:64:4")

t1= mripipes.createRpnMathTool(a,'$x,100,*,$y,+');
t1.setVerbose()

out= mripipes.createMRIFileOutputTool(a,'prod')
out.setVerbose()
out.setDebug()
out.getSink(0).connect(t1.getSource(0))
t1.getSink(0).connect(r1.getSourceByName('images'))
if not a.init():
    print("Initialization failed!")
else:
    print("Initialization OK!")
    if not a.execute():
        print("Execution failed!")
    
