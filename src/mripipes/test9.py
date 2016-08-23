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

a= mripipes.createArena()
infile= mripipes.createMRIFileInputTool(a,'prod')
infile.setDebug()
out= mripipes.createMRIFileOutputTool(a,'prod2')
out.setDebug()
#t1= mripipes.createPassthruTool(a)
t1= mripipes.createSubsetTool(a,"y",5,3)
t1.setDebug()
#t1= mripipes.createSubsetTool(a,"z",1,0)
#t1= mripipes.createSubsetTool(a,"x",2,0)
out.getSink(0).connect(infile.getSourceByName('orphans'))
t1.getSink(0).connect(infile.getSourceByName('images'))
out.getSink(1).connect(t1.getSource(0))

if not a.init():
    print "Initialization failed!"
else:
    print "Initialization OK!"
    if not a.execute():
        print "Execution failed!"
    
