from __future__ import print_function
import sys
import os
sys.path.append(os.environ['FIASCO'])
import mripipes

file1= 'random1_stretch'
#file2= 'ident_qxts'
file2= 'random2'

a= mripipes.createArena();
r1= mripipes.createMRIFileInputTool(a,file1)
r2= mripipes.createMRIFileInputTool(a,file2)
t1= mripipes.createMatmultTool(a)
t1.setVerbose()
t1.setDebug()
out= mripipes.createMRIFileOutputTool(a,'prod')
t1.getSink(0).connect(r1.getSourceByName('images'))
t1.getSink(1).connect(r2.getSourceByName('images'))
out.getSink(0).connect(r1.getSourceByName('orphans'))
out.getSink(1).connect(t1.getSource(0))
if not a.init():
    print("Initialization failed!")
else:
    if not a.execute():
        print("Execution failed!")
    
