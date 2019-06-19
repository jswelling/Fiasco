from __future__ import print_function
import sys
import os
sys.path.append(os.environ['FIASCO'])
import mripipes

a= mripipes.createArena();
t1= mripipes.createMRIFileInputTool(a,'test')
nSources= t1.getNSources()
for i in range(nSources):
    src= t1.getSource(i)
    print(src.getName())
t2= mripipes.createMRIFileOutputTool(a,"junk")
t3= mripipes.createPassthruTool(a);
t2.getSink(0).connect(t1.getSourceByName('factor1'))
t3.getSink(0).connect(t1.getSourceByName('factor2'))
t2.getSink(1).connect(t3.getSource(0))
t2.getSink(2).connect(t1.getSourceByName('orphans'))
a.init();
a.execute()

