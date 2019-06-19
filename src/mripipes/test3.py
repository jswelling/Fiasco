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
t2= mripipes.createPassthruTool(a)
t2.getSink(0).connect(t1.getSourceByName('factor1'))
t3= mripipes.createDevnullTool(a)
t3.getSink(0).connect(t2.getSource(0))
a.init();
a.execute()

