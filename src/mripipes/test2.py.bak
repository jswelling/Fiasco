import sys
import os
sys.path.append(os.environ['FIASCO'])
import mripipes

t1= mripipes.createMRIFileInputTool('test')
nSources= t1.getNDataSources()
for i in range(nSources):
    src= t1.getDataSource(i)
    print src.getName()
t2= mripipes.createDevnullTool()
t2.getDataSink(0).connect(t1.getSourceByName('factor1'))
t1.init()
t2.init()
t2.execute()

