from __future__ import print_function
import sys
import os
sys.path.append(os.environ['FIASCO'])
import mripipes

t1= mripipes.createBaseTool()
s1= mripipes.createBaseDataSource();
s1.setName('images')
t1.addDataSource(s1)

s2= mripipes.createBaseDataSource();
s2.setName('missing')
t1.addDataSource(s2)

print(t1.getDataSource(0).getName())
print(t1.getDataSource(1).getName())
print(t1.getDataSource(0).getName())

mystery= t1.getSourceByName('missing')
print(mystery.getName())
