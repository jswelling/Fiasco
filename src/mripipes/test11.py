from __future__ import print_function
#! /bin/env python
import sys
import os
sys.path.append(os.environ['FIASCO'])
import mripipes

def initCallback(dim,dimstr,fastBlk,upstream_extent,extent1,extent2,slowBlk,hookDict):
    print("initCallback Args: %s <%s> %d %d %d %d %d %s"%\
          (dim,dimstr,fastBlk,upstream_extent,extent1,extent2,slowBlk,hookDict))
    hookDict['fastBlk']=fastBlk*extent1
    hookDict['slowBlk']=slowBlk
    hookDict['extent']= extent2
    if (len(hookDict['blockOffsetList'])%extent2 != 0):
        print("Length of block offset list does not match extent!")
        return 0
    return 1

def myCallback(size,offset,hookDict):
    blockOffsetList= hookDict['blockOffsetList']
    extent= hookDict['extent']
    fastBlk= hookDict['fastBlk']
    n_fast_blks= offset/fastBlk
    fast_blk_offset= offset%fastBlk
    n_full_extents= n_fast_blks/extent
    extent_offset= n_fast_blks%extent
##     print "callback: %d %d %d %d"%\
##           (n_fast_blks,fast_blk_offset,n_full_extents,extent_offset)
    return (fastBlk-fast_blk_offset,
            blockOffsetList[extent_offset]+fast_blk_offset)

#file1= 'crand1'
#file1= 'crand1_stretch2'
#file1= 'test1_paste'
#file2= 'crand2'
file1= 'crand1_stretch'
file2= 'crand3_stretch'

hookDict= {'blockOffsetList':[1,2,3,4]}

a= mripipes.createArena();
r1= mripipes.createZeroSrcTool(a,"xyz","3:9:4")

t2= mripipes.createBlockMapTool(a,"y","q",5,4,(initCallback,myCallback,hookDict))
t2.setDebug()
t2.setVerbose()

out= mripipes.createMRIFileOutputTool(a,'junk')
out.setDebug()
out.setVerbose()
out.getSink(0).connect(t2.getSource(0))
t2.getSink(0).connect(r1.getSourceByName('images'))
a.init()
a.execute()
    
