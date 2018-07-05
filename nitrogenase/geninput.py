import numpy
import itertools

def init():
   norb = 76
   # Classifier
   orbs = [0]*13
   orbs[0]  = [0,1]				# _end
   orbs[1]  = [2,3,4,5,6]                       # _fe
   orbs[2]  = [7,8,9,10,11,12,13,14,15]	        # _s 
   orbs[3]  = [16,17,18,19,20]                  # _fe
   orbs[4]  = [21,22,23,24,25]                  # _fe
   orbs[5]  = [26,27,28,29,30]                  # _fe
   orbs[6]  = [31,32,33,34,35,36,37,38,39,40,41,42,43] # _s 
   orbs[7]  = [44,45,46,47,48]                  # _fe
   orbs[8]  = [49,50,51,52,53]                  # _fe
   orbs[9]  = [54,55,56,57,58]                  # _fe
   orbs[10] = [59,60,61,62,63,64,65,66,67]	# _s  
   orbs[11] = [68,69,70,71,72]                  # _mo [1,1,1,0,0]
   orbs[12] = [73,74,75]                        # _end
   dlst = [0,2,6,10,12]
   felst = [1,3,4,5,7,8,9]
   #=====================================================
   # 	      Fe(II)up
   feOccList = {0:[1.,1.,1.,0.,1.,0.,1.,0.,1.,0.],\
   # 	      Fe(II)dn
    	     1:[1.,1.,0.,1.,0.,1.,0.,1.,0.,1.],\
   # 	      Fe(III)up
   	     2:[1.,0.,1.,0.,1.,0.,1.,0.,1.,0.],\
   # 	      Fe(III)dn
   	     3:[0.,1.,0.,1.,0.,1.,0.,1.,0.,1.]}
   #=====================================================
   return norb,orbs,dlst,felst,feOccList

def configuration(conf):
   norb,orbs,dlst,felst,feOccList = init()
   pos3Up = []
   pos2Up = []
   pos3Dn = []
   pos2Dn = []
   for idx,item in enumerate(conf):
      if item == (3,1):
	 pos3Up.append(felst[idx])
      elif item == (2,1):
	 pos2Up.append(felst[idx])
      elif item == (3,0):
	 pos3Dn.append(felst[idx])
      elif item == (2,0):
	 pos2Dn.append(felst[idx])
   # Conf
   occ = numpy.zeros(norb*2)
   # Doubly occupied
   for dx in dlst:
      occ[[2*i+ispin for i in orbs[dx] for ispin in [0,1]]] = 1.0
   # Mo [ms=-0.5]
   occ[2*orbs[11][0]+1] = 1.0
   occ[2*orbs[11][1]+1] = 1.0
   occ[2*orbs[11][2]+0] = 1.0
   # Setup occ for irons
   # Fe(II)up
   for ife2up in pos2Up:
      occ[[2*i+ispin for i in orbs[ife2up] for ispin in [0,1]]] = feOccList[0]
   # Fe(II)dn
   for ife2dn in pos2Dn:
       occ[[2*i+ispin for i in orbs[ife2dn] for ispin in [0,1]]] = feOccList[1]
   # Fe(III)up
   for ife3up in pos3Up:
       occ[[2*i+ispin for i in orbs[ife3up] for ispin in [0,1]]] = feOccList[2]
   # Fe(III)dn
   for ife3dn in pos3Dn:
       occ[[2*i+ispin for i in orbs[ife3dn] for ispin in [0,1]]] = feOccList[3]
   # Save
   conf = numpy.array(occ)
   ncls = sum([2*len(orbs[i]) for i in dlst])
   na = numpy.sum(conf[0::2])
   nb = numpy.sum(conf[1::2]) 
   print ' (na,nb)=',(na,nb),' q[metal],2m=',na+nb-ncls,na-nb,\
         ' (3Up,2Up,3Dn,2Dn)=',pos3Up,pos2Up,pos3Dn,pos2Dn
   return conf

from zmpo_dmrg.source.mpsmpo import mps_class
from zmpo_dmrg.source import mpo_dmrg_io
import h5py
def dumpMPS(conf,fname):
   k = len(conf)
   n = int(sum(conf))
   na = numpy.sum(conf[0::2])
   nb = numpy.sum(conf[1::2])
   spin = (na-nb)/2.0
   # MPS
   mps0 = mps_class.class_mps(k)
   mps0.hfstate(n,conf)
   mps0 = mps0.merge([[2*i,2*i+1] for i in range(k/2)])
   mps0.prt()
   with h5py.File(fname,'w') as f:
      mpo_dmrg_io.dumpMPS(f,mps0.sites,icase=0)
      f['qnum%s'%(k/2)] = numpy.array([[n,spin]])
   return 0

# import os
# import shutil
# prefix = './tmp'
# for idx in range(len(configurations)):
# 
#    conf = configurations[idx]
#    dname = prefix+str(idx)
#    try: 
#       shutil.rmtree(dname)
#    except:  
#       pass
#    os.mkdir(dname)
#    shutil.copy('./mole.h5',dname+'/')
# 
#    # File
#    f0 = open('./template.py')
#    f1 = open(dname+'/main.py','w')
# 
#    p0 = 'mol.tmpdir'
#    l0 = 'mol.tmpdir = "/scratch/global/zhendong/femo_bs7a_'+str(idx)+'/"'
#    p1 = 'conf = '
#    l1 = 'conf = '+str(conf)+'\n'
#    for line in f0.readlines():
#    
#       if p0 in line:
#          f1.writelines(l0)
#       elif p1 in line:	
#          f1.writelines(l1)
#       else:
#          f1.writelines(line)
# 
#    f0.close()
#    f1.close()
#    
#    # File
#    f0 = open('./script.cmd')
#    scriptname = 'tmp_'+str(idx)+'.cmd'
#    f1 = open(dname+'/'+scriptname,'w')
# 
#    p0 = 'export SCRATCHDIR'
#    l0 = 'export SCRATCHDIR="/scratch/global/zhendong/femo_bs7a_'+str(idx)+'"\n'
#    for line in f0.readlines():
#    
#       if p0 in line:
#          f1.writelines(l0)
#       else:
#          f1.writelines(line)
# 
#    f0.close()
#    f1.close()
# 
#    # Submit jobs
#    print 'idx=',idx
#    os.chdir(dname)
#    os.system('sbatch '+'tmp_'+str(idx)+'.cmd')
#    os.chdir('..')
# 
# print 'finished.'
