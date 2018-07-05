import numpy
import itertools

def init():
   norb = 71
   # Classifier
   orbs = [0]*13
   orbs[0]  = [0,1]                                # _end
   orbs[1]  = [2,3,4,5,6]                          # _fe 
   orbs[2]  = [7,8,9,10,11,12,13,14,15,16]         # _s  
   orbs[3]  = [17,18,19,20,21]                     # _fe 
   orbs[4]  = [22,24,25,26,27]                     # _fe 
   orbs[5]  = [23,28,29,30,31]                     # _fe 
   orbs[6]  = [32,33,34,35,36,37,38]               # _s  
   orbs[7]  = [39,40,41,42,47]                     # _fe 
   orbs[8]  = [43,44,45,46,48]                     # _fe 
   orbs[9]  = [49,50,51,52,53]                     # _fe 
   orbs[10] = [54,55,56,57,58,59,60,61,62,63]      # _s  
   orbs[11] = [64,65,66,67,68]                     # _fe 
   orbs[12] = [69,70]                              # _end
   # doubly occupancy
   dlst = [0,2,6,10,12]
   print numpy.sum([2.0*len(orbs[i]) for i in dlst]) 
   felst = [1,3,4,5, 7,8,9,11]
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
   #>>># Mo [ms=-0.5]
   #>>>occ[2*orbs[11][0]+1] = 1.0
   #>>>occ[2*orbs[11][1]+1] = 1.0
   #>>>occ[2*orbs[11][2]+0] = 1.0
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

def dumpMPSwithQnums(conf,fname):
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
   #================================== 
   # DUMP Qnums also ...
   #================================== 
   qnums = [numpy.array([[0.0,0.0]])]
   ne = 0.0
   ms = 0.0
   for i in range(k/2):
      ne += conf[2*i]+conf[2*i+1]
      ms += (conf[2*i]-conf[2*i+1])/2.0
      qnum = numpy.array([[ne,ms]])
      qnums.append(qnum) 
   with h5py.File(fname,'w') as f:
      mpo_dmrg_io.dumpMPS(f,[mps0.sites,qnums],icase=1)
   #================================== 
   return 0

if __name__ == '__main__':
   init()
