import numpy
from nitrogenase.spinAnalysis import spinTools
import h5py

k = 76
n = 113
s = 1.5
# Classifier
orbs = [None]*13
orbs[0]  = [0,1]                             # _end
orbs[1]  = [2,3,4,5,6]                       # _fe
orbs[2]  = [7,8,9,10,11,12,13,14,15]         # _s 
orbs[3]  = [16,17,18,19,20]                  # _fe
orbs[4]  = [21,22,23,24,25]                  # _fe
orbs[5]  = [26,27,28,29,30]                  # _fe
orbs[6]  = [31,32,33,34,35,36,37,38,39,40,41,42,43] # _s 
orbs[7]  = [44,45,46,47,48]                  # _fe
orbs[8]  = [49,50,51,52,53]                  # _fe
orbs[9]  = [54,55,56,57,58]                  # _fe
orbs[10] = [59,60,61,62,63,64,65,66,67]      # _s  
orbs[11] = [68,69,70,71,72]                  # _mo [1,1,1,0,0]
orbs[12] = [73,74,75]                        # _end
groups = [orbs[1],orbs[3],orbs[4],orbs[5],\
	  orbs[7],orbs[8],orbs[9],orbs[11]]
#groups = [orbs[i] for i in range(13)]
ng = len(groups)

def genCorr(ftwopdm='twopdm0.npy',fname='ss.h5'):
   with h5py.File(fname,'w') as f: 

      rdm2 = numpy.load(ftwopdm)
      print 'idx=',idx,rdm2.shape
      
      # G[ijkl]=1/2*<ais1+ajs2^+aks2als1>
      twopdm = 2.0*rdm2
      onepdm = numpy.einsum('ijjl->il',rdm2)/(n-1)*2.0
      
      nelec = numpy.einsum('ii',onepdm)
      print 'nelec/n=',nelec,n
      print 'n2-n=',numpy.einsum('ijji',twopdm),n**2-n
      print '<s2>=',spinTools.local_spinsquare(onepdm,twopdm)
      assert abs(nelec-n)<1.e-3
      
      si2 = numpy.zeros(ng)
      for i,ig in enumerate(groups):
         rdm1 = onepdm[numpy.ix_(ig,ig)]
         rdm2 = twopdm[numpy.ix_(ig,ig,ig,ig)]
         sij = spinTools.local_spinsquare(rdm1,rdm2)
         #print 'i=',i,'ig=',ig,sij
         si2[i] = sij
      print     
      for i in range(len(groups)):
         print 'igroup=',i,'s2exp=',si2[i],'seff=',spinTools.from_spinsquare_to_spin(si2[i])
 
      sisj = numpy.zeros((ng,ng))
      for i,ig in enumerate(groups):
         for j,jg in enumerate(groups):
            if i>=j: continue 
            bas1 = ig+jg
            rdm1 = onepdm[numpy.ix_(bas1,bas1)].copy()
            rdm2 = twopdm[numpy.ix_(bas1,bas1,bas1,bas1)].copy()
            sij = spinTools.sisj(rdm1,rdm2,len(ig),len(jg))
            sisj[i,j] = sij
            sisj[j,i] = sij
         sisj[i,i] = si2[i] 
      print
      print 's2sum=',numpy.sum(sisj)
      print sisj
      f['sisj'] = sisj
      
      spdm = spinTools.szHS(onepdm,twopdm,n,s)
      print
      print 'tr(Sz)=',numpy.trace(spdm)
      
      sz = numpy.zeros(ng)
      for i,ig in enumerate(groups):
         rdm1 = spdm[numpy.ix_(ig,ig)]
         print 'i,ig=',i,ig,' rdm1diag=',numpy.diag(rdm1)
         sz[i] = numpy.trace(rdm1)
      print 'sz=',sz,numpy.sum(sz)
      print
      f['sz'] = sz
      
      ne = numpy.zeros(ng)
      for i,ig in enumerate(groups):
         rdm1 = onepdm[numpy.ix_(ig,ig)]
         print 'i,ig=',i,ig,' rdm1diag=',numpy.diag(rdm1)
         ne[i] = numpy.trace(rdm1)
      print 'N=',ne,numpy.sum(ne)
      f['ne'] = ne 
      
      print
      print '='*80
      print

# 
# # S = \sum_i -si^2*log(si^2)
# def vonNeumannEntropy(sigs,thresh=1.e-12):
#    ssum = 0.
#    for sig2 in sigs:
#      if sig2 < thresh: continue
#      ssum += -sig2*numpy.log(sig2)
#    return ssum
# 
#   sisj2 = sisj - numpy.einsum('i,j->ij',sz,sz)
#   print 'sisj-si*sj'
#   print sisj2
#
#   # Sp - entanglement entropy
#   for ig in range(4):
#      print '\nig=',ig
#      for i in groups[ig]:
#         rdm1s = spdm[i,i] # 0.5*(nu-nd)
#         rdm1 = onepdm[i,i] # nu+nd
#         rdm2 = twopdm[i,i,i,i] # 2*nund due to the summation over spins
#         ud = 0.5*rdm2
#         u = 0.5*rdm1 + rdm1s 
#         d = 0.5*rdm1 - rdm1s 
#         vv = 1-u-d+ud
#         uv = u-ud
#         vd = d-ud
# 	 assert abs(vv+uv+vd+ud-1)<1.e-8
#	 svon = vonNeumannEntropy([vv,uv,vd,ud])	
#         style = 'i=%3d'+' '+3*'%8.2f'+' '+6*'%8.2f'
#         print style%(i,rdm1,rdm1s,rdm2,vv,uv,vd,ud,svon,numpy.log(4))

def loadRDMs(fname,ftwopdm):
   with open(fname,'r') as f:
      line = f.readline()
      n = int(line)
      print 'no. of orbs = ',n
      rdm2 = numpy.zeros((n,n,n,n))
      for i in range(n):
       #print 'i=',i
       for j in range(n):
        for k in range(n):
         for l in range(n):
            line = f.readline()
            data = line.split()
            rdm2[i,j,k,l] = float(data[-1])
   numpy.save(ftwopdm,rdm2)
   print 'finished loading'
   return 0


if __name__ == '__main__':
   import time
   for idx in range(35):
      print 'idx=',idx
      ftwopdm = './rdms/twopdm'+str(idx)+'.npy'
      fname   = './rdms/ss'+str(idx)+'.h5'
      ftxt    = './tmp'+str(idx)+'/node_new/spatial_twopdm.0.0.txt'
      t0 = time.time()
      #loadRDMs(ftxt,ftwopdm)
      genCorr(ftwopdm,fname)
      t1 = time.time()
      print 'dt=',t1-t0
