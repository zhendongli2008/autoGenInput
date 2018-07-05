import geninit
import geninput

fe2 = 2
fe3 = 3
up = 1
dw = 0

#S=3:
seed = [(fe3,dw)]*2 + [(fe2,up)]*5+[(fe2,dw)]*1
#S=4:
#seed = [(fe3,up)] + [(fe3,dw)] + [(fe2,up)]*4 + [(fe2,dw)]*2

# Simply add a screening function is OK to fix position of fe3 !!!
def ifsave(item):
   # third position
   if item[2][0] != fe3 or item[3][0] != fe3: 
      return False
   else:
      return True

dic,keys,flst = geninit.gen(seed,ifsave)
dic,keys,flst = geninit.updown(flst)

confs_eval = []
idx = 0
for key in keys:
   confs = dic[key]
   confs_eval.append(confs[0])
   print idx,confs
   idx += 1
   occ = geninput.configuration(confs[0])
exit()

print '\nGenOCCs:'
for idx,conf in enumerate(confs_eval):
   print 'idx=',idx,conf
   #Input occ for ZMPO code
   occ = geninput.configuration(conf)
   geninput.dumpMPS(occ)
