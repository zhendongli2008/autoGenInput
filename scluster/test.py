import geninit
import geninput

fe2 = 2
fe3 = 3
up = 1
dw = 0

seed = [(fe3,dw)]*2 + [(fe2,up)]*3+[(fe2,dw)]*3

# Simply add a screening function is OK to fix position of fe3 !!!
def ifsave(item):
   # third position
   if item[0][0] != fe3 or item[7][0] != fe3: 
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

print '\nGenOCCs:'
for idx,conf in enumerate(confs_eval):
   print 'idx=',idx,conf
   #Input occ for ZMPO code
   occ = geninput.configuration(conf)
   geninput.dumpMPS(occ)
exit()
