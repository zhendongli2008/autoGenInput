import geninit
import geninput

fe2 = 2
fe3 = 3
up = 1
dw = 0
seed = [(fe2,up)]*4+[(fe2,dw)]*4 

dic,keys,flst = geninit.gen(seed,itr=1)
dic,keys,flst = geninit.updown(flst)
print dic

confs_eval = []
for key in keys:
   confs = dic[key]
   confs_eval.append(confs[0])

print '\nGenOCCs:'
for idx,conf in enumerate(confs_eval):
   print 'idx=',idx,conf
   #Input occ for ZMPO code
   occ = geninput.configuration(conf)
   geninput.dumpMPS(occ)
