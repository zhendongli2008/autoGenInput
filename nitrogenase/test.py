import geninit
import geninput

fe2 = 2
fe3 = 3
up = 1
dw = 0
seed = [(fe2,up)]*2+[(fe2,dw)]+[(fe3,up)]*2+[(fe3,dw)]*2

dic,keys,flst = geninit.gen(seed)
print keys
#print flst
exit()

dic,keys,flst = geninit.updown(flst)
print len(flst)
print len(keys)
print keys
exit()

confs_eval = []
for key in keys:
   confs = dic[key]
   confs_eval.append(confs[0])

for idx,conf in enumerate(confs_eval):
   print 'idx=',idx,conf
   occ = geninput.configuration(conf)
   geninput.dumpMPS(occ)
   exit()
