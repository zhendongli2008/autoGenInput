import itertools
import qsym

def classification(dic):
   flst = []
   keys = sorted(dic.keys())
   kdx = 0
   counter = 0
   for idx,key in enumerate(keys):
      ln = len(dic[key])
      counter += ln
      print '*idx=',idx,' key=',key,' len=',ln,' counter=',counter
   print
   for idx,key in enumerate(keys):
      ln = len(dic[key])
      counter += ln
      print '*idx=',idx,' key=',key,' len=',ln,' counter=',counter
      flst += dic[key] 
      for jdx,conf in enumerate(dic[key]):
         print ' jdx/kdx=',(jdx,kdx),' conf=',conf
         kdx += 1
   print 'total=',len(flst)
   return keys,flst

def gen(seed,ifsave):
   # Generate all of them by permutations
   perms = list(itertools.permutations(range(8),8)) 
   # Brute-force generation
   unique = []
   for iperm in perms:
      new = [seed[i] for i in iperm]
      if new not in unique and ifsave(new) == 1:
         unique.append(new)
   # Classification
   dic = {}
   for conf in unique:
      ql = qsym.qnum(conf[:4]) # left four FEs
      qr = qsym.qnum(conf[4:]) # right three FEs + MO
      key = str([ql,qr])
      if key not in dic:
         dic[key] = [conf]
      else:
         dic[key].append(conf)
   # sort 
   keys,flst = classification(dic)
   return dic,keys,flst

def updown(flst):
   dic = {}
   for conf in flst:
      key = str(map(lambda x:x[1],conf))
      if key not in dic:
         dic[key] = [conf]
      else:
         dic[key].append(conf)
   # sort 
   keys,flst = classification(dic)
   return dic,keys,flst
