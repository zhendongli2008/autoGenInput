def single(charge,spin):
   if charge == 2 and spin == 1:
      q = +2.0
      m = +2.0
   elif charge == 2 and spin == 0:
      q = +2.0
      m = -2.0
   elif charge == 3 and spin == 1:
      q = +3.0
      m = +2.5
   elif charge == 3 and spin == 0:
      q = +3.0
      m = -2.5
   return q,m

def qnum(qs):
  q = 0.
  m = 0.
  for iq in qs:
     q1,m1 = single(iq[0],iq[1])
     q += q1
     m += m1
  return q,m
