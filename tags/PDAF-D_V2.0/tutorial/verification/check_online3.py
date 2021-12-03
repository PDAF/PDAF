#!/use/bin/env python2

import numpy as np
import argparse as ap

if __name__ == "__main__":
   parser = ap.ArgumentParser()
   parser.add_argument('refdir')
   parser.add_argument('refdir2')
   args = parser.parse_args()

   limit=1.e-9

   diffs = []
   maxcheck = 0.0
   imax = 0

   for i in range (0,10):
      if i==0:
         fname = 'stateB_step18_ana.txt'
      else:
         fname = 'ensB_0'+str(i)+'_step18_ana.txt'
      field = np.loadtxt(args.refdir+'/'+fname)
      field_ref = np.loadtxt('verification/'+args.refdir2+'/'+fname) 		

      f_ref = field_ref.reshape(-1)
      f = field.reshape(-1)

      diff1 = f_ref-f

      mdiff = max(abs(diff1).reshape(-1))
      diffs.append(mdiff)
      if i==0:
         imax = i
         maxcheck = mdiff
      else:
         if mdiff>maxcheck:
            imax = i
            maxcheck = mdiff


   adiffs=np.array(diffs)
   if (max(abs(adiffs.reshape(-1)))<limit):
      print ('\033[92mCheck %10.2e   %s at member %3i  OK\033[0m'% (max(abs(adiffs).reshape(-1)), args.refdir, imax))
   else:
      print ('\033[91mCheck ---> %10.2e   %s at member %3i  WARNING!\033[0m'% (max(abs(adiffs).reshape(-1)), args.refdir, imax))
