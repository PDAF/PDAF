#!/use/bin/env python2

import numpy as np
import argparse as ap

if __name__ == "__main__":
   parser = ap.ArgumentParser()
   parser.add_argument('filename')
   parser.add_argument('refdir')
   args = parser.parse_args()

   fname = args.filename
   field = np.loadtxt(fname)
   field_ref = np.loadtxt(args.refdir+fname) 		

   #   diff = field_ref[len(field)-1,1]-field[len(field)-1,1]
   diff1 = field_ref[:,1]-field[:,1]
   diff2 = field_ref[:,2]-field[:,2]

   limit=1.e-9

   if (max(abs(diff1))<limit and max(abs(diff2))<limit):
      print max(abs(diff1)), max(abs(diff2)), fname
   else:
      print '---> ', max(abs(diff1)), max(abs(diff2)), fname, 'WARNING!'
