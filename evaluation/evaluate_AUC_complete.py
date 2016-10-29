import subprocess
import sys
import numpy as np
from cindex_measure import cindex

temp = np.loadtxt(sys.argv[1])

Y_ = temp[:,0]
Y = temp[:,1]
ans = cindex(Y,Y_)
print('AUC',ans,'\n')
