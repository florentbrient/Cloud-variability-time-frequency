import numpy as np

# Calculate confidence interval of slope (scatterplot)
# Florent Brient
# florent.brient@gmail.com


def intver(F0,slope,int,stda,stdb) :

  # Number of points to plot
  Nm = 1000
  pp = np.zeros((5,len(F0)))
  pp[0,:] = slope*F0+int
  pp[1,:] = (slope+stdx)*F0+int
  pp[2,:] = (slope-stdx)*F0+int
  pp[3,:] = (slope)*F0+int+stdb
  pp[4,:] = (slope)*F0+int-stdb

  pmin=np.min(pp,axis=0)
  pmax=np.max(pp,axis=0)

  return pmin,pmax




  



