# -*- coding: utf-8 -*-
 
from scipy import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d  as p3

import subprocess                 # For issuing commands to the OS.
import os
import sys  
from glob import glob

x=[]
y=[]
z=[]
fps=10
path=str(sys.argv[1])
i=0
#print sorted(glob(path+'/outFile*.txt'))
for name in sorted(glob(path+'outFile*.txt')):
  i+=1
  f=open(name,'r')
  print "Doing "+name+" now"
  #Read data in
  x=[]
  y=[]
  z=[]
  for line in f:
    data=line.split(',')
    x.append(float(data[0]))
    y.append(float(data[1]))
    z.append(float(data[2]))
  f.close()
  #Set Max Val on first run through
  if i==1:
    maxval=max(x)
  
  #Plot Data
  fig=plt.figure() 
  ax = p3.Axes3D(fig)  
  ax.plot3D(x,y,z,'.y')
  #ax.plot3D(xdark,ydark,zdark,'.k')
  ax.set_xlim3d(-maxval,maxval)
  ax.set_ylim3d(-maxval,maxval)
  ax.set_zlim3d(-maxval,maxval)
  filename='M'+str('%05d' % i) + '.png'
  #print len(x)
  plt.savefig(filename, dpi=100)
  plt.close()








#for i in range(0,len(Xe),len(Xe)/(250)):
  #plt.figure(2)
  #plt.plot(Xe[i],Ye[i],'bo',label='Earth')
  #plt.plot(Xj[i],Yj[i],'ro',label='Jupiter')
  #plt.ylabel('Y (gigameters)')
  #plt.xlabel('X (gigameters)')
  #plt.title('Earth and Jupiter Orbits, Coupled') 
  #plt.legend(numpoints=1)
  #plt.axis('equal')
  #filename='M'+str('%03d' % i) + '.png'
  #plt.savefig(filename, dpi=100)
  #plt.close()


command = ('mencoder',
           'mf://M*.png',
           '-mf',
           'type=png:w=800:h=600:fps='+str(fps),
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4',
           '-oac',
           'copy',
           '-o',
           'output.avi')

os.spawnvp(os.P_WAIT, 'mencoder', command)

print "\n\nabout to execute:\n%s\n\n" % ' '.join(command)
subprocess.check_call(command)

print "\n\n The movie was written to 'output.avi'"

print "\n\n Deleting M*.png now.\n\n"
os.system('rm M*.png')
