# -*- coding: utf-8 -*-
from scipy import *
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d  as p3
#import random



def density(star_number):
  """Returns how the star density varries with repect to the raduis"""
  
  return range(star_number)

#def z_location():
  ##loc=random.randomint()
  #return random.normal(loc=0,scale=20)

#def speed(r):
  #"""  """
  ##r_scale=5
  #r_opt=10
  #B=50
  #v_c=7.2
  #if r<(.7*r_opt):
    #x=r/r_opt
    #v_sq=B*(1.97*x**1.22)/((x**2+0.78**2)**1.43)
    #return v_sq**(.5)
  #else:
    #return .2*r+v_c

def star_speed(r,fit):
  """  """
  i=0
  temp=0
  while(r>temp):
    temp=fit[i][0]
    i+=1
  i+=-1
  #v=(fit[i+1][1]-fit[i][1])/(fit[i+1][0]-fit[i][0])*(r-fit[i][0])+fit[i][1]
  
  return (fit[i+1][1]-fit[i][1])/(fit[i+1][0]-fit[i][0])*(r-fit[i][0])+fit[i][1]

x=[]
y=[]
z=[]
xv=[]
yv=[]
zv=[]
rt=[]
sp=[]
mass=[]

xt=[]
yt=[]
valt=[]
#r_t=arange(0,50000,10)
#sp=[]
#Read in Speed Fit File

text_file = open("MWVR_inferred.out", "r")
fit=[]
for line in text_file.readlines():
  if not line[0]=='#':
    #fit.append([1000*float(s) for s in line.rstrip().split()])
    temp=[float(s) for s in line.rstrip().split()]
    temp[0]=temp[0]*1000
    fit.append(temp)

#for i in r_t:
  #sp.append(speed(i,fit))
#plt.figure()
#plt.plot(r_t,sp)
#text_file.close()



###Create Arms  
#for i in range(n):
  #for j in range(m):
    #t=j+(pi/8)*random.random()
    #X_t=(0.2*(i**2)+0.001)*sin(t)
    #Y_t=(0.3*(i**2)+0.001)*cos(t)
    
    #a=(50*i)*pi/(2.0*n)
    
    #x.append(X_t*cos(a)-Y_t*sin(a)+100*random.random()-50)
    #y.append(Y_t*cos(a)+X_t*sin(a)+100*random.random()-50)
    #z.append(z_location())

n=50.0
m=100.0

for i in range(int(n)):
  #print str(i)+"     "+str(int(m*exp(-2*i/n)))+"      "+str(m*exp(-2*i/n))
  for j in range(int(m*exp(-1.0*i/n))):
    #t=j#+(pi/16)*random.random()
    b=(0.2*(i**2)+0.0001)
    c=(0.4*(i**2)+0.0001)
    
    X_t= lambda t:b*sin(t)
    Y_t= lambda t:c*cos(t)
    
    a=(3*i)*pi/(2.0*n)
    
    x_r=lambda t:X_t(t)*cos(a)-Y_t(t)*sin(a)
    y_r=lambda t:Y_t(t)*cos(a)+X_t(t)*sin(a)
    
    x.append(32*(x_r(j)*(1+0.2*(random.random()-1))))
    y.append(32*(y_r(j)*(1+0.2*(random.random()-1))))
    
    z.append(random.normal(loc=0,scale=150))
    rt.append((x[-1]**2+y[-1]**2)**0.5)
    speed=star_speed(rt[-1],fit)
    speed+=(speed*0.05*random.random()-speed*0.05)
    sp.append(speed)
    
    #find unit vectors for function
    #dt=j*0.002
    dxdt=b*cos(j)*cos(a)+c*sin(j)*sin(a)
    dydt=-c*sin(j)*cos(a)+b*cos(j)*sin(a)
    
    val=(dxdt**2+dydt**2)**-0.5
    xu=val*dxdt
    yu=val*dydt

    #print dydx
    xt.append(xu)
    yt.append(yu)
    
    xv.append(xu*speed)
    yv.append(yu*speed)
    valt.append((xv[-1]**2+yv[-1]**2)**0.5)
    zv.append(1)
    
    
    #rt.append((x[-1]**2+y[-1]**2+z[-1]**2)**0.5)
 
#plt.figure()
#plt.hist(xv,100)
#plt.title('dx')
#plt.figure()
#plt.hist(yv,100)
#plt.title('dy')
#plt.figure()
#plt.plot(rt,valt,'.k')
#plt.title('Val')

##Create Bulge;  30 Million Stars Compared to 200 Billion in Disk 
n_c=len(x)*0.015
#print n_c
for i in range(int(n_c)):
  r=4200*random.power(3)
  theta=random.random()*pi
  thi=random.random()*2*pi
  
  rt.append(r)
  x.append(r*cos(thi)*sin(theta))
  y.append(r*sin(thi)*sin(theta))
  z.append(r*cos(theta))
  
  #Velocity
  speed=star_speed(r,fit)
  speed+=(speed*0.05*random.random()-speed*0.05)
  sp.append(speed)

  rdot=speed*(0.01)*(1-random.random())
  a=2*pi*random.random()
  thetadot=cos(a)*speed
  thidot=sin(a)*speed
  
  xv.append(cos(thi)*sin(theta)*rdot-r*sin(thi)*sin(theta)*thidot+r*cos(thi)*cos(theta)*thetadot)
  yv.append(sin(thi)*sin(theta)*rdot+r*cos(thi)*sin(theta)*thidot+r*sin(thi)*cos(theta)*thetadot)
  zv.append(cos(theta)*rdot-r*sin(theta)*thetadot)
  

if len(x)%2==1:
  x.pop()
  y.pop()
  z.pop()
  xv.pop()
  yv.pop()
  zv.pop()
out=open('out.txt','w')
out.write('%012d'%len(x)+'\n')
for i in range(len(x)):
  out.write('% 018.9f'%x[i]+',''% 018.9f'%y[i]+','+'% 018.9f'%z[i]+','+'% 018.9f'%xv[i]+','+'% 018.9f'%yv[i]+','+'% 018.9f'%zv[i]+',''% 018.9f'%1+'\n')

out.close()
#plt.figure()
#plt.plot(rt,sp,'.k')



maxval=max(x)*0.8
maxval_r=(maxval/0.8)*1.5
maxval2=(maxval_r)**2
xdark=[]
ydark=[]
zdark=[]
#Create Dark Matter Halo
for i in arange(-maxval_r,maxval_r,maxval_r/(2**4)):
  for j in arange(-maxval_r,maxval_r,maxval_r/(2**4)):
    for k in arange(-maxval_r,maxval_r,maxval_r/(2**4)):
      if (i**2+j**2+k**2)<maxval2:
        xdark.append(i)
        ydark.append(j)
        zdark.append(k)

out=open('out_dark.txt','w')

for i in range(len(x)):
  out.write('% 018.9f'%xdark[i]+',''% 018.9f'%ydark[i]+','+'% 018.9f'%zdark[i]+'\n')

out.close()




#print len(xdark)
#print len(ydark)
#print len(zdark)
#maxval=max(x)*0.8
fig=plt.figure() 
ax = p3.Axes3D(fig)  
ax.plot3D(x,y,z,'.y')
ax.plot3D(xdark,ydark,zdark,'.k')
ax.set_xlim3d(-maxval_r,maxval_r)
ax.set_ylim3d(-maxval_r,maxval_r)
ax.set_zlim3d(-maxval_r,maxval_r)
print len(x)
plt.show()