# -*- coding: utf-8 -*-
from scipy import *
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d  as p3
#import random


x=[]
y=[]
z=[]
xv=[]
yv=[]
zv=[]
rt=[]

mass=[]

xt=[]
yt=[]
valt=[]
#r_t=arange(0,50000,10)
#sp=[]



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




#Read in Speed Fit File
text_file = open("MWVR_inferred.out", "r")
fit=[]
for line in text_file.readlines():
  if not line[0]=='#':
    #fit.append([1000*float(s) for s in line.rstrip().split()])
    temp=[float(s) for s in line.rstrip().split()]
    temp[0]=temp[0]*1000
    temp[1]=temp[1]*0.00102268944 
    #temp[
    fit.append(temp)
    #print temp
text_file.close()

#Adjust these numbers to change start density
#Basically n is the number of ellipse to be created and m is the density of stars along those arms. 
#Warning: If m gets to large with respect to n you must adjust the spin rate and increase n to keep it physically real
#Warning: If n gets to large the galaxy will become to large, you must adjust scalling factor bellow
#for about 12,000 stars use n=75 m=250                                                                                                              
n=75.0
m=50.0
##Create Arms


for i in range(int(n)):
  #print str(i)+"     "+str(int(m*exp(-2*i/n)))+"      "+str(m*exp(-2*i/n))
  for j in range(int(m*exp(-1.0*i/n))):
    #Adjust number that multiplies i/n above to change scaling 
    #t=j#+(pi/16)*random.random()
    
    #Spin rates
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
    rt.append((x[-1]**2+y[-1]**2+z[-1]**2)**0.5)
   
    speed=star_speed(rt[-1],fit)
    speed+=(speed*0.05*random.random()-speed*0.05)
    
    
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
    zv.append(random.normal(scale=0.017)+0.000017*z[-1])
    valt.append((xv[-1]**2+yv[-1]**2+zv[-1]**2)**0.5)
    #rt.append((x[-1]**2+y[-1]**2+z[-1]**2)**0.5)

 

##Create Bulge;  30 Million Stars Compared to 200 Billion in Disk 
n_c=len(x)*0.015
#n_c=1000
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


  rdot=speed*(0.01)*(1-random.random())
  a=2*pi*random.random()
  thetadot=cos(a)*speed/r
  thidot=sin(a)*speed/r
  
  xv.append(cos(thi)*sin(theta)*rdot-r*sin(thi)*sin(theta)*thidot+r*cos(thi)*cos(theta)*thetadot)
  yv.append(sin(thi)*sin(theta)*rdot+r*cos(thi)*sin(theta)*thidot+r*sin(thi)*cos(theta)*thetadot)
  zv.append(cos(theta)*rdot-r*sin(theta)*thetadot)
  valt.append((xv[-1]**2+yv[-1]**2+zv[-1]**2)**0.5)


#if len(x)%2==1:
  #x.pop()
  #y.pop()
  #z.pop()
  #xv.pop()
  #yv.pop()
  #zv.pop()
  
  
#Find mass of stars approx

totalmass=6*10**11
MassPerStar=totalmass/len(x)
#print MassPerStar

out=open('out.txt','w')
out.write('%012d'%len(x)+'\n')
for i in range(len(x)):
  out.write('% 018.9f'%x[i]+',''% 018.9f'%y[i]+','+'% 018.9f'%z[i]+','+'% 018.9f'%xv[i]+','+'% 018.9f'%yv[i]+','+'% 018.9f'%zv[i]+',''% 020.7f'%MassPerStar+'\n')

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
for i in arange(-maxval_r,maxval_r,maxval_r/(2**3.5)):
  for j in arange(-maxval_r,maxval_r,maxval_r/(2**3.5)):
    for k in arange(-maxval_r,maxval_r,maxval_r/(2**3.5)):
      if (i**2+j**2+k**2)<maxval2:
        xdark.append(i)
        ydark.append(j)
        zdark.append(k)

out=open('out_dark.txt','w')
DarkMassPer=(totalmass*5)/len(xdark)
#print DarkMassPer
out.write('%012d'%len(xdark)+'\n')
for i in range(len(xdark)):
  out.write('% 018.9f'%xdark[i]+',''% 018.9f'%ydark[i]+','+'% 018.9f'%zdark[i]+','+'% 020.7f'%DarkMassPer+'\n')

out.close()

#print max(rt)
#maxval_r=maxval

#print "number of dark matter pieces:"+str(len(xdark))

fig=plt.figure() 
ax = p3.Axes3D(fig)  
ax.plot3D(x,y,z,'.y')
ax.plot3D(xdark,ydark,zdark,'.k')
ax.set_xlim3d(-maxval_r,maxval_r)
ax.set_ylim3d(-maxval_r,maxval_r)
ax.set_zlim3d(-maxval_r,maxval_r)
print len(x)
plt.show()