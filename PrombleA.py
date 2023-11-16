import numpy as np
import math
n=5
eig1=np.array([2*(math.cos(i*math.pi/(n+1))-1) for i in range(1,n+1)])
eig2=np.array([2*(math.cos(i*math.pi/(n))-1) for i in range(1,n)])
eig1=eig1[::-1]
eig2=eig2[::-1]

def Compute_w(eig1,eig2):
  temp=np.zeros(n)
  for i in range(0,n-1):
    temp[i]=eig1[i+1]
  w=np.ones(n)
  for i in range(0,n):
    for j in range(0,n-1):
      w[i]=w[i]*(eig1[i]-eig2[j])/(eig1[i]-temp[j])
    temp[i]=eig1[i]
  return w

def Compute_a_b_p(eig1,eig2):
  a=np.zeros(n)
  b=np.zeros(n-1)
  

  pkml=np.zeros(n)
  pk=np.zeros(n)
  
  
  w=Compute_w(eig1,eig2)
  s=np.sum(w)
  a[0]=np.sum(eig1*w)/s
  
  for i in range(0,n):
    pkml[i]=1
    pk[i]=eig1[i]-a[0]
    
  
  
  for k in range(1,n):
    s_tmp=s
    s=0
    t=0
    for i in range(0,n):
      p=w[i]*pk[i]*pk[i]
      s+=p
      t+=eig1[i]*p
    b[k-1]=s/s_tmp
    a[k]=t/s
    for i in range(0,n):
      p=pk[i]
      pk[i]=(eig1[i]-a[k])*p-b[k-1]*pkml[i]
      pkml[i]=p
  for i in range(0,n-1):
    b[i]=math.sqrt(b[i])
  return a,b

a,b=Compute_a_b_p(eig1,eig2)