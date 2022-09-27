import cv2
import numpy as np
import flag

def encrypt(img,key):

  img="tutu.png"
  im=cv2.imread(img)[:,:,(2,1,0)]
  # im=cv2.imread(img)
  print(im)
  # print(cv2.imread(img)[:,:,(2,1,0)])
  [w,h,dim]=im.shape
  a0=key[0]
  p0=key[1]
  u0=key[2]
  v0=key[3]
  w0=key[4]
  x0=key[5]
  y0=key[6]
  z0=key[7]
  q0=key[8]
  pixels = im.flatten(order = 'C')#将数组变为一维
  ai=[]
  for i in range(3*w*h):
    if 0<=a0<p0:
      a0=a0/p0
    elif a0<0.5:
      a0=(a0-p0)*(0.5-p0)
    else:
      a0=1-a0
    ai.append(a0)
  dic=list(zip(ai,pixels))
  dic.sort(key=lambda x:x[0])#x[0]升序排列
  pixels=list(list(zip(*dic))[1])
  R=pixels[:w*h]
  G=pixels[w*h:2*w*h]
  B=pixels[2*w*h:]
  print(R)
  print(B)

#这几个参数是IDS
  t=100
  f=10
  r=28
  g=8/3
  Y,Z,Q=Lorenz(y0,z0,q0,f,r,g,t+w*h)#PCS
  Y=Y[t:]
  Z=Z[t:]
  Q=Q[t:]
  Y_R=list(zip(Y,R))
  Y_R.sort(key=lambda x:x[0])
  R=list(list(zip(*Y_R))[1])
  Z_G=list(zip(Z,G))
  Z_G.sort(key=lambda x:x[0])
  G=list(list(zip(*Z_G))[1])
  Q_B=list(zip(Q,B))
  Q_B.sort(key=lambda x:x[0])
  B=list(list(zip(*Q_B))[1])

  a=36
  b=3
  c=28
  d=16
  k=0.2
  t=100
  U,V,W,X=Chen(u0,v0,w0,x0,a,b,c,d,k,t+3*w*h)#KDS

  U=U[t:]
  V=V[t:]
  W=W[t:]
  X=X[t:]
  #a1
  for i in range(3*w*h):
    rule='ACGT'
    if(int(U[i]%1/0.05) in [0,4,8,10,19]):
      rule='AGCT'
    elif(int(U[i]%1/0.05) in [1,6,12,14,17]):
      rule='ACGT'
    elif(int(U[i]%1/0.05) in [2,7,11,13,16]):
      rule='GATC'
    elif(int(U[i]%1/0.05) in [3,5,9,15,18]):
      rule='CATG'
    if(i/(w*h)<1):
      R[i]=DNA_Encode(R[i],rule)
    elif(i/(w*h)<2):
      G[i-w*h]=DNA_Encode(G[i-w*h],rule)
    else:
      B[i-2*w*h]=DNA_Encode(B[i-2*w*h],rule)#DS

  start=[]
  times=[]
  for i in V:
    start.append(int(i*pow(10,12))%8)
  for i in W:
    times.append(int(i*pow(10,12))%8)

  startR=start[:w*h]
  startG=start[w*h:2*w*h]
  startB=start[2*w*h:]
  timesR=times[:w*h]
  timesG=times[w*h:2*w*h]
  timesB=times[2*w*h:]
  rules=['ACGT','CATG','GTAC','TCGA','CTAG','AGCT','TGCA','GATC']#KDS
  for i in range(w*h):
    s=startR[i]
    for j in range(timesR[i]):
      R[i]=DNA_XOR(R[i],rules[s])
      s=(s+1)%8

  for i in range(w*h):
    s=startG[i]
    for j in range(timesG[i]):
      G[i]=DNA_XOR(G[i],rules[s])
      s=(s+1)%8

  for i in range(w*h):
    s=startB[i]
    for j in range(timesB[i]):
      B[i]=DNA_XOR(B[i],rules[s])
      s=(s+1)%8

#这里的RGB已经变成了DS


  for i in range(3*w*h):
    rule='ACGT'
    if(int(X[i]%1/0.05) in [0,4,8,10,19]):
      rule='GTAC'
    elif(int(X[i]%1/0.05) in [1,6,12,14,17]):
      rule='TGCA'
    elif(int(X[i]%1/0.05) in [2,7,11,13,16]):
      rule='CTAG'
    elif(int(X[i]%1/0.05) in [3,5,9,15,18]):
      rule='TCGA'
    if(i/(w*h)<1):
      R[i]=DNA_Decode(R[i],rule)
    elif(i/(w*h)<2):
      G[i-w*h]=DNA_Decode(G[i-w*h],rule)
    else:
      B[i-2*w*h]=DNA_Decode(B[i-2*w*h],rule)
#step10
  encrypt_img=np.array((R+G+B)).reshape((512,512,3),order='C')
  return encrypt_img


def Lorenz(x0,y0,z0,p,q,r,T):
  h=0.01
  x=[]
  y=[]
  z=[]
  for t in range(T):
    xt=x0+h*p*(y0-x0)
    yt=y0+h*(q*x0-y0-x0*z0)
    zt=z0+h*(x0*y0-r*z0)
    x0,y0,z0=xt,yt,zt
    x.append(x0)
    y.append(y0)
    z.append(z0)
  return x,y,z

def Chen(u0,v0,w0,x0,a,b,c,d,k,T):
  h=0.001
  u=[]
  v=[]
  w=[]
  x=[]
  for t in range(T):
    ut=u0+h*(a*(v0-u0))
    vt=v0+h*(-u0*w0+d*u0+c*u0-x0)
    wt=w0+h*(u0*v0-b*w0)
    xt=u0+k
    #u0、v0、w0,x0统一更新
    u0,v0,w0,x0=ut,vt,wt,xt
    u.append(u0)
    v.append(v0)
    w.append(w0)
    x.append(x0)
  return u,v,w,x


def DNA_Encode(pixel,rule):
  base=''
  bits=bin(pixel)[2:].zfill(8)
  for k in range(4):
    b=bits[k*2:2*k+2]
    if b=='00':
      base+=rule[0]
    elif b=='01':
      base+=rule[1]
    elif b=='10':
      base+=rule[2]
    else:
      base+=rule[3]
  return base

def DNA_Decode(base,rule):
  pixel=''
  for k in base:
    if k==rule[0]:
      pixel+='00'
    elif k==rule[1]:
      pixel+='01'
    elif k==rule[2]:
      pixel+='10'
    else:
      pixel+='11'
  return int(pixel,2)

def DNA_XOR(base1,base2):
  pixel=DNA_Decode(base1,'AGCT')^DNA_Decode(base2,'AGCT')
  return DNA_Encode(pixel,'AGCT')

def main():
  key = [0.49226688, 0.28059747, 0.87321577, 0.63073925, 0.66753483, 0.49983341, 0.37095885, 0.12800098, 0.14163127, 0.23561871]
  img_encrypt=encrypt(flag,key)
  cv2.imwrite('./123.tiff',img_encrypt)


def reversedna(CIDS,KDS):




if __name__ == '__main__':
  main()