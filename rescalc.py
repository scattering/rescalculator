from __future__ import division
import sys
import numpy as N
import matplotlib
matplotlib.use('Agg')
import pylab
import math
import unittest
from matplotlib.patches import Ellipse
from matplotlib.ticker import NullFormatter, MultipleLocator
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import MaxNLocator
import copy,sys

import lattice_calculator

eps=1e-3
pi=N.pi

def sign(x):
     if x>0:
          ret=1
     if x<0:
          ret=-1
     if x==0:
          ret=0
     return ret

def blkdiag(g):
     "Returns a block diagonal matrix, given a list of square matrices, g"
     glen=len(g)
     n=0
     for i in range(glen):
     #     print g[i]
          n=n+g[i].shape[0]
     #     print n


     gout=N.zeros((n,n))
     offset=0
     for i in range(glen):
          currblock=g[i]
          lenx,leny=currblock.shape
     #    print i
     #    print lenx, leny
          for x in range(lenx):
               for y in range(leny):
                    gout[x+offset,y+offset]=currblock[x,y]
     #           print gout
          offset=offset+lenx
     return gout


def similarity_transform(A,B):
     G=N.dot(B,A.transpose())
     G2=N.dot(A,G)
     return G2

class rescalculator:
     def __init__(self, mylattice):
          self.lattice_calculator=mylattice
          return

     def ResMat(self, Q, W,EXP):
          CONVERT1=0.4246609*N.pi/60/180;
          CONVERT2=2.072;
          npts=len(EXP)
          Qnpts=N.size(Q)
          RM=N.zeros((4, 4, npts),'d');
          R0=N.zeros((npts,1),'d');
          RM_=N.zeros((4, 4),'d');
          D=N.zeros((8, 13),'d');
          d=N.zeros((4, 7),'d');
          T=N.zeros((4, 13),'d');
          t=N.zeros((2, 7),'d');
          A=N.zeros((6, 8),'d');
          C=N.zeros((4, 8),'d');
          B=N.zeros((4, 6),'d');

          for ind in range(npts):
#            %Assign default values and decode parameters
               moncor=1;
               if 'moncor' in EXP[ind]:
                    moncor = EXP[ind]['moncor']
               alpha = EXP[ind]['hcol']*CONVERT1
               beta =  EXP[ind]['vcol']*CONVERT1
               mono=EXP[ind]['mono']
               etam = mono['mosaic']*CONVERT1
               etamv=etam
               if 'vmosaic' in mono:
                    etamv = mono['vmosaic']*CONVERT1
               ana=EXP[ind]['ana']
               etaa = ana['mosaic']*CONVERT1
               etaav=etaa
               if 'vmosaic' in ana:
                    etaav = ana['vmosaic']*CONVERT1
               sample=EXP[ind]['sample']
               infin=-1
               if 'infin' in EXP[ind]:
                    infin = EXP[ind]['infin']
               efixed=EXP[ind]['efixed']
               epm=1
               if 'dir1' in EXP[ind]:
                    epm= EXP[ind]['dir1']
               ep=1
               if 'dir2' in EXP[ind]:
                    ep= EXP[ind]['dir2']
               monitorw=1.0
               monitorh=1.0
               beamw=1.0
               beamh=1.0
               monow=1.0
               monoh=1.0
               monod=1.0
               anaw=1.0
               anah=1.0
               anad=1.0
               detectorw=1.0
               detectorh=1.0
               sshape=N.eye(3)
               L0=1.0
               L1=1.0
               L1mon=1.0
               L2=1.0
               L3=1.0
               monorv=1e6
               monorh=1e6
               anarv=1e6
               anarh=1e6
               if 'beam' in EXP[ind]:
                    beam=EXP[ind]['beam']
                    if 'width' in beam:
                         beamw=beam['width']**2
                    if 'height' in beam:
                         beamh=beam['height']**2
               bshape=N.diag([beamw,beamh])
               if 'monitor' in EXP[ind]:
                    monitor=EXP[ind]['monitor']
                    if 'width' in monitor:
                         monitorw=monitor['width']**2
                    monitorh=monitorw
                    if 'height' in monitor:
                         monitorh=monitor['height']**2;
               monitorshape=N.diag([monitorw,monitorh])
               if 'detector' in EXP[ind]:
                    detector=EXP[ind]['detector']
                    if 'width' in detector:
                         detectorw=detector['width']**2
                    if 'height' in detector:
                         detectorh=detector['height']**2
               dshape=N.diag([detectorw,detectorh])
               if 'width' in mono:
                    monow=mono['width']**2
               if 'height' in mono:
                    monoh=mono['height']**2
               if 'depth' in mono:
                    monod=mono.depth**2
               mshape=N.diag([monod,monow,monoh])
               if 'width' in ana:
                    anaw=ana['width']**2
               if 'height' in ana:
                    anah=ana['height']**2
               if 'depth' in ana:
                    anad=ana['depth']**2
               ashape=N.diag([anad,anaw,anah])
               if 'shape' in sample:
                    sshape=sample['shape']
               if 'arms' in EXP[ind]:
                    arms=EXP[ind]['arms']
                    L0=arms[0]
                    L1=arms[2]
                    L2=arms[3]
                    L3=arms[4]
                    L1mon=L1
                    if len(arms)>3:
                         L1mon=arms[4]
               if 'rv' in mono:
                    monorv=mono['rv']
               if 'rh' in mono:
                    monorh=mono['rh']
               if 'rv' in ana:
                    anarv=ana['rv']
               if 'rh' in ana:
                    anarh=ana['rh']
               method=0;
               if 'method' in EXP[ind]:
                    method=EXP[ind]['method']
               myinstrument=lattice_calculator.Instrument()
               taum=myinstrument.get_tau(mono['tau'])
               taua=myinstrument.get_tau(ana['tau'])

               horifoc=-1;
               if 'horifoc' in EXP[ind]:
                    horifoc=EXP[ind]['horifoc']
               if horifoc==1:
                    alpha[2]=alpha[2]*N.sqrt(8*N.log(2)/12)
#            %---------------------------------------------------------------------------------------------
#            %Calculate angles and energies
               w=W[ind]
               q=Q[ind]
               ei=efixed
               ef=efixed
               if infin>0:
                    ef=efixed-w
               else:
                    ei=efixed+w
               ki = N.sqrt(ei/CONVERT2)
               kf = N.sqrt(ef/CONVERT2)
               thetam=N.arcsin(taum/(2*ki))*sign(epm)
               thetaa=N.arcsin(taua/(2*kf))*sign(ep)
               s2theta=-N.arccos( (ki**2+kf**2-q**2)/(2*ki*kf))# %2theta sample
               thetas=s2theta/2
               phi=N.arctan2(-kf*N.sin(s2theta), ki-kf*N.cos(s2theta)) #%Angle from ki to Q

     #           %---------------------------------------------------------------------------------------------
#            %Calculate beam divergences defined by neutron guides
               pi=N.pi
               if alpha[0]<0:
                    alpha[0]=-alpha[0]*2*0.427/ki*pi/180
               if alpha[1]<0:
                    alpha[1]=-alpha[1]*2*0.427/ki*pi/180
               if alpha[2]<0:
                    alpha[2]=-alpha[2]*2*0.427/ki*pi/180
               if alpha[3]<0:
                    alpha[3]=-alpha[3]*2*0.427/ki*pi/180

               if beta[0]<0:
                    beta[0]=-beta[0]*2*0.427/ki*pi/180
               if beta[1]<0:
                    beta[1]=-beta[1]*2*0.427/ki*pi/180
               if beta[2]<0:
                    beta[2]=-beta[2]*2*0.427/ki*pi/180
               if beta[3]<0:
                    beta[3]=-beta[3]*2*0.427/ki*pi/180

#            %---------------------------------------------------------------------------------------------
#            %Rededine sample geometry
               psi=thetas-phi;# %Angle from sample geometry X axis to Q
               rot=N.zeros((3,3))
               rot[0,0]=N.cos(psi)
               rot[1,1]=N.cos(psi)
               rot[0,1]=N.sin(psi)
               rot[1,0]=-N.sin(psi)
               rot[2,2]=1
               sshape=N.dot(rot.transpose(),N.dot(sshape,rot)) #matrix multiplication?
#            %---------------------------------------------------------------------------------------------
#            %Definition of matrix G
               G=1.0/N.array([alpha[0],alpha[1],beta[0],beta[1],alpha[2],alpha[3],beta[2],beta[3]])**2
               G=N.diag(G)
#            %---------------------------------------------------------------------------------------------
#            %Definition of matrix F
               F=1.0/N.array([etam,etamv,etaa,etaav])**2
               F=N.diag(F);
#            %---------------------------------------------------------------------------------------------
#            %Definition of matrix A
               A[0,0]=ki/2/N.tan(thetam)
               A[0,1]=-A[0,0]
               A[3,4]=kf/2/N.tan(thetaa)
               A[3,5]=-A[3,4]
               A[1,1]=ki
               A[2,3]=ki
               A[4,4]=kf
               A[5,6]=kf
#            %---------------------------------------------------------------------------------------------
#            %Definition of matrix C
               C[0,0]=1.0/2
               C[0,1]=1.0/2
               C[2,4]=1.0/2
               C[2,5]=1.0/2
               C[1,2]=1.0/(2*N.sin(thetam))
               C[1,3]=-C[1,2];# %mistake in paper
               C[3,6]=1.0/(2*N.sin(thetaa))
               C[3,7]=-C[3,6]
#            %---------------------------------------------------------------------------------------------
#            %Definition of matrix B
               B[0,0]=N.cos(phi)
               B[0,1]=N.sin(phi)
               B[0,3]=-N.cos(phi-s2theta)
               B[0,4]=-N.sin(phi-s2theta)
               B[1,0]=-B[0,1]
               B[1,1]=B[0,0]
               B[1,3]=-B[0,4]
               B[1,4]=B[0,3]
               B[2,2]=1.0
               B[2,5]=-1.0
               B[3,0]=2*CONVERT2*ki
               B[3,3]=-2*CONVERT2*kf
#            %---------------------------------------------------------------------------------------------
#            %Definition of matrix S
               Sinv=blkdiag([bshape,mshape,sshape,ashape,dshape])# %S-1 matrix
               S=N.linalg.inv(Sinv)
#            %---------------------------------------------------------------------------------------------
#            %Definition of matrix T
               T[0,0]=-1./(2*L0); # %mistake in paper
               T[0,2]=N.cos(thetam)*(1./L1-1./L0)/2
               T[0,3]=N.sin(thetam)*(1./L0+1./L1-2./(monorh*N.sin(thetam)))/2
               T[0,5]=N.sin(thetas)/(2*L1)
               T[0,6]=N.cos(thetas)/(2*L1)
               T[1,1]=-1./(2*L0*N.sin(thetam))
               T[1,4]=(1./L0+1./L1-2*N.sin(thetam)/monorv)/(2*N.sin(thetam))
               T[1,7]=-1./(2*L1*N.sin(thetam))
               T[2,5]=N.sin(thetas)/(2*L2)
               T[2,6]=-N.cos(thetas)/(2*L2)
               T[2,8]=N.cos(thetaa)*(1./L3-1./L2)/2
               T[2,9]=N.sin(thetaa)*(1./L2+1./L3-2/(anarh*N.sin(thetaa)))/2
               T[2,11]=1./(2*L3)
               T[3,7]=-1./(2*L2*N.sin(thetaa))
               T[3,10]=(1./L2+1./L3-2*N.sin(thetaa)/anarv)/(2*N.sin(thetaa))
               T[3,12]=-1./(2*L3*N.sin(thetaa))
#            %---------------------------------------------------------------------------------------------
#            %Definition of matrix D
#            % Lots of index mistakes in paper for matix D
               D[0,0]=-1./L0
               D[0,2]=-N.cos(thetam)/L0
               D[0,3]=N.sin(thetam)/L0
               D[2,1]=D[0,0]
               D[2,4]=-D[0,0]
               D[1,2]=N.cos(thetam)/L1
               D[1,3]=N.sin(thetam)/L1
               D[1,5]=N.sin(thetas)/L1
               D[1,6]=N.cos(thetas)/L1
               D[3,4]=-1./L1
               D[3,7]=-D[3,4]
               D[4,5]=N.sin(thetas)/L2
               D[4,6]=-N.cos(thetas)/L2
               D[4,8]=-N.cos(thetaa)/L2
               D[4,9]=N.sin(thetaa)/L2
               D[6,7]=-1./L2
               D[6,10]=-D[6,7]
               D[5,8]=N.cos(thetaa)/L3
               D[5,9]=N.sin(thetaa)/L3
               D[5,11]=1./L3
               D[7,10]=-D[5,11]
               D[7,12]=D[5,11]
#            %---------------------------------------------------------------------------------------------
#            %Definition of resolution matrix M
               if method==1:
                    Minv=N.dot(N.dot(B,A),N.dot(N.linalg.inv(D.linalg.inv(D*N.linalg.inv(S+T.transpose()*F*T)*D.transpose())+G),\
                                                N.dot(A.transpose(),B.transpose())))#; %Popovici
                    #Minv=B*A*((D*(S+T'*F*T)^(-1)*D')^(-1)+G)^(-1)*A'*B'; %Popovici
               else:
                    #%Horizontally focusing analyzer if needed
                    #print 'intermediate C.T*F*C'
                    #print N.dot(C.transpose(),N.dot(F,C))
                    #print 'product'
                    #print N.dot(F,C)
                    HF_int=N.linalg.inv(G+N.dot(C.transpose(),N.dot(F,C)))
                    HF=similarity_transform(A,HF_int)
                    #print 'HF'
                    #print HF
                    if horifoc>0:
                         HF=N.linalg.inv(HF);
                         HF[4,4]=(1.0/(kf*alpha[2]))**2
                         HF[4,3]=0
                         HF[3,4]=0
                         HF[3,3]=(N.tan(thetaa)/(etaa*kf))**2
                         HF=N.linalg.inv(HF)
                    Minv=similarity_transform(B,HF)#; %Cooper-Nathans
               M=N.linalg.inv(Minv)
               #print 'A'
               #print A
               #print 'B'
               #print B
               #print 'C'
               #print C
               #print 'D'
               #print D
               #print 'T'
               #print T
               #print 'M'
               #print M
               #print 'Minv'
               #print Minv
               #print 'G'
               #print G
               #print 'F'
               #print F
               RM_[0,0]=M[0,0]
               RM_[1,0]=M[1,0]
               RM_[0,1]=M[0,1]
               RM_[1,1]=M[1,1]

               RM_[0,2]=M[0,3]
               RM_[2,0]=M[3,0]
               RM_[2,2]=M[3,3]
               RM_[2,1]=M[3,1]
               RM_[1,2]=M[1,3]

               RM_[0,3]=M[0,2]
               RM_[3,0]=M[2,0]
               RM_[3,3]=M[2,2]
               RM_[3,1]=M[2,1]
               RM_[1,3]=M[1,2]
     #           %--------------------------------------------------------------------------------------------
     #           %Calculation of prefactor, normalized to source
               #print 'RM_'
               #print RM_
               Rm=ki**3/N.tan(thetam)
               Ra=kf**3/N.tan(thetaa)
               #print 'Rm'
               #print Rm
               #print 'Ra'
               #print Ra
               if method==1:
                    R0_=Rm*Ra*(2*pi)**4/(64*pi**2*N.sin(thetam)*N.sin(thetaa))\
                       *N.sqrt(N.linalg.det(F)/N.linalg.det\
                               (N.linalg.inv(N.dot(D,N.dot(N.linalg.inv(S+N.dot(T.T,N.dot(F,T))),D.T)))+G)) #%Popovici
               else:
                    R0_=Rm*Ra*(2*pi)**4/(64*pi**2*N.sin(thetam)*N.sin(thetaa))\
                       *N.sqrt( N.linalg.det(F)/N.linalg.det(G+N.dot(C.transpose(),N.dot(F,C)))) #%Cooper-Nathans
                    #print 'RO_'
                    #print R0_
#            %---------------------------------------------------------------------------------------------
#            %Normalization to flux on monitor
               if moncor==1:
                    #print 'Moncor'
                    #print R0_
                    #print 'RM in Moncor'
                    #print RM_
                    g=G[0:4][:,0:4]
                    f=F[0:2][:,0:2]
                    c=C[0:2][:,0:4]
                    #print 'f'
                    #print f
                    #print 'g'
                    #print g
                    #print 'c'
                    #print c

                    t[0,0]=-1./(2*L0)#  %mistake in paper
                    t[0,2]=N.cos(thetam)*(1./L1mon-1./L0)/2
                    t[0,3]=N.sin(thetam)*(1./L0+1./L1mon-2./(monorh*N.sin(thetam)))/2
                    t[0,6]=1./(2*L1mon)
                    t[1,1]=-1./(2*L0*N.sin(thetam))
                    t[1,4]=(1./L0+1./L1mon-2*N.sin(thetam)/monorv)/(2*N.sin(thetam))
                    sinv=blkdiag([bshape,mshape,monitorshape])# %S-1 matrix
                    s=N.linalg.inv(sinv)
                    d[0,0]=-1./L0
                    d[0,2]=-N.cos(thetam)/L0
                    d[0,3]=N.sin(thetam)/L0
                    d[2,1]=D[0,0]
                    d[2,4]=-D[0,0]
                    d[1,2]=N.cos(thetam)/L1mon
                    d[1,3]=N.sin(thetam)/L1mon
                    d[1,5]=0
                    d[1,6]=1./L1mon
                    d[3,4]=-1./L1mon
                    if method==1:
                         Rmon=Rm*(2*pi)**2/(8*pi*N.sin(thetam))*N.sqrt(N.linalg.det(f)/\
                                                                       N.linalg.det(N.linalg.inv(N.dot(d,N.dot(N.linalg.inv(s+t.transpose()*f*t),d.transpose())))+g)) #%Popovici
                    else:
                         Rmon=Rm*(2*pi)**2/(8*pi*N.sin(thetam))*N.sqrt(N.linalg.det(f)/N.linalg.det(g+N.dot(c.transpose(),N.dot(f,c)))) #%Cooper-Nathans
                         #print 'Cooper Nathans'
                         #print 'Rmon', Rmon
                    R0_=R0_/Rmon
                    R0_=R0_*ki #%1/ki monitor efficiency
                    #print 'R01', R0_
                    #print 'Rmon', Rmon

#            %---------------------------------------------------------------------------------------------
#            %Transform prefactor to Chesser-Axe normalization
               R0_=R0_/(2*pi)**2*N.sqrt(N.linalg.det(RM_))
#            %---------------------------------------------------------------------------------------------
#            %Include kf/ki part of cross section
               R0_=R0_*kf/ki
#            %---------------------------------------------------------------------------------------------
#            %Take care of sample mosaic if needed [S. A. Werner & R. Pynn, J. Appl. Phys. 42, 4736, (1971)]
               #print 'R0_ before mosaic', R0_
               #print RM_
               if 'mosaic' in sample:
                    etas = sample['mosaic']*CONVERT1;
                    #print 'sample',sample
                    if 'vmosaic' in sample:
                         #print 'vertical'
                         etasv = sample['vmosaic']*CONVERT1
                         #print 'etasv',etasv
                    else:
                         etasv=etas

                    #print 'etas',etas, 'etasv',etasv
                    R0_=R0_/N.sqrt((1.+(q*etas)**2*RM_[3,3])*(1.0+(q*etasv)**2*RM_[1,1]))
                    Minv=N.linalg.inv(RM_)
                    Minv[1,1]=Minv[1,1]+q**2*etas**2
                    Minv[3,3]=Minv[3,3]+q**2*etasv**2
                    RM_=N.linalg.inv(Minv)
                    #print 'mosaic'
                    #print 'R0_', R0_
                    #print 'RM_'
                    #print RM_
#            %---------------------------------------------------------------------------------------------
#            %Take care of analyzer reflectivity if needed [I. Zaliznyak, BNL]
               if ('thickness' in ana) & ('Q' in ana):
                    KQ = ana['Q']
                    KT = ana['thickness']
                    toa=(taua/2)/N.sqrt(kf**2-(taua/2)**2)
                    smallest=alpha[3]
                    if alpha[3]>alpha[2]:
                         smallest=alpha[2]
                    Qdsint=KQ*toa
                    dth=(N.arange(201)/200)*N.sqrt(2*N.log(2))*smallest
                    wdth=N.exp(-dth**2/2./etaa**2)
                    sdth=KT*Qdsint*wdth/etaa/N.sqrt(2.*pi)
                    rdth=1./(1+1./sdth)
                    reflec=rdth.sum()/wdth.sum()
                    R0_=R0_*reflec
#            %---------------------------------------------------------------------------------------------
               #print 'ind ', ind
               #print 'shape ', R0.shape
               R0[ind]=R0_
               RM[:,:,ind]=RM_[:,:]
          return R0, RM

     def ResMatS(self,H,K,L,W,EXP):
# [len,H,K,L,W,EXP]=CleanArgs(H,K,L,W,EXP);
          x=self.lattice_calculator.x
          y=self.lattice_calculator.y
          z=self.lattice_calculator.z
          Q=self.lattice_calculator.modvec(H,K,L,'latticestar')
          uq=N.zeros((3,self.lattice_calculator.npts),'d')
          uq[0,:]=H/Q;  #% Unit vector along Q
          uq[1,:]=K/Q;
          uq[2,:]=L/Q;
          xq=self.lattice_calculator.scalar(x[0,:],x[1,:],x[2,:],uq[0,:],uq[1,:],uq[2,:],'latticestar');
          yq=self.lattice_calculator.scalar(y[0,:],y[1,:],y[2,:],uq[0,:],uq[1,:],uq[2,:],'latticestar');
          zq=0; # %scattering vector assumed to be in (self.orient1,self.orient2) plane;
          tmat=N.zeros((4,4,self.lattice_calculator.npts)); #%Coordinate transformation matrix
          tmat[3,3,:]=1;
          tmat[2,2,:]=1;
          tmat[0,0,:]=xq;
          tmat[0,1,:]=yq;
          tmat[1,1,:]=xq;
          tmat[1,0,:]=-yq;

          RMS=N.zeros((4,4,self.lattice_calculator.npts));
          rot=N.zeros((3,3));
          EXProt=EXP;

#        %Sample shape matrix in coordinate system defined by scattering vector
          for i in range(len(EXP)):
               sample=EXP[i]['sample'];
               if 'shape' in sample:
                    rot[0,0]=tmat[0,0,i];
                    rot[1,0]=tmat[1,0,i];
                    rot[0,1]=tmat[0,1,i];
                    rot[1,1]=tmat[1,1,i];
                    rot[2,2]=tmat[2,2,i];
                    EXProt[i]['sample']['shape']=N.dot(rot,N.dot(sample['shape'],rot.T));

          R0,RM= self.ResMat(Q,W,EXProt)
          #print 'RM ',RM.shape
          #print 'npts ',self.lattice_calculator.npts
          #print 'tmat ',tmat.shape
          for i in range(self.lattice_calculator.npts):
               RMS[:,:,i]=N.dot((tmat[:,:,i]).transpose(),N.dot(RM[:,:,i],tmat[:,:,i]))

          mul=N.zeros((4,4));
          e=N.eye(4,4);
          for i in range(len(EXP)):
               if 'Smooth' in EXP[i]:
                    if 'X' in (EXP[i]['Smooth']):
                         mul[0,0]=1./(EXP[i]['Smooth']['X']**2/8/N.log(2));
                         mul[1,1]=1./(EXP[i]['Smooth']['Y']**2/8/N.log(2));
                         mul[2,2]=1./(EXP[i]['Smooth']['E']**2/8/N.log(2));
                         mul[3,3]=1./(EXP[i]['Smooth']['Z']**2/8/N.log(2));
                         R0[i]=R0[i]/N.sqrt(N.linalg.det(e/RMS[:,:,i]))*N.sqrt(N.linalg.det(e/mul+e/RMS[:,:,i]));
                         RMS[:,:,i]=e/(e/mul+e/RMS[:,:,i]);
          return R0.T, RMS

     def ResPlot(self,H,K,L,W,EXP):
          """Plot resolution ellipse for a given scan"""
          center=N.round(H.shape[0]/2)
          if center<1:
               center=0
          if center>H.shape[0]:
               center=H.shape[0]
          #EXP=[EXP[center]]
          Style1=''
          Style2='--'

          XYAxesPosition=[0.1, 0.6, 0.3, 0.3]
          XEAxesPosition=[0.1, 0.1, 0.3, 0.3]
          YEAxesPosition=[0.6, 0.6, 0.3, 0.3]
          TextAxesPosition=[0.45, 0.0, 0.5, 0.5]
          GridPoints=101

          [R0,RMS]=self.ResMatS(H,K,L,W,EXP)
          #[xvec,yvec,zvec,sample,rsample]=self.StandardSystem(EXP);
          self.lattice_calculator.StandardSystem()
          #print 'shape ',self.lattice_calculator.x.shape
          qx=self.lattice_calculator.scalar(self.lattice_calculator.x[0,:],self.lattice_calculator.x[1,:],self.lattice_calculator.x[2,:],H,K,L,'latticestar')
          qy=self.lattice_calculator.scalar(self.lattice_calculator.y[0,:],self.lattice_calculator.y[1,:],self.lattice_calculator.y[2,:],H,K,L,'latticestar')
          qw=W;

          #========================================================================================================
          #find reciprocal-space directions of X and Y axes

          o1=self.lattice_calculator.orientation.orient1.T#[:,0] #EXP['orient1']
          o2=self.lattice_calculator.orientation.orient2.T#[:,0] #EXP['orient2']
          pr=self.lattice_calculator.scalar(o2[0,:],o2[1,:],o2[2,:],self.lattice_calculator.y[0,:],self.lattice_calculator.y[1,:],self.lattice_calculator.y[2,:],'latticestar')
          o2[0]=self.lattice_calculator.y[0,:]*pr
          o2[1]=self.lattice_calculator.y[1,:]*pr
          o2[2]=self.lattice_calculator.y[2,:]*pr

          if N.abs(o2[0,center])<1e-5:
               o2[0,center]=0.0
          if N.absolute(o2[1,center])<1e-5:
               o2[1,center]=0.0
          if N.absolute(o2[2,center])<1e-5:
               o2[2,center]=0.0

          if N.abs(o1[0,center])<1e-5:
               o1[0,center]=0.0
          if N.absolute(o1[1,center])<1e-5:
               o1[1,center]=0.0
          if N.absolute(o1[2,center])<1e-5:
               o1[2,center]=0.0

          #%========================================================================================================
          #%determine the plot range
          XWidth=max(self.fproject(RMS,0))
          YWidth=max(self.fproject(RMS,1))
          WWidth=max(self.fproject(RMS,2))
          XMax=(max(qx)+XWidth*1.5)
          XMin=(min(qx)-XWidth*1.5)
          YMax=(max(qy)+YWidth*1.5)
          YMin=(min(qy)-YWidth*1.5)
          WMax=(max(qw)+WWidth*1.5)
          WMin=(min(qw)-WWidth*1.5)
          #print 'qx ',qx
          #print 'qy ',qy
          #print 'XWidth ',XWidth
          #print 'YWidth ',YWidth
          fig=pylab.figure()
          #%========================================================================================================
          #% plot XY projection


          proj,sec=self.project(RMS,2)
          (a,b,c)=N.shape(proj)
          mat=N.copy(proj)
          #print 'proj ', proj.shape
          a1=[];b1=[];theta=[];a1_sec=[];b1_sec=[];theta_sec=[];e=[]; e_sec=[]
          for i in range(c):
               matm=N.matrix(mat[:,:,i])
               w,v=N.linalg.eig(matm)
               vm=N.matrix(v)
               vmt=vm.T
               mat_diag=vmt*matm*vm
               a1.append(1.0/N.sqrt(mat_diag[0,0]))
               b1.append(1.0/N.sqrt(mat_diag[1,1]))
               #thetar=N.arccos(vm[0,0])#*N.sign(vm[0,0])
               thetar=-N.arctan2(vm[0,1],vm[0,0])
               theta.append(math.degrees(thetar))

          mat_sec=N.copy(sec)
          #print 'proj ', proj
          (a,b,c)=N.shape(sec)
          for i in range(c):
               matm_sec=N.matrix(mat_sec[:,:,i])
               w_sec,v_sec=N.linalg.eig(matm_sec)
               vm_sec=N.matrix(v)
               vmt_sec=vm_sec.T
               mat_diag_sec=vmt_sec*matm_sec*vm_sec
               a1_sec.append(1.0/N.sqrt(mat_diag_sec[0,0]))
               b1_sec.append(1.0/N.sqrt(mat_diag_sec[1,1]))
               #thetar_sec=N.arccos(vm_sec[0,0])*N.sign(vm_sec[0,0])
               thetar_sec=-N.arctan2(vm_sec[0,1],vm_sec[0,0])
               theta_sec.append(math.degrees(thetar_sec))
               x0y0=N.array([qx[i],qy[i]])
               e.append(Ellipse(x0y0,width=2*a1[i],height=2*b1[i],angle=theta[i]))
               e_sec.append(Ellipse(x0y0,width=2*a1_sec[i],height=2*b1_sec[i],angle=theta_sec[i]))

          #print 'a1_sec ',a1_sec
          #print 'b1_sec ',b1_sec
          #print 'theta_sec ',theta_sec
          #print 'mat_diag_sec ',mat_diag_sec




          rsample='latticestar'
          oxmax=XMax/self.lattice_calculator.modvec(o1[0],o1[1],o1[2],rsample)
          oxmin=XMin/self.lattice_calculator.modvec(o1[0],o1[1],o1[2],rsample)
          oymax=YMax/self.lattice_calculator.modvec(o2[0],o2[1],o2[2],rsample)
          oymin=YMin/self.lattice_calculator.modvec(o2[0],o2[1],o2[2],rsample)
          #print 'a1 ',a1
          #print 'b1 ',b1
          #print 'theta ',theta
          #print 'mat_diag ',mat_diag
          #x0y0=N.array([1.0,0.0])
          #make right y-axis
          ax2 = fig.add_subplot(2,2,1)
          pylab.subplots_adjust(hspace=0.6,wspace=0.3)
          #ax2.set_xlim(oxmin, oxmax)
          ax2.set_ylim(oymin[center], oymax[center])
          ax2.yaxis.tick_right()
          ax2.yaxis.set_label_position('right')
          ax2.xaxis.set_major_formatter(pylab.NullFormatter())
          ax2.xaxis.set_major_locator(pylab.NullLocator())
          ylabel=r'Q$_y$' +'(units of ['+str(o2[0,center])+' '+str(o2[1,center])+' '+str(o2[2,center])+'])'
          ax2.set_ylabel(ylabel)
          #ax2.set_zorder(3)
          #make top x-axis
          if 1:
               ax3 = fig.add_axes(ax2.get_position(), frameon=False,label='x-y top')
               ax3.xaxis.tick_top()
               ax3.xaxis.set_label_position('top')
               ax3.set_xlim(oxmin[center], oxmax[center])
               ax3.yaxis.set_major_formatter(NullFormatter())
               ax3.yaxis.set_major_locator(pylab.NullLocator())
               xlabel=r'Q$_x$' +'(units of ['+str(o1[0,center])+' '+str(o1[1,center])+' '+str(o1[2,center])+'])'
               ax3.set_xlabel(xlabel)
               #ax3.set_zorder(2)

          #make bottom x-axis, left y-axis
          if 1:
               ax = fig.add_axes(ax2.get_position(), frameon=False,label='x-y')
               ax.yaxis.tick_left()
               ax.yaxis.set_label_position('left')
               ax.xaxis.tick_bottom()
               #ax.xaxis.set_label_position('bottom')
               for i in range(c):
                    ax.add_artist(e[i])
                    e[i].set_clip_box(ax.bbox)
                    e[i].set_alpha(0.5)
                    e[i].set_facecolor('red')
                    ax.add_artist(e_sec[i])
                    e_sec[i].set_clip_box(ax.bbox)
                    e_sec[i].set_alpha(0.7)
                    e_sec[i].set_facecolor('blue')

               ax.set_xlim(XMin, XMax)
               ax.set_ylim(YMin, YMax)
               xlabel=r'Q$_x$ ('+r'$\AA^{-1}$)'
               ax.set_xlabel(xlabel)
               ylabel=r'Q$_y$ ('+r'$\AA^{-1}$)'
               ax.set_ylabel(ylabel)
               #ax.set_zorder(1)

          #%========================================================================================================
          #% plot XE projection


          proj,sec=self.project(RMS,1);
          (a,b,c)=N.shape(proj)
          mat=N.copy(proj)
          #print 'proj ', proj
          a1=[];b1=[];theta=[];a1_sec=[];b1_sec=[];theta_sec=[];e=[]; e_sec=[]
          for i in range(c):
               matm=N.matrix(mat[:,:,i])
               w,v=N.linalg.eig(matm)
               vm=N.matrix(v)
               vmt=vm.T
               mat_diag=vmt*matm*vm
               a1.append(1.0/N.sqrt(mat_diag[0,0]))
               b1.append(1.0/N.sqrt(mat_diag[1,1]))
               #thetar=N.arccos(vm[0,0])*N.sign(vm[0,0])
               thetar=-N.arctan2(vm[0,1],vm[0,0])
               theta.append(math.degrees(thetar))
          mat_sec=N.copy(sec)
          #print 'proj ', proj
          (a,b,c)=N.shape(sec)
          for i in range(c):
               matm_sec=N.matrix(mat_sec[:,:,i])
               w_sec,v_sec=N.linalg.eig(matm_sec)
               vm_sec=N.matrix(v)
               vmt_sec=vm_sec.T
               mat_diag_sec=vmt_sec*matm_sec*vm_sec
               a1_sec.append(1.0/N.sqrt(mat_diag_sec[0,0]))
               b1_sec.append(1.0/N.sqrt(mat_diag_sec[1,1]))
               #thetar_sec=N.arccos(vm_sec[0,0])*N.sign(vm_sec[0,0])
               thetar_sec=-N.arctan2(vm_sec[0,1],vm_sec[0,0])
               theta_sec.append(math.degrees(thetar_sec))
               x0y0=N.array([qx[i],qw[i]])
               e.append(Ellipse(x0y0,width=2*a1[i],height=2*b1[i],angle=theta[i]))
               e_sec.append(Ellipse(x0y0,width=2*a1_sec[i],height=2*b1_sec[i],angle=theta_sec[i]))
          rsample='latticestar'
          oxmax=XMax/self.lattice_calculator.modvec(o1[0],o1[1],o1[2],rsample)
          oxmin=XMin/self.lattice_calculator.modvec(o1[0],o1[1],o1[2],rsample)
          oymax=WMax
          oymin=WMin
          #print 'a1 ',a1
          #print 'b1 ',b1
          #print 'theta ',theta
          #print 'mat_diag ',mat_diag
          #x0y0=N.array([1.0,0.0])
          #make right y-axis
          ax2 = fig.add_subplot(2,2,3)
          #ax2.set_xlim(oxmin, oxmax)
          ax2.set_ylim(oymin, oymax)
          ax2.yaxis.tick_right()
          ax2.yaxis.set_label_position('right')
          ax2.xaxis.set_major_formatter(pylab.NullFormatter())
          ax2.xaxis.set_major_locator(pylab.NullLocator())
          ax2.yaxis.set_major_locator(pylab.NullLocator())
          ylabel=r'E' +'(meV)'
          #ax2.set_ylabel(ylabel)
          #ax2.set_zorder(3)
          #make top x-axis
          if 1:
               ax3 = fig.add_axes(ax2.get_position(), frameon=False,label='x-E top')
               ax3.xaxis.tick_top()
               ax3.xaxis.set_label_position('top')
               ax3.set_xlim(oxmin[center], oxmax[center])
               ax3.yaxis.set_major_formatter(NullFormatter())
               ax3.yaxis.set_major_locator(pylab.NullLocator())
               xlabel=r'Q$_x$' +'(units of ['+str(o1[0,center])+' '+str(o1[1,center])+' '+str(o1[2,center])+'])'
               ax3.set_xlabel(xlabel)
               #ax3.set_zorder(2)

          #make bottom x-axis, left y-axis
          if 1:
               ax = fig.add_axes(ax2.get_position(), frameon=False,label='x-E')
               ax.yaxis.tick_left()
               ax.yaxis.set_label_position('left')
               ax.xaxis.tick_bottom()
               #ax.xaxis.set_label_position('bottom')

               for i in range(c):
                    ax.add_artist(e[i])
                    e[i].set_clip_box(ax.bbox)
                    e[i].set_alpha(0.5)
                    e[i].set_facecolor('red')
                    ax.add_artist(e_sec[i])
                    e_sec[i].set_clip_box(ax.bbox)
                    e_sec[i].set_alpha(0.7)
                    e_sec[i].set_facecolor('blue')
               ax.set_xlim(XMin, XMax)
               ax.set_ylim(WMin, WMax)
               xlabel=r'Q$_x$ ('+r'$\AA^{-1}$)'
               ax.set_xlabel(xlabel)
               ylabel=r'E (meV)'
               ax.set_ylabel(ylabel)
               #ax.yaxis.set_major_formatter(NullFormatter())
               #ax.yaxis.set_major_locator(pylab.NullLocator())

               #ax.set_zorder(1)
          #pylab.show()
          #%========================================================================================================
          #% plot YE projection


          proj,sec=self.project(RMS,0);
          (a,b,c)=N.shape(proj)
          mat=N.copy(proj)
          #print 'proj ', proj
          a1=[];b1=[];theta=[];a1_sec=[];b1_sec=[];theta_sec=[];e=[]; e_sec=[]
          for i in range(c):
               matm=N.matrix(mat[:,:,i])
               w,v=N.linalg.eig(matm)
               vm=N.matrix(v)
               vmt=vm.T
               mat_diag=vmt*matm*vm
               a1.append(1.0/N.sqrt(mat_diag[0,0]))
               b1.append(1.0/N.sqrt(mat_diag[1,1]))
               #thetar=N.arccos(vm[0,0])*N.sign(vm[0,0])
               thetar=-N.arctan2(vm[0,1],vm[0,0])
               theta.append(math.degrees(thetar))
          #print 'proj ', proj
          (a,b,c)=N.shape(sec)
          mat_sec=N.copy(sec)
          for i in range(c):
               matm_sec=N.matrix(mat_sec[:,:,i])
               w_sec,v_sec=N.linalg.eig(matm_sec)
               vm_sec=N.matrix(v)
               vmt_sec=vm_sec.T
               mat_diag_sec=vmt_sec*matm_sec*vm_sec
               a1_sec.append(1.0/N.sqrt(mat_diag_sec[0,0]))
               b1_sec.append(1.0/N.sqrt(mat_diag_sec[1,1]))
               #thetar_sec=N.arccos(vm_sec[0,0])*N.sign(vm_sec[0,0])
               thetar_sec=-N.arctan2(vm_sec[0,1],vm_sec[0,0])
               theta_sec.append(math.degrees(thetar_sec))
               x0y0=N.array([qy[i],qw[i]])
               e.append(Ellipse(x0y0,width=2*a1[i],height=2*b1[i],angle=theta[i]))
               e_sec.append(Ellipse(x0y0,width=2*a1_sec[i],height=2*b1_sec[i],angle=theta_sec[i]))

          print 'qx ',qx
          print 'qw ',qw
          print 'theta_sec ',theta_sec
          rsample='latticestar'
          oxmax=YMax/self.lattice_calculator.modvec(o2[0],o2[1],o2[2],rsample)
          oxmin=YMin/self.lattice_calculator.modvec(o2[0],o2[1],o2[2],rsample)
          oymax=WMax
          oymin=WMin
          #print 'a1 ',a1
          #print 'b1 ',b1
          #print 'theta ',theta
          #print 'mat_diag ',mat_diag
          #x0y0=N.array([1.0,0.0])
          #make right y-axis
          ax2 = fig.add_subplot(2,2,4)
          #ax2.set_xlim(oxmin, oxmax)
          ax2.set_ylim(oymin, oymax)
          ax2.yaxis.tick_right()
          ax2.yaxis.set_label_position('right')
          ax2.xaxis.set_major_formatter(pylab.NullFormatter())
          ax2.xaxis.set_major_locator(pylab.NullLocator())
          ylabel=r'E' +'(meV)'
          ax2.set_ylabel(ylabel)
          #ax2.set_zorder(3)
          #make top x-axis
          if 1:
               ax3 = fig.add_axes(ax2.get_position(), frameon=False,label='y-E top')
               ax3.xaxis.tick_top()
               ax3.xaxis.set_label_position('top')
               ax3.set_xlim(oxmin[center], oxmax[center])
               ax3.yaxis.set_major_formatter(NullFormatter())
               ax3.yaxis.set_major_locator(pylab.NullLocator())
               xlabel=r'Q$_y$' +'(units of ['+str(o2[0,center])+' '+str(o2[1,center])+' '+str(o2[2,center])+'])'
               ax3.set_xlabel(xlabel)
               #ax3.set_zorder(2)

          #make bottom x-axis, left y-axis
          if 1:
               ax = fig.add_axes(ax2.get_position(), frameon=False,label='y-E')
               ax.yaxis.tick_left()
               ax.yaxis.set_label_position('left')
               ax.xaxis.tick_bottom()
               #ax.xaxis.set_label_position('bottom')

               for i in range(c):
                    ax.add_artist(e[i])
                    e[i].set_clip_box(ax.bbox)
                    e[i].set_alpha(0.5)
                    e[i].set_facecolor('red')
                    ax.add_artist(e_sec[i])
                    e_sec[i].set_clip_box(ax.bbox)
                    e_sec[i].set_alpha(0.7)
                    e_sec[i].set_facecolor('blue')

               ax.set_xlim(YMin, YMax)
               ax.set_ylim(WMin, WMax)
               xlabel=r'Q$_y$ ('+r'$\AA^{-1}$)'
               ax.set_xlabel(xlabel)
               #ylabel=r'E (meV)'
               #ax.set_ylabel(ylabel)
               #ax.yaxis.set_major_formatter(NullFormatter())
               ax.yaxis.set_major_locator(pylab.NullLocator())

               #ax.set_zorder(1)


          #pylab.show()
          pylab.savefig("resout.png")


          #self.PlotEllipse(proj,qx,qw,Style1);
          #self.PlotEllipse(sec,qx,qw,Style2);
          #pylab.show()
          return



     def fproject(self,mat_in,i):
          """return hwhm of projection"""
          if (i==0):
               v=2
               j=1
          if (i==1):
               v=0
               j=2
          if (i==2):
               v=0
               j=1
          mat=N.array(mat_in)
          (a,b,c)=N.shape(mat)
          proj=N.zeros((2,2,c),'d')
          proj[0,0,:]=mat[i,i,:]-mat[i,v,:]*mat[i,v,:]/mat[v,v,:]
          proj[0,1,:]=mat[i,j,:]-mat[i,v,:]*mat[j,v,:]/mat[v,v,:]
          proj[1,0,:]=mat[j,i,:]-mat[j,v,:]*mat[i,v,:]/mat[v,v,:]
          proj[1,1,:]=mat[j,j,:]-mat[j,v,:]*mat[j,v,:]/mat[v,v,:]
          hwhm=proj[0,0,:]-proj[0,1,:]*proj[0,1,:]/proj[1,1,:]
          hwhm=N.sqrt(2*N.log(2))/N.sqrt(hwhm)
          return hwhm


     def PlotEllipse(self,mat_in,x0,y0,style):
          """plot ellipse"""
          mat=N.array(mat_in)
          (a,b,c)=N.shape(mat)
               #phi=0:2*pi/3000:2*pi;
          phi=N.arange(0,2*pi,2*pi/3000)
          for i in range(c):
               r=N.sqrt(2*N.log(2)/(mat[0,0,i]*N.cos(phi)*N.cos(phi)+mat[1,1,i]*N.sin(phi)*N.sin(phi)\
                                    +2*mat[0,1,i]*N.cos(phi)*N.sin(phi)))
               x=r*N.cos(phi)+x0[i];
               y=r*N.sin(phi)+y0[i];
               pylab.plot(x,y,style);
          return

     def project(self,mat_in,v):
          """return projection and cross section matrices"""
          if v == 2:
               i=0;j=1
          if v == 0:
               i=1;j=2
          if v == 1:
               i=0;j=2
          mat=N.array(mat_in)
          (a,b,c)=N.shape(mat)
          proj=N.zeros((2,2,c),'d')
          sec=N.zeros((2,2,c),'d')
          proj[0,0,:]=mat[i,i,:]-mat[i,v,:]*mat[i,v,:]/mat[v,v,:]
          proj[0,1,:]=mat[i,j,:]-mat[i,v,:]*mat[j,v,:]/mat[v,v,:]
          proj[1,0,:]=mat[j,i,:]-mat[j,v,:]*mat[i,v,:]/mat[v,v,:]
          proj[1,1,:]=mat[j,j,:]-mat[j,v,:]*mat[j,v,:]/mat[v,v,:]
          sec[0,0,:]=mat[i,i,:]
          sec[0,1,:]=mat[i,j,:]
          sec[1,0,:]=mat[j,i,:]
          sec[1,1,:]=mat[j,j,:]
          return proj,sec

     def cooper_nathans_crystal_transform(self,H,K,L):
          "Provides a transformation matrix from the CooperNathan's Coordinate  system to the CrystalCoordinate system"
          x=self.lattice_calculator.x
          y=self.lattice_calculator.y
          z=self.lattice_calculator.z
          Q=self.lattice_calculator.modvec(H,K,L,'latticestar')
          uq=N.zeros((3,self.lattice_calculator.npts),'Float64')
          uq[0,:]=H/Q;  #% Unit vector along Q
          uq[1,:]=K/Q;
          uq[2,:]=L/Q;
          xq=self.lattice_calculator.scalar(x[0,:],x[1,:],x[2,:],uq[0,:],uq[1,:],uq[2,:],'latticestar');
          yq=self.lattice_calculator.scalar(y[0,:],y[1,:],y[2,:],uq[0,:],uq[1,:],uq[2,:],'latticestar');
          zq=0; # %scattering vector assumed to be in (self.orient1,self.orient2) plane;
          tmat=N.zeros((4,4,self.lattice_calculator.npts)); #%Coordinate transformation matrix
          tmat[3,3,:]=1;
          tmat[2,2,:]=1;
          tmat[0,0,:]=xq;
          tmat[0,1,:]=yq;
          tmat[1,1,:]=xq;
          tmat[1,0,:]=-yq;
          return tmat


     def calc_correction(self,H,K,L,W,EXP,qscan=None):
          "Returns the correction factor such that I_int=F^2*correction"
          Q=self.lattice_calculator.modvec(H,K,L,'latticestar')
          R0,RM=self.ResMat(Q,W,EXP)
          th_correction=[]
          tth_correction=[]
          if qscan != None:
               qx=[]
               qy=[]
               qmod=[]
               for i in range(len(qscan)):
                    qh=N.array([qscan[i][0]],'Float64')
                    qk=N.array([qscan[i][1]],'Float64')
                    ql=N.array([qscan[i][2]],'Float64')
                    tmat=(self.cooper_nathans_crystal_transform(qh,qk,ql))
                    qcooper=N.matrix(tmat[0:3,0:3,i])*N.matrix(qscan[i],'Float64').T
                    qscan_mod=self.lattice_calculator.modvec(qh,qk,ql,'latticestar')
                    #print 'qcooper ',qcooper[0,0]
                    qx.append(qcooper[0,0]/qscan_mod)
                    qy.append(qcooper[1,0]/qscan_mod)

          q_correction=[]
          #print 'qx ',qx[0][0]
          for i in range(RM.shape[2]):
               Myy=RM[1,1,i]
               Mxx=RM[0,0,i]
               Mxy=RM[0,1,i]
               th_correction.append(N.sqrt(N.linalg.det(RM[:,:,i]))/N.sqrt(Myy)/Q[i])

               #print 'RM',RM[:,:,i]
               #print 'det', N.sqrt(N.linalg.det(RM[:,:,i]))
               #print 'sqrt',N.sqrt(Myy)
               #print 'Q', Q[i]
               tth_correction.append(N.sqrt(N.linalg.det(RM[:,:,i]))/N.sqrt(Mxx)/Q[i])
               if qscan!=None:
                    q_correction.append(N.sqrt(N.linalg.det(RM[:,:,i]))/N.sqrt(Mxx*qx[i]**2+2*Mxy*qx[i]*qy[i]+Myy*qy[i]**2)/Q[i])
          corrections={}
          corrections['th_correction']=th_correction
          corrections['tth_correction']=tth_correction
          if qscan!=None:
               corrections['q_correction']=q_correction
          return corrections

     def CalcWidths(self,H,K,L,W,EXP):
          """Calculates the full width of Bragg Peaks"""
          Q=self.lattice_calculator.modvec(H,K,L,'latticestar')
          #print 'Q ',Q
          npts=self.lattice_calculator.npts
          #print 'npts ',npts
          center=N.round(H.shape[0]/2)
          if center<1:
               center=0
          if center>H.shape[0]:
               center=H.shape[0]
          #EXP=[EXP[center]]

          [R0c,RMc]=self.ResMat(Q,W,EXP)
          [R0,RMS]=self.ResMatS(H,K,L,W,EXP)
          self.lattice_calculator.StandardSystem()
          qx=self.lattice_calculator.scalar(self.lattice_calculator.x[0,:],self.lattice_calculator.x[1,:],self.lattice_calculator.x[2,:],H,K,L,'latticestar')
          qy=self.lattice_calculator.scalar(self.lattice_calculator.y[0,:],self.lattice_calculator.y[1,:],self.lattice_calculator.y[2,:],H,K,L,'latticestar')
          qw=W;

          #========================================================================================================
          #find reciprocal-space directions of X and Y axes

          o1=self.lattice_calculator.orientation.orient1.T #EXP['orient1']
          o2=self.lattice_calculator.orientation.orient2.T #EXP['orient2']
#        print 'o2shape ',o2.shape
#        print 'yshape ',self.lattice_calculator.y.shape
          pr=self.lattice_calculator.scalar(o2[0,:],o2[1,:],o2[2,:],self.lattice_calculator.y[0,:],self.lattice_calculator.y[1,:],self.lattice_calculator.y[2,:],'latticestar')
          o2[0]=self.lattice_calculator.y[0,:]*pr
          o2[1]=self.lattice_calculator.y[1,:]*pr
          o2[2]=self.lattice_calculator.y[2,:]*pr


          XWidth=max(self.fproject(RMS,0))
          YWidth=max(self.fproject(RMS,1))
          WWidth=max(self.fproject(RMS,2))

          proj,sec=self.project(RMS,1);
          (a,b,c)=N.shape(proj)
          mat=N.copy(proj)
          for i in range(c):
               matm=N.matrix(mat[:,:,i])
               w,v=N.linalg.eig(matm)
               vm=N.matrix(v)
               vmt=vm.T
               mat_diag=vmt*matm*vm
          a1=1.0/N.sqrt(mat_diag[0,0])
          b1=1.0/N.sqrt(mat_diag[1,1])
          thetar=N.arccos(vm[0,0])*N.sign(vm[0,0])
          theta=math.degrees(thetar)
          XWidth=2*self.fproject(RMS,0)
          YWidth=2*self.fproject(RMS,1)
          WWidth=2*self.fproject(RMS,2)
          ZWidth=2*N.sqrt(2*N.log(2))/N.sqrt(RMS[3,3,:])
          XBWidth=2*N.sqrt(2*N.log(2))/N.sqrt(RMS[0,0,:])
          YBWidth=2*N.sqrt(2*N.log(2))/N.sqrt(RMS[1,1,:])
          WBWidth=2*N.sqrt(2*N.log(2))/N.sqrt(RMS[2,2,:])

          XWidthc=1*self.fproject(RMc,0)
          YWidthc=1*self.fproject(RMc,1)
          ZWidthc=1*N.sqrt(2*N.log(2))/N.sqrt(RMc[3,3,:])
          XBWidthc=1*N.sqrt(2*N.log(2))/N.sqrt(RMc[0,0,:])
          YBWidthc=1*N.sqrt(2*N.log(2))/N.sqrt(RMc[1,1,:])

          CONVERT1=0.4246609*N.pi/60/180;
          CONVERT2=2.072;
          for ind in range(npts):
               infin=-1
               if 'infin' in EXP[ind]:
                    infin = EXP[ind]['infin']
               efixed=EXP[ind]['efixed']
               w=W[ind]
               q=Q[ind]
               ei=efixed
               ef=efixed
               if infin>0:
                    ef=efixed-w
               else:
                    ei=efixed+w;
               ki = N.sqrt(ei/CONVERT2)
               kf = N.sqrt(ef/CONVERT2)
          S2=N.arccos((ki**2+kf**2-Q**2)/(2*ki*kf))
                    #			todeg=1./reselpsq/!dtor
                    #			; Q^2 = ki^2 + kf^2 + 2 ki kf cos(A4)
                    #			; 2 Q dQ = 2 ki kf (-)sin(a4) d A4
                    #			; d a4 = -Q dQ / ( ki kf sin(A4) )
                    #			todega4=(reselpsq/ki/kf/sin(a4*!dtor))/!dtor

          tthwidth_sec=XBWidthc*2*Q/ki/kf/N.sin(S2)
          tthwidth_proj=XWidthc*2*Q/ki/kf/N.sin(S2)
          thwidth_sec=(YBWidthc*2)/Q  # recall for small theta, arc=r*theta
          thwidth_proj=(YWidthc*2)/Q


##        #M1,M2,S1,S2,A1,A2=mylattice.SpecGoTo(H+XWidth,K,L,W,EXP)
##        uq=N.zeros((3,self.lattice_calculator.npts),'Float64')
##        #Ql=N.sqrt(H**2+K**2+L**2)
##        uq[0,:]=H/Q  #% Unit vector along Q
##        uq[1,:]=K/Q
##        uq[2,:]=L/Q
##        #projection
##        Hup=uq[0,:]*XWidthc
##        Kup=uq[1,:]*XWidthc
##        Lup=uq[2,:]*XWidthc
##        M1,M2,S1,S2,A1,A2=mylattice.SpecGoTo(H+Hup,K+Kup,L+Lup,W,EXP)
##        M1c,M2c,S1c,S2c,A1c,A2c=mylattice.SpecGoTo(H,K,L,W,EXP)
##        tthwidth_proj=N.absolute(S2-S2c)*2
##        #cross section
##        Hup=uq[0,:]*XBWidthc
##        Kup=uq[1,:]*XBWidthc
##        Lup=uq[2,:]*XBWidthc
##        print 'Hup ', Hup
##        print 'Hadded ', H+Hup
##        print 'Kup ', Kup
##        print 'Kadded ', K+Kup
##        print 'Lup ', Lup
##        print 'Ladded ', L+Lup
##
##        M1,M2,S1,S2,A1,A2=self.lattice_calculator.SpecGoTo(H+Hup,K+Kup,L+Lup,W,EXP)
##        M1c,M2c,S1c,S2c,A1c,A2c=self.lattice_calculator.SpecGoTo(H,K,L,W,EXP)
##        tthwidth_sec=N.absolute(S2-S2c)*2
##        print 'S2 ',math.degrees(S2)
##        print 'S2c ',math.degrees(S2c)
##        print 'tthwidth_sec ',tthwidth_sec

##        XWidth=fproject (RMS,1);
##        YWidth=fproject (RMS,2);
##        WWidth=fproject (RMS,3);
##        ZWidth=sqrt(2*log(2))./sqrt(RMS(4,4,:));
##
##        XBWidth=sqrt(2*log(2))./sqrt(RMS(1,1,:));
##        YBWidth=sqrt(2*log(2))./sqrt(RMS(2,2,:));
##        WBWidth=sqrt(2*log(2))./sqrt(RMS(3,3,:));
##
##
##        ResVol=(2*pi)^2/sqrt(det(RMS(:,:,center)));


          Widths={}
          Widths['XWidth']=XWidth
          Widths['YWidth']=YWidth
          Widths['ZWidth']=ZWidth
          Widths['WWidth']=WWidth
          Widths['XBWidth']=XBWidth
          Widths['YBWidth']=YBWidth
          Widths['WBWidth']=WBWidth
          Widths['tthwidth_proj']=N.degrees(tthwidth_proj)
          Widths['tthwidth_sec']=N.degrees(tthwidth_sec)
          Widths['thwidth_proj']=N.degrees(thwidth_proj)
          Widths['thwidth_sec']=N.degrees(thwidth_sec)
          return Widths


class TestLattice(unittest.TestCase):

     def setUp(self):
          a=N.array([2*pi],'d')
          b=N.array([2*pi],'d')
          c=N.array([2*pi],'d')
          alpha=N.array([90],'d')
          beta=N.array([90],'d')
          gamma=N.array([90],'d')
          orient1=N.array([[1,0,0]],'d')
          orient2=N.array([[0,1,1]],'d')
          self.fixture = lattice(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,\
                                 orient1=orient1,orient2=orient2)

     def test_astar(self):
          self.assertAlmostEqual(self.fixture.astar[0],1.0,2,'astar Not equal to '+str(1.0))
     def test_bstar(self):
          self.assertAlmostEqual(self.fixture.bstar[0],1.0,2,'bstar Not equal to '+str(1.0))
     def test_cstar(self):
          self.assertAlmostEqual(self.fixture.cstar[0],1.0,2,'cstar '+str(self.fixture.cstar[0])+' Not equal to '+str(1.0))
     def test_alphastar(self):
          self.assertAlmostEqual(self.fixture.alphastar[0],pi/2,2,'alphastar Not equal to '+str(pi/2))
     def test_betastar(self):
          self.assertAlmostEqual(self.fixture.betastar[0],pi/2,2,'betastar Not equal to '+str(pi/2))
     def test_gammastar(self):
          self.assertAlmostEqual(self.fixture.gammastar[0],pi/2,2,'gammastar Not equal to '+str(pi/2))
     def test_V(self):
          self.assertAlmostEqual(self.fixture.V[0],248.0502,2,'V Not equal to '+str(248.0502))
     def test_Vstar(self):
          self.assertAlmostEqual(self.fixture.Vstar[0],1.0,2,'Vstar Not equal to '+str(1.0))
     def test_g(self):
          #print self.fixture.g
          self.assertAlmostEqual((self.fixture.g[:,:,0][0,0]),39.4784*(N.eye(3)[0,0]) ,2,'g Not equal to '+str(39.4784 ))
     def test_gstar(self):
          #print self.fixture.gstar
          self.assertAlmostEqual(self.fixture.gstar[:,:,0][0,0],1.0*N.eye(3)[0,0] ,2,'gstar Not equal to '+str(1.0 ))

     def test_StandardSystem_x(self):
     #       #print self.fixture.gstar
          self.assertAlmostEqual(self.fixture.x[0],1.0 ,2,'Standard System x Not equal to '+str(1.0 ))




#    def test_zeroes(self):
#        self.assertEqual(0 + 0, 0)
#        self.assertEqual(5 + 0, 5)
#        self.assertEqual(0 + 13.2, 13.2)
#
#    def test_positive(self):
#        self.assertEqual(123 + 456, 579)
#        self.assertEqual(1.2e20 + 3.4e20, 3.5e20)
#
#    def test_mixed(self):
#        self.assertEqual(-19 + 20, 1)
#        self.assertEqual(999 + -1, 998)
#        self.assertEqual(-300.1 + -400.2, -700.3)
#


def get_tokenized_line(myfile,returnline=['']):
     lineStr=myfile.readline()
     returnline[0]=lineStr.rstrip()
     strippedLine=lineStr.lower().rstrip()
     tokenized=strippedLine.split()

     return tokenized

if __name__=="__main__":
     
     if 0:
          if 1:
               print 'trying'
               try:
                    infile=sys.argv[1]
               except:
                    infile=r'rescalc_in.txt'
               myfile=open(infile,'r')
               myFlag=True
               returnline=['']
               i=0
               EXP={}
               EXP['ana']={}
               EXP['mono']={}
               EXP['sample']={}
               EXP['method']=0
               while myFlag:
                    toks=get_tokenized_line(myfile,returnline=returnline)
                    #print 'toks',toks
                    if toks==[]:
                         break
                    if toks[0].startswith('#'):
                         continue
                    else:
                         if i==0:
                              a=N.array([float(toks[0])])
                              b=N.array([float(toks[1])])
                              c=N.array([float(toks[2])])
                              alpha=N.radians(N.array([float(toks[3])]))
                              beta=N.radians(N.array([float(toks[4])]))
                              gamma=N.radians(N.array([float(toks[5])]))
                         elif i==1:
                              EXP['hcol']=N.array([float(toks[0]),
                                                   float(toks[1]),
                                                   float(toks[2]),
                                                   float(toks[3])],'d')
                         elif i==2:
                              EXP['vcol']=N.array([float(toks[0]),
                                                   float(toks[1]),
                                                   float(toks[2]),
                                                   float(toks[3])],'d')
                         elif i==3:
                              EXP['mono']['mosaic']=float(toks[0])
                         elif i==4:
                              EXP['ana']['mosaic']=float(toks[0])
                         elif i==5:
                              EXP['sample']['mosaic']=float(toks[0])
                         elif i==6:
                              EXP['sample']['vmosaic']=float(toks[0])
                         elif i==7:
                              EXP['efixed']=float(toks[0])
                         elif i==8:
                              EXP['infix']=1 if toks[0]=='Ei' else -1
                         elif i==9:
                              EXP['mono']['tau']=str(toks[0])
                         elif i==10:
                              EXP['ana']['tau']=str(toks[0])
                         elif i==11:
                              orient1=N.array([toks],'d')
                         elif i==12:
                              orient2=N.array([toks],'d')
                         elif i==13:
                              H=N.array(toks,'d')
                         elif i==14:    
                              K=N.array(toks,'d')
                         elif i==15:     
                              L=N.array(toks,'d')
                         elif i==16:     
                              W=N.array(toks,'d')             
                         i=i+1
               orientation=lattice_calculator.Orientation(orient1,orient2)
               mylattice=lattice_calculator.Lattice(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,\
                                               orientation=orientation)
               setup=[EXP]
               myrescal=rescalculator(mylattice)
               newinput=lattice_calculator.CleanArgs(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,orient1=orient1,orient2=orient2,\
                                                     H=H,K=K,L=L,W=W,setup=setup)
               neworientation=lattice_calculator.Orientation(newinput['orient1'],newinput['orient2'])
               mylattice=lattice_calculator.Lattice(a=newinput['a'],b=newinput['b'],c=newinput['c'],alpha=newinput['alpha'],\
                                                    beta=newinput['beta'],gamma=newinput['gamma'],orientation=neworientation)
               myrescal.__init__(mylattice)
               Q=myrescal.lattice_calculator.modvec(H,K,L,'latticestar')
               print 'Q', Q
               #sys.exit()
               R0,RM=myrescal.ResMat(Q,W,setup)
               print 'RM '
               print RM.transpose()
               print 'R0 ',R0
               #exit()
               #print "shapes"
               #print H.shape
               #print K.shape
               #print L.shape
               #print W.shape
               R0,RMS=myrescal.ResMatS(H,K,L,W,setup)
               myrescal.ResPlot(H, K, L, W, setup)
               print 'RMS'
               print RMS.transpose()[0]
               #print myrescal.calc_correction(H,K,L,W,setup,qscan=[[1,0,1],[1,0,1]])
               print myrescal.calc_correction(H,K,L,W,setup)
               print myrescal.CalcWidths(H,K,L,W,setup)
               print 'setup length ',len(setup)
               sys.exit()
          if 0:
               sys.exit()
          
     if 0:
          a=N.array([2*pi],'d')
          b=N.array([2*pi],'d')
          c=N.array([2*pi],'d')
          alpha=N.array([90],'d')
          beta=N.array([90],'d')
          gamma=N.array([90],'d')
     #       orient1=N.array([[0,1,1]],'d')
          orient1=N.array([[1,0,0]],'d')
          orient2=N.array([[0,1,1]],'d')
          mylattice=lattice_calculator.lattice(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,\
                                               orient1=orient1,orient2=orient2)
          H=N.array([1],'d');K=N.array([0],'d');L=N.array([0],'d');W=N.array([0],'d')
          EXP={}
          EXP['ana']={}
          EXP['ana']['tau']='pg(002)'
          EXP['mono']={}
          EXP['mono']['tau']='pg(002)';
          EXP['ana']['mosaic']=30
          EXP['mono']['mosaic']=30
          EXP['sample']={}
          EXP['sample']['mosaic']=10
          EXP['sample']['vmosaic']=10
          EXP['hcol']=N.array([40, 10, 20, 80],'d')
          EXP['vcol']=N.array([120, 120, 120, 120],'d')
          EXP['infix']=-1 #positive for fixed incident energy
          EXP['efixed']=14.7
          EXP['method']=0
          setup=[EXP]
          myrescal=rescalculator(mylattice)
          #R0,RMS=myrescal.ResMatS(H,K,L,W,setup)
          R0,RM= self.ResMat(Q,W,setup)
          #myrescal.ResPlot(H, K, L, W, setup)
          print 'RMS'
          print RMS.transpose()[0]
          print myrescal.calc_correction(H,K,L,W,setup)
#        mylattice.StandardSystem()
##    x1=N.array([1.0, 1.0], 'd'); y1=N.array([1.0, 1.0], 'd'); z1=N.array([1.0, 1.0], 'd'); x2=x1; y2=y1; z2=z1;
##    print 'scalar ', mylattice.scalar(x1,y1,z1,x2,y2,z2,'lattice')
##    print 'me ',mylattice.gtensor('lattice')
##    myinstrument=instrument();
##    print myinstrument.get_tau('pg(004)')
#    unittest.main()


     if 1:
          a=N.array([6, 6.],'d')
          b=N.array([7., 7.],'d')
          c=N.array([8.,8],'d')
          alpha=N.array([pi/2,pi/2],'d')
          beta=N.array([pi/2,pi/2],'d')
          gamma=N.array([pi/2,pi/2],'d')
     #       orient1=N.array([[0,1,1]],'d')
          orient1=N.array([[1,0,0],[1,0,0]],'d')
          orient2=N.array([[0,0,1],[0,0,1]],'d')
          orientation=lattice_calculator.Orientation(orient1,orient2)
          mylattice=lattice_calculator.Lattice(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,\
                                               orientation=orientation)
          H=N.array([1,1],'d');K=N.array([0,0.0],'d');L=N.array([5.0,5.0],'d');W=N.array([2,2],'d')
          EXP={}
          EXP['ana']={}
          EXP['ana']['tau']='pg(002)'
          EXP['mono']={}
          EXP['mono']['tau']='pg(002)';
          EXP['ana']['mosaic']=25
          EXP['mono']['mosaic']=25
          EXP['sample']={}
          EXP['sample']['mosaic']=25
          EXP['sample']['vmosaic']=25
          EXP['hcol']=N.array([40, 40, 40, 40],'d')
          EXP['vcol']=N.array([120, 120, 120, 120],'d')
          EXP['infix']=-1 #positive for fixed incident energy
          EXP['efixed']=14.7
          EXP['method']=0
          setup=[EXP]
          myrescal=rescalculator(mylattice)
          newinput=lattice_calculator.CleanArgs(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,orient1=orient1,orient2=orient2,\
                                                H=H,K=K,L=L,W=W,setup=setup)
          neworientation=lattice_calculator.Orientation(newinput['orient1'],newinput['orient2'])
          mylattice=lattice_calculator.Lattice(a=newinput['a'],b=newinput['b'],c=newinput['c'],alpha=newinput['alpha'],\
                                               beta=newinput['beta'],gamma=newinput['gamma'],orientation=neworientation)
          myrescal.__init__(mylattice)
          Q=myrescal.lattice_calculator.modvec(H,K,L,'latticestar')
          print 'Q', Q
          R0,RM=myrescal.ResMat(Q,W,setup)
          print 'RM '
          print RM.transpose()
          print 'R0 ',R0
          #exit()
          R0,RMS=myrescal.ResMatS(H,K,L,W,setup)
          myrescal.ResPlot(H, K, L, W, setup)
          print 'RMS'
          print RMS.transpose()[0]
          #print myrescal.calc_correction(H,K,L,W,setup,qscan=[[1,0,1],[1,0,1]])
          print myrescal.calc_correction(H,K,L,W,setup)
          print myrescal.CalcWidths(H,K,L,W,setup)
          print 'setup length ',len(setup)
          
