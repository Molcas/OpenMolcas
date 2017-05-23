************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine qpasy(alpha,dfac,npi,l,lambu,lmahi,lmbhi,ltot1,xka,
     &  xkb,prd,dum,qsum)
c
c  compute type 2 radial integrals, scaled by exp(-arc2)/sqrt(pi),
c  using the partially asymptotic method,
c  for lama=l to lmahi, lamb=l to lmbhi, n=lama+lamb-l-l
c
      implicit real*8 (a-h,o-z)
      parameter (a0=0.0d0, accrcy=1.0d-13, a1s4=0.25d0, a1s2=0.5d0,
     1  a1=1.0d0)
      dimension dfac(*),qsum(ltot1,lambu,*)
c
      sqalpi=a1/sqrt(alpha)
      alf1=a1
      if(xka.gt.xkb) go to 42
c
c  xka is smaller: set up parameters for qcomp using xkb
c
      xk=xkb*sqalpi
      t=a1s4*xk*xk
      prde=prd*exp(t-dum)*sqalpi**(npi+l)
      if(l.ge.2) prde=prde*xka**(l-1)
      tk=xka*xka/(alpha+alpha)
      do 30 lama=l,lmahi
      la=lama-1
      prefac=prde
      do 28 lamb=l,lmbhi
      lb=lamb-1
      n=((1-l-l)+lama)+lamb
c     # run power series using xka, obtaining initial
c     # q(n,l) values from qcomp, then recurring upwards
c     # j=0 term in sum
      nprime=npi+n+la-1
      qold1=qcomp(alf1,dfac,nprime,lb,t,xk)/dfac(la+la+3)
      sum=qold1
      if(tk.eq.a0) go to 24
c     # j=1 term in sum
      nprime=nprime+2
      qnew=qcomp(alf1,dfac,nprime,lb,t,xk)/dfac(la+la+3)
      f1=(la+la+3)
      qold2=(tk/f1)*qold1
      qold1=(tk/f1)*qnew
      sum=sum+qold1
      j=1
c     # increment j for next term
   22 j=j+1
      nprime=nprime+2
      f1=(nprime+nprime-5)
      f2=((lb-nprime+4)*(lb+nprime-3))
      qnew=(t+a1s2*f1)*qold1+a1s4*f2*qold2
      f1=(j*(la+la+j+j+1))
      qold2=(tk/f1)*qold1
      qold1=(tk/f1)*qnew
      sum=sum+qold1
      if(qold1.gt.accrcy*sum) go to 22
   24 qsum(n,lamb,lama)=qsum(n,lamb,lama)+prefac*sum
      prefac=prefac*sqalpi
   28 continue
      prde=prde*(xka/alpha)
   30 continue
      return
c
c  xkb is smaller: set up parameters for qcomp using xka
c
   42 xk=xka*sqalpi
      t=a1s4*xk*xk
      prde=prd*exp(t-dum)*sqalpi**(npi+l)
      if(l.ge.2) prde=prde*xkb**(l-1)
      tk=xkb*xkb/(alpha+alpha)
      do 60 lama=l,lmahi
      la=lama-1
      prefac=prde
      do 58 lamb=l,lmbhi
      lb=lamb-1
      n=((1-l-l)+lama)+lamb
c     # run power series using xkb, obtaining initial
c     # q(n,l) values from qcomp, then recurring upwards
c     # j=0 term in sum
      nprime=npi+n+lb-1
      qold1=qcomp(alf1,dfac,nprime,la,t,xk)/dfac(lb+lb+3)
      sum=qold1
      if(tk.eq.a0) go to 54
c     # j=1 term in sum
      nprime=nprime+2
      qnew=qcomp(alf1,dfac,nprime,la,t,xk)/dfac(lb+lb+3)
      f1=(lb+lb+3)
      qold2=(tk/f1)*qold1
      qold1=(tk/f1)*qnew
      sum=sum+qold1
      j=1
c     # increment j for next term
   52 j=j+1
      nprime=nprime+2
      f1=(nprime+nprime-5)
      f2=((la-nprime+4)*(la+nprime-3))
      qnew=(t+a1s2*f1)*qold1+a1s4*f2*qold2
      f1=(j*(lb+lb+j+j+1))
      qold2=(tk/f1)*qold1
      qold1=(tk/f1)*qnew
      sum=sum+qold1
      if(qold1.gt.accrcy*sum) go to 52
   54 qsum(n,lamb,lama)=qsum(n,lamb,lama)+prefac*sum
      prefac=prefac*(xkb/alpha)
   58 continue
      prde=prde*sqalpi
   60 continue
      return
      end
