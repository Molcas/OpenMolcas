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
      function qcomp(alpha,dfac,n,l,t,xk)
c
c  compute q(n,l) scaled by sqrt(pi)*exp(-t) to prevent overflows
c  arguments are alpha, xk, and t=xk**2/(4*alpha)
c  no restriction on the magnitude of t
c  increase dfac array to raise n, l restrictions
c
      implicit real*8 (a-h,o-z)
      parameter (am1=-1.0d0, a0=0.0d0, accpow=1.0d-14,
     1  accasy=1.0d-10, a1rtpi=0.56418958354775629d0, a1=1.0d0,
     2  a2=2.0d0, a4=4.0d0)
      dimension dfac(*)
      dimension tmin(0:8)
      data tmin/31.0d0,28.0d0,25.0d0,23.0d0,22.0d0,20.0d0,19.0d0,
     1  18.0d0,15.0d0/
c
      if(mod(n+l,2).ne.0.or.n.le.l) go to 30
c
c  use alternating series (n+l.le.22.and.l.le.10)
c
      if(l.eq.0) then
        xkp=a1
      else
        xkp=(xk/(alpha+alpha))**l
      endif
      prefac=xkp * dfac(n+l+1)/
     1  ((alpha+alpha)**((n-l)/2) * sqrt(a4*alpha) * dfac(l+l+3))
      num=l-n+2
      xden=(l+l+3)
      term=a1
      sum=term
      xc=am1
   10 if(num.ne.0) then
        fnum=num
        term=term*fnum*t/(xden*xc)
        xc=xc+am1
        sum=sum+term
        num=num+2
        xden=xden+a2
        go to 10
      endif
c
      qcomp=prefac*sum
      return
c
   30 if(t.lt.tmin(min(n,8))) go to 60
c
c  use asymptotic series (arbitrary n,l)
c
      xkp=(xk/(alpha+alpha))**(n-2)
      prefac=xkp/((alpha+alpha)*sqrt(a4*alpha))
      sum=a1
      term=a1
      fac1=(l-n+2)
      fac2=(1-l-n)
      xc=a1
   40 term=term*fac1*fac2/(a4*xc*t)
      if(term.eq.a0) go to 50
      sum=sum+term
      if(abs(term/sum).lt.accasy) go to 50
      fac1=fac1+a2
      fac2=fac2+a2
      xc=xc+a1
      go to 40
   50 qcomp=prefac*sum
      return
c
c  use power series (n+l.le.22.and.l.le.10)
c
   60 if(l.eq.0) then
        xkp=a1
      else
        xkp=(xk/(alpha+alpha))**l
      endif
      prefac=exp(-t)*xkp/(alpha+alpha)**((n-l+1)/2)
      if(mod(n+l,2).eq.0) then
        prefac=prefac/sqrt(a4*alpha)
      else
        prefac=a1rtpi*prefac
      endif
      xnum=(l+n-1)
      xden=(l+l+1)
      term=dfac(l+n+1)/dfac(l+l+3)
      sum=term
      xj=a0
   70 xnum=xnum+a2
      xden=xden+a2
      xj=xj+a1
      term=term*t*xnum/(xj*xden)
      sum=sum+term
      if((term/sum).gt.accpow) go to 70
      qcomp=prefac*sum
      return
      end
