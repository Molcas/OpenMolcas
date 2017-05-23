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
      subroutine ptwt(abess,arc2,bbess,dfac,npi,l,lambu,ltot1,lmahi,
     &  lmbhi,alpha,ptpow,q2,rc,rka,rkb,prd,hpt,hwt,qsum)
c
c  compute type 2 radial integrals, scaled by exp(-arc2)/sqrt(pi),
c  using the points and weights method,
c  for lama=l to lmahi, lamb=l to lmbhi, n=lama+lamb-l-l
c
      implicit real*8 (a-h,o-z)
      parameter (a500=500.0d0, a50000=50000.0d0)
      dimension abess(*),bbess(*),dfac(*),hpt(*),hwt(*),ptpow(*),
     &  q2(lambu,*),qsum(ltot1,lambu,*)
c
      call wzero(lambu*lmahi,q2,1)
      if(arc2.gt.a50000) then
        npt=5
        idif=0
      else if(arc2.gt.a500) then
        npt=10
        idif=5
      else
        npt=20
        idif=15
      endif
      sqalp=sqrt(alpha)
      prd=prd/sqalp
      do 90 i=1,npt
        pt=rc+hpt(i+idif)/sqalp
        call ssibfn(lmahi-1,rka*pt,abess)
        call ssibfn(lmbhi-1,rkb*pt,bbess)
        if((npi+l+l-2).eq.0) then
          ptpow(1)=prd
        else
          ptpow(1)=prd*pt**(npi+l+l-2)
        endif
        do 70 n=2,ltot1
          ptpow(n)=(pt*pt)*ptpow(n-1)
   70   continue
        do 88 lama=l,lmahi
          do 86 lamb=l,lmbhi
            n=((1-l-l)+lama)+lamb
            q2(lamb,lama)=q2(lamb,lama)+(hwt(i+idif)*abess(lama))*
     1        bbess(lamb)*ptpow(n)
   86     continue
   88   continue
   90 continue
      fctr=rkb**(l-1)
      do 92 lamb=l,lmbhi
        bbess(lamb)=fctr/dfac(lamb+lamb+1)
        fctr=rkb*fctr
   92 continue
      fctr=rka**(l-1)
      do 96 lama=l,lmahi
        do 94 lamb=l,lmbhi
          n=((1-l-l)+lama)+lamb
          qsum(n,lamb,lama)=qsum(n,lamb,lama)+(fctr/dfac(lama+lama+1))*
     1      bbess(lamb)*q2(lamb,lama)
   94   continue
        fctr=rka*fctr
   96 continue
      return
      end
