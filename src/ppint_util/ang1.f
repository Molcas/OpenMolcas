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
      subroutine ang1(ang,dfac,nanb,lalb,mamb,lamu,lmf,lml,lmx,lmy,lmz,
     &  ltot1,xab,yab,zab,xk,yk,zk,zlm)
c
c  compute type 1 angular integrals
c
      implicit real*8 (a-h,o-z)
      parameter (a0=0.0d0, a1=1.0d0)
      dimension ang(ltot1,*), dfac(*), lmf(*), lml(*), lmx(*), lmy(*),
     &  lmz(*), xab(*), yab(*), zab(*), zlm(*)
c
      call wzero(ltot1*lamu,ang,1)
      do 96 n=1,nanb
      if(xab(n).eq.a0) go to 96
      do 94 l=1,lalb
      if(yab(l).eq.a0) go to 94
      do 92 m=1,mamb
      if(zab(m).eq.a0) go to 92
      nlm=((n-2)+l)+m
      lamlo=mod(nlm-1,2)+1
      lamhi=min(nlm,lamu)
      if(lamlo.gt.lamhi) go to 92
      do 90 lam=lamlo,lamhi,2
      l2=lam+lam-1
      angt=a0
      loc=(lam-1)**2
      do 80 mu1=1,l2
      istart=lmf(loc+mu1)
      if(mod(n,2).eq.mod(lmx(istart),2).or.
     1   mod(l,2).eq.mod(lmy(istart),2).or.
     2   mod(m,2).eq.mod(lmz(istart),2)) go to 80
      pre=a0
      aint=a0
      iend=lml(loc+mu1)
      do 70 i=istart,iend
        indx=lmx(i)
        indy=lmy(i)
        indz=lmz(i)
        if(indx.eq.0) then
          xkp=a1
        else
          xkp=xk**indx
        endif
        if(indy.eq.0) then
          ykp=a1
        else
          ykp=yk**indy
        endif
        if(indz.eq.0) then
          zkp=a1
        else
          zkp=zk**indz
        endif
        pre=pre+zlm(i)*xkp*ykp*zkp
        aint=aint+zlm(i) * dfac(n+indx) * dfac(l+indy) * dfac(m+indz)/
     1    dfac((n+indx)+(l+indy)+(m+indz))
   70 continue
      angt=angt+pre*aint
   80 continue
      ang(nlm,lam)=ang(nlm,lam)+((xab(n)*yab(l))*zab(m))*angt
   90 continue
   92 continue
   94 continue
   96 continue
      return
      end
