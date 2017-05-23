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
      subroutine rad1(aa,aarr1,alpt,arp2,ccr,dfac,fctr2,kcrl,kcru,
     &  lamu,ltot1,ncr,qsum,rk,tol,zcr)
c
c  compute type 1 radial integrals
c
      implicit real*8 (a-h,o-z)
      parameter (a0=0.0d0, a2=2.0d0)
      dimension ccr(*), dfac(*), ncr(*), qsum(ltot1,*), zcr(*)
c
      do 40 kcr=kcrl,kcru
        npi=ncr(kcr)
        alpha=aa+zcr(kcr)
c       # exponential factor from q functions included in dum
        dum=aarr1+zcr(kcr)*arp2/alpha
        if(dum.gt.tol) go to 40
        prd=fctr2*ccr(kcr)*exp(-dum)
        if(rk.eq.a0) then
          t=a0
          do 20 n=1,ltot1-mod(ltot1-1,2),2
            qsum(n,1)=qsum(n,1)+prd*
     &                qcomp(alpha,dfac,npi+n-1,0,t,rk)
   20     continue
        else
          t=alpt/alpha
          do 30 lam=1,lamu
            qsum(lam,lam)=qsum(lam,lam)+prd*
     &                    qcomp(alpha,dfac,npi+lam-1,lam-1,t,rk)
   30     continue
        endif
   40 continue
c
      if(rk.ne.a0) then
        f2lam3=(lamu+lamu-3)
        do 60 lam=lamu-2,1,-1
          do 50 n=lam+2,lamu-mod(lamu-lam,2),2
            qsum(n,lam)=qsum(n,lam+2)+
     &                  (f2lam3/rk)*qsum(n-1,lam+1)
   50     continue
          f2lam3=f2lam3-a2
   60   continue
      endif
      return
      end
