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
      subroutine rad2_molcas( a,          ccr,        ipt,        kcrl,
     &                 kcru,       l,          lambu,      lmahi,
     &                 lmalo,      lmbhi,      lmblo,      ltot1,
     &                 ncr,        qsum,       rka,        rkb,
     &                 zcr,        lit,        ljt,        ca,
     &                 cb,         tai,        taj,        aa,
     &                 aarr2,      fctr2 )
c
c  compute type 2 radial integrals.
c
c  28-nov-90 new version replaced old version. -rmp
c  19-jan-97 put bessel formula into a separate subroutine, qbess. -rmp
c
      implicit real*8 (a-h,o-z)
      parameter (a0=0.0d0, eps1=1.0d-15, a1=1.0d0, a2=2.0d0,
     1  a4=4.0d0, a50=50.0d0)
      dimension a(*), ccr(*), ipt(*), ncr(*), qsum(ltot1,lambu,*),
     &  zcr(*)
c
      tol=20.*log(10.d0)
      call wzero(ltot1*lambu*lmahi,qsum,1)
c     # sum over potential terms
      do 68 kcr=kcrl,kcru
      npi=ncr(kcr)
      alpha=aa+zcr(kcr)
      rc=(rka+rkb)/(alpha+alpha)
      arc2=alpha*rc*rc
      dum=aarr2+zcr(kcr)*arc2/aa
      if(dum.gt.tol) go to 68
      prd=fctr2*ccr(kcr)*exp(-dum)
c
      if(rka.eq.a0.and.rkb.eq.a0) then
c       # rka=0 and rkb=0
        rk=a0
        t=a0
        qsum(ltot1,1,1)=qsum(ltot1,1,1)+prd*
     &             qcomp(alpha,a(ipt(11)),npi+ltot1-1,0,t,rk)
      elseif(rka.eq.a0) then
c       # rka=0 and rkb>0
        rk=rkb
        t=arc2
        do 12 lamb=l,lmbhi
          qsum(lamb-l+lit,lamb,1)=qsum(lamb-l+lit,lamb,1)+prd*
     1             qcomp(alpha,a(ipt(11)),npi+lamb-l+lit-1,lamb-1,t,rk)
   12   continue
      elseif(rkb.eq.a0) then
c       # rka>0 and rkb=0
        rk=rka
        t=arc2
        do 16 lama=l,lmahi
          qsum(lama-l+ljt,1,lama)=qsum(lama-l+ljt,1,lama)+prd*
     1             qcomp(alpha,a(ipt(11)),npi+lama-l+ljt-1,lama-1,t,rk)
   16   continue
      elseif(npi.eq.2) then
c       # rka>0 and rkb>0; use bessel function formula.
c       # to be applicable for a set of integrals, must have nu.le.l and
c       # nu.eq.(integer), where nu=l+1-npi/2, so it is used here only
c       # for the npi=2 case.  It can't be used at all for npi=(odd) and
c       # only for partial sets for npi=0
        nu=l
            call qbess( alpha,      a(ipt(40)), a(ipt(41)), a(ipt(42)),
     &      a(ipt(12)), a(ipt(43)), a(ipt(44)), a(ipt(45)), a(ipt(11)),
     &      l,          lambu,      lmahi,      lmbhi,      ltot1,
     &      nu,         prd,        qsum,       rka,        rkb,
     &      a(ipt(46)) )
      elseif(arc2.ge.a50) then
c       # rka>0 and rkb>0; use pts and wts method
c       # estimate radial integrals and compare to threshold
        qlim=abs(prd)/(max(a1,(rc+rc)*rka)*max(a1,(rc+rc)*rkb))*
     1    sqrt(a4*(tai+tai)**lit*(taj+taj)**ljt*sqrt(tai*taj)/alpha)
        if(rc.lt.ca) then
          nlim=npi
          qlim=qlim*ca**(lit-1)
        else
          nlim=npi+(lit-1)
        endif
        if(rc.lt.cb) then
          qlim=qlim*cb**(ljt-1)
        else
          nlim=nlim+(ljt-1)
        endif
        if(qlim*rc**nlim.ge.eps1) then
          call ptwt(a(ipt(47)),arc2,a(ipt(48)),a(ipt(11)),npi,l,lambu,
     &      ltot1,lmahi,lmbhi,alpha,a(ipt(49)),a(ipt(50)),rc,rka,rkb,
     &      prd,a(ipt(9)),a(ipt(10)),qsum)
        endif
      else
c       # rka>0 and rkb>0; use partially asymptotic method
        call qpasy(alpha,a(ipt(11)),npi,l,lambu,lmahi,lmbhi,ltot1,rka,
     &    rkb,fctr2*ccr(kcr),dum+arc2,qsum)
      endif
   68 continue
c
      if(rka.eq.a0.and.rkb.ne.a0) then
c       # rka=0 and rkb>0
        f2lmb3=dble(2*lmbhi-3)
        do 76 lamb=lmbhi-2,lmblo,-1
          nlo=abs(lamb-l+1)+lit+1
          nhi=ljt-mod((ljt-1)-abs(lamb-l),2)+lit-1
          do 72 n=nlo,nhi,2
            qsum(n,lamb,1)=qsum(n,lamb+2,1)+
     &                     (f2lmb3/rkb)*qsum(n-1,lamb+1,1)
   72     continue
          f2lmb3=f2lmb3-a2
   76   continue
      elseif(rka.ne.a0.and.rkb.eq.a0) then
c       # rka>0 and rkb=0
        f2lma3=dble(2*lmahi-3)
        do 84 lama=lmahi-2,lmalo,-1
          nlo=abs(lama-l+1)+ljt+1
          nhi=lit-mod((lit-1)-abs(lama-l),2)+ljt-1
          do 80 n=nlo,nhi,2
            qsum(n,1,lama)=qsum(n,1,lama+2)+
     &                     (f2lma3/rka)*qsum(n-1,1,lama+1)
   80     continue
          f2lma3=f2lma3-a2
   84   continue
      elseif(rka.ne.a0.and.rkb.ne.a0) then
c       # rka>0 and rkb>0
        f2lma3=dble(lmahi+lmahi+1)
        do 96 lama=lmahi,lmalo,-1
          ldifa1=abs(l-lama)+1
          f2lmb3=dble(2*lmbhi+1)
          do 92 lamb=lmbhi,lmblo,-1
            ldifb=abs(l-lamb)
            nlo=ldifa1+ldifb
            nhi=(ltot1-mod(lit-ldifa1,2))-mod((ljt-1)-ldifb,2)
            do 88 n=nlo,nhi,2
              if(n-(lama+lamb).eq.(1-l-l)) go to 88
              if(lama.gt.(lmahi-2).or.n.le.(abs(l-lama-2)+ldifb)) then
c               # lamb recursion
                qsum(n,lamb,lama)=qsum(n,lamb+2,lama)+
     1                            (f2lmb3/rkb)*qsum(n-1,lamb+1,lama)
              else
c               # lama recursion
                qsum(n,lamb,lama)=qsum(n,lamb,lama+2)+
     1                            (f2lma3/rka)*qsum(n-1,lamb,lama+1)
              endif
   88       continue
            f2lmb3=f2lmb3-a2
   92     continue
          f2lma3=f2lma3-a2
   96   continue
      endif
      return
      end
