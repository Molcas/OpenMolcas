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
      subroutine pseud2_molcas( a, anga, angb, ccr, gout, ipt,
     &  lambu, ltot1, mproju, ncr, nkcrl, nkcru, qsum,
     &  zcr, lit, ljt, ai, aj, xi, yi, zi, xj, yj, zj, xc, yc, zc,
     &  kcrs, lcru, lproju1,crda,crdb)
      implicit real*8 (a-h,o-z)
      dimension llt(7,2)
c
c  compute type 2 core potential integrals
c
      dimension a(*),anga(lit,mproju,*),angb(ljt,mproju,*),ccr(*),
     &  crda(lit,3),crdb(ljt,3),gout(*),ipt(*),ncr(*),
     &  nkcrl(lproju1,*),nkcru(lproju1,*),qsum(ltot1,lambu,*),zcr(*)
      data llt /1,2,5,11,21,36,57,1,4,10,20,35,56,84/
      data      a0   ,a4     /        0.0d0,      4.0d0       /
*     data eps1,a0,a1,a4,a50 /1.0d-15,0.0d0,1.0d0,4.0d0,50.0d0/

      call pseud2_molcas_internal(a)
*
*     This is to allow type punning without an explicit interface
      contains
      subroutine pseud2_molcas_internal(a)
      use iso_c_binding
      real*8, target :: a(*)
      integer, pointer :: ia1(:),ia13(:),ia14(:),ia15(:),ia16(:),ia17(:)
      tol=20.*log(10.d0)
      fctr2=a4
      itl=llt(lit,1)
      itu=llt(lit,2)
      jtl=llt(ljt,1)
      jtu=llt(ljt,2)
      tai=ai+ai
      taj=aj+aj
      aa=ai+aj
      do 120 i=1,3
      crda(1,i)=1.d0
  120 crdb(1,i)=1.d0
      xka=xc-xi
      yka=yc-yi
      zka=zc-zi
      ca=sqrt(xka*xka+yka*yka+zka*zka)
      if(lit.eq.1) go to 220
      crda(2,1)=xka
      crda(2,2)=yka
      crda(2,3)=zka
      if(lit.eq.2) go to 220
      do 210 i=1,3
      do 210 j=3,lit
  210 crda(j,i)=crda(2,i)*crda(j-1,i)
  220 xkb=xc-xj
      ykb=yc-yj
      zkb=zc-zj
      cb=sqrt(xkb*xkb+ykb*ykb+zkb*zkb)
      if(ljt.eq.1) go to 240
      crdb(2,1)=xkb
      crdb(2,2)=ykb
      crdb(2,3)=zkb
      if(ljt.eq.2) go to 240
      do 230 i=1,3
      do 230 j=3,ljt
  230 crdb(j,i)=crdb(2,i)*crdb(j-1,i)
  240 continue
      aarr2=(ai*aj/aa)*(ca-cb)**2
c*
      if(ca.eq.a0) then
        rka=a0
        lmau=1
      else
        xka=-xka/ca
        yka=-yka/ca
        zka=-zka/ca
        rka=tai*ca
        lmau=lcru+(lit-1)
      endif
      if(cb.eq.a0) then
        rkb=a0
        lmbu=1
      else
        xkb=-xkb/cb
        ykb=-ykb/cb
        zkb=-zkb/cb
        rkb=taj*cb
        lmbu=lcru+(ljt-1)
      endif
      if((ca.eq.a0).and.(cb.eq.a0)) then
        lhi=min(lcru,lit,ljt)
        llo=mod((lit-1),2)+1
        if(llo.ne.mod((ljt-1),2)+1.or.llo.gt.lhi) return
        inc=2
      elseif(ca.eq.a0) then
        lhi=min(lcru,lit)
        llo=mod((lit-1),2)+1
        if(llo.gt.lhi) return
        inc=2
      elseif(cb.eq.a0) then
        lhi=min(lcru,ljt)
        llo=mod((ljt-1),2)+1
        if(llo.gt.lhi) return
        inc=2
      else
        lhi=lcru
        llo=1
        inc=1
      endif
      do 88 l=llo,lhi,inc
      mhi=l+l-1
      lmalo=max(l-(lit-1),1)
      lmahi=min(lmau,l+(lit-1))
      lmblo=max(l-(ljt-1),1)
      lmbhi=min(lmbu,l+(ljt-1))
c
c  compute radial integrals
c
      kcrl=nkcrl(l+1,kcrs)
      kcru=nkcru(l+1,kcrs)
            call rad2_molcas( a,          ccr,        ipt,        kcrl,
     &                 kcru,       l,          lambu,      lmahi,
     &                 lmalo,      lmbhi,      lmblo,      ltot1,
     &                 ncr,        qsum,       rka,        rkb,
     &                 zcr,        lit,        ljt,        ca,
     &                 cb,         tai,        taj,        aa,
     &                 aarr2,      fctr2 )
c
c  compute angular integrals and combine with radial integrals
c
      ijt=0
      do 84 it=itl,itu
      call c_f_pointer(c_loc(a(ipt(1))),ia1,[1])
      call c_f_pointer(c_loc(a(ipt(1))),ia1,[1])
      call c_f_pointer(c_loc(a(ipt(13))),ia13,[1])
      call c_f_pointer(c_loc(a(ipt(14))),ia14,[1])
      call c_f_pointer(c_loc(a(ipt(15))),ia15,[1])
      call c_f_pointer(c_loc(a(ipt(16))),ia16,[1])
      call c_f_pointer(c_loc(a(ipt(17))),ia17,[1])
      call ang2_molcas(anga,a(ipt(12)),crda,a(ipt(11)),it,
     &  l,lit,lmalo,lmahi,
     &  ia13,ia14,ia15,ia16,ia17,
     &  ia1,mproju,xka,yka,zka,a(ipt(18)))
      nullify(ia1,ia13,ia14,ia15,ia16,ia17)
      do 80 jt=jtl,jtu
      ijt=ijt+1
      s=a0
      call c_f_pointer(c_loc(a(ipt(1))),ia1,[1])
      call c_f_pointer(c_loc(a(ipt(1))),ia1,[1])
      call c_f_pointer(c_loc(a(ipt(13))),ia13,[1])
      call c_f_pointer(c_loc(a(ipt(14))),ia14,[1])
      call c_f_pointer(c_loc(a(ipt(15))),ia15,[1])
      call c_f_pointer(c_loc(a(ipt(16))),ia16,[1])
      call c_f_pointer(c_loc(a(ipt(17))),ia17,[1])
      call ang2_molcas(angb,a(ipt(12)),crdb,a(ipt(11)),jt,
     &  l,ljt,lmblo,lmbhi,
     &  ia13,ia14,ia15,ia16,ia17,
     &  ia1,mproju,xkb,ykb,zkb,a(ipt(18)))
      nullify(ia1,ia13,ia14,ia15,ia16,ia17)
      do 76 lama=lmalo,lmahi
        ldifa1=abs(l-lama)+1
        nlmau=lit-mod(lit-ldifa1,2)
        do 72 lamb=lmblo,lmbhi
          ldifb=abs(l-lamb)
          nlmbu=(ljt-1)-mod((ljt-1)-ldifb,2)
          nlo=ldifa1+ldifb
          nhi=nlmau+nlmbu
          do 68 n=nlo,nhi,2
            nlmalo=max(ldifa1,n-nlmbu)
            nlmahi=min(nlmau,n-ldifb)
            angp=a0
            do 60 m=1,mhi
              do 56 nlma=nlmalo,nlmahi,2
                angp=angp+anga(nlma,m,lama)*angb((n+1)-nlma,m,lamb)
   56         continue
   60       continue
            s=s+angp*qsum(n,lamb,lama)
   68     continue
   72   continue
   76 continue
      gout(ijt)=gout(ijt)+s
   80 continue
   84 continue
   88 continue
      return
      end subroutine pseud2_molcas_internal
*
      end
