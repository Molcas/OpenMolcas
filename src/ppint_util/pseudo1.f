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
      subroutine pseud1_molcas(a,ang,ccr,gout,ipt,lmnv,ltot1,ncr,
     &  nkcrl,nkcru,qsum,xab,yab,zab,zcr,lit,ljt,ai,aj,
     &  xi,yi,zi,xj,yj,zj,xc,yc,zc,kcrs,lproju1,crda,crdb)
c
c  compute type 1 core potential integrals
c
      implicit real*8 (a-h,o-z)
      dimension llt(7,2)
      parameter (a0=0.0d0,a4=4.0d0)
      dimension a(*),ang(ltot1,*),ccr(*),crda(lit,3),crdb(ljt,3),
     &  gout(*),ipt(*),lmnv(3,*),ncr(*),nkcrl(lproju1,*),
     &  nkcru(lproju1,*),qsum(ltot1,*),xab(*),yab(*),zab(*),zcr(*)
      data llt /1,2,5,11,21,36,57,1,4,10,20,35,56,84/

      tol=20.*log(10.d0)
      fctr2=a4
      itl=llt(lit,1)
      itu=llt(lit,2)
      jtl=llt(ljt,1)
      jtu=llt(ljt,2)
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
      xij=0.5d0*(xi+xj)
      yij=0.5d0*(yi+yj)
      zij=0.5d0*(zi+zj)
      xijm=0.5d0*(xi-xj)
      yijm=0.5d0*(yi-yj)
      zijm=0.5d0*(zi-zj)
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      aaa=(ai-aj)/aa
      xk=xij+aaa*xijm-xc
      yk=yij+aaa*yijm-yc
      zk=zij+aaa*zijm-zc
      aarr1=(ai*aj/aa)*rr
      taa=aa+aa
c*
      rp2=xk*xk+yk*yk+zk*zk
      if(rp2.eq.a0) then
        rp=a0
        arp2=a0
        alpt=a0
        rk=a0
        lamu=1
      else
        rp=sqrt(rp2)
        xk=xk/rp
        yk=yk/rp
        zk=zk/rp
        arp2=aa*rp2
        alpt=aa*arp2
        rk=taa*rp
        lamu=ltot1
      endif
c
c  compute radial integrals and sum over potential terms
c
      call wzero((ltot1*lamu),qsum,1)
      kcrl=nkcrl(1,kcrs)
      kcru=nkcru(1,kcrs)
      call rad1(aa,aarr1,alpt,arp2,ccr,
     &  a(ipt(11)),fctr2,kcrl,kcru,
     &  lamu,ltot1,ncr,qsum,rk,tol,zcr)
      ijt=0
      do 90 it=itl,itu
        na1=lmnv(1,it)+1
        la1=lmnv(2,it)+1
        ma1=lmnv(3,it)+1
        do 80 jt=jtl,jtu
          ijt=ijt+1
          s=a0
          nb1=lmnv(1,jt)+1
          lb1=lmnv(2,jt)+1
          mb1=lmnv(3,jt)+1
c
c  compute angular integrals
c
          call facab(a(ipt(12)),na1,nb1,crda(1,1),crdb(1,1),xab)
          call facab(a(ipt(12)),la1,lb1,crda(1,2),crdb(1,2),yab)
          call facab(a(ipt(12)),ma1,mb1,crda(1,3),crdb(1,3),zab)
          call ang1(ang,a(ipt(11)),na1+nb1-1,la1+lb1-1,ma1+mb1-1,lamu,
     &      a(ipt(13)),a(ipt(14)),a(ipt(15)),a(ipt(16)),a(ipt(17)),
     &      ltot1,xab,yab,zab,xk,yk,zk,a(ipt(18)))
c
c  combine angular and radial integrals
c
          do 70 lam=1,lamu
            nhi=ltot1-mod(ltot1-lam,2)
            do 60 n=lam,nhi,2
              s=s+ang(n,lam)*qsum(n,lam)
   60       continue
   70     continue
          gout(ijt)=gout(ijt)+s
   80   continue
   90 continue
      return
      end
