************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996-2006, T. Thorsteinsson and D. L. Cooper           *
************************************************************************
      subroutine dev2a_2_cvb(v1,v2,cfrom,hessorb,oaa2,aa1,nprorb,
     > i1alf,i1bet,iafrm,ibfrm,iato,ibto,phato,phbto,
     > iapr,ixapr,ibpr,ixbpr,npvb,
     > nda,ndb,n1a,n1b,nam1,nbm1,norb,commut,sc,absym)
c  Calculates V1 EijEkl CFROM and V2 EijEkl CFROM
      implicit real*8 (a-h,o-z)
      logical commut,sc,absym
      dimension v1(nda,ndb),v2(nda,ndb),cfrom(nda,ndb)
      dimension hessorb(nprorb,nprorb)
      dimension i1alf(n1a,norb),i1bet(n1b,norb)
      dimension iafrm(norb,nda),ibfrm(norb,ndb)
      dimension iato(norb,0:nam1),ibto(norb,0:nbm1)
      dimension phato(norb,nam1),phbto(norb,nbm1)
      dimension iapr(npvb),ixapr(nda+1),ibpr(npvb),ixbpr(ndb+1)
      save zero,two
      data zero/0.0d0/,two/2d0/

      do 9379 iprm1=1,nprorb
      i1=(iprm1-1)/(norb-1)+1
      i2=iprm1-(i1-1)*(norb-1)
      if(i2.ge.i1) i2=i2+1
      do 9379 iprm2=1,iprm1
      i3=(iprm2-1)/(norb-1)+1
      i4=iprm2-(i3-1)*(norb-1)
      if(i4.ge.i3) i4=i4+1
      if(i1.ne.i4 .or. i2.eq.i3)then
        iorb=i3
        jorb=i4
        korb=i1
        lorb=i2
       else
        iorb=i1
        jorb=i2
        korb=i3
        lorb=i4
      endif
      res1=zero
      res2=zero
      if(commut.and..not.sc)then
c  1) Alpha excitation
      do 1100 ia=1,n1a
      iaxtmp=i1alf(ia,iorb)
      jax=iato(jorb,iaxtmp)
      if(jax.ne.0)then
        iax=iato(iorb,iaxtmp)
        tcof=phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
c  I -> J
        do 1200 ixa=ixapr(iax),ixapr(iax+1)-1
        ibx=iapr(ixa)
c 2. alpha k -> l
        itmp=iafrm(korb,jax)
        kax=iato(lorb,itmp)
        if(kax.ne.0)then
          phase=phato(korb,itmp)*phato(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1+v1(kax,ibx)*term
          res2=res2+v2(kax,ibx)*term
        endif
c 2. alpha l -> k
        itmp=iafrm(lorb,iax)
        lax=iato(korb,itmp)
        if(lax.ne.0)then
          phase=phato(lorb,itmp)*phato(korb,itmp)
          term=tcof*phase*cfrom(lax,ibx)
          res1=res1-v1(jax,ibx)*term
          res2=res2-v2(jax,ibx)*term
        endif
c 2. beta  k -> l
        itmp=ibfrm(korb,ibx)
        kbx=ibto(lorb,itmp)
        if(kbx.ne.0)then
          phase=phbto(korb,itmp)*phbto(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1+v1(jax,kbx)*term
          res2=res2+v2(jax,kbx)*term
        endif
c 2. beta  l -> k
        itmp=ibfrm(lorb,ibx)
        lbx=ibto(korb,itmp)
        if(lbx.ne.0)then
          phase=phbto(lorb,itmp)*phbto(korb,itmp)
          term=tcof*phase*cfrom(iax,lbx)
          res1=res1-v1(jax,ibx)*term
          res2=res2-v2(jax,ibx)*term
        endif
1200    continue
c  J -> I
        do 1300 ixa=ixapr(jax),ixapr(jax+1)-1
        ibx=iapr(ixa)
c 2. alpha k -> l
        itmp=iafrm(korb,jax)
        kax=iato(lorb,itmp)
        if(kax.ne.0)then
          phase=phato(korb,itmp)*phato(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1-v1(kax,ibx)*term
          res2=res2-v2(kax,ibx)*term
        endif
c 2. alpha l -> k
        itmp=iafrm(lorb,iax)
        lax=iato(korb,itmp)
        if(lax.ne.0)then
          phase=phato(lorb,itmp)*phato(korb,itmp)
          term=tcof*phase*cfrom(lax,ibx)
          res1=res1+v1(jax,ibx)*term
          res2=res2+v2(jax,ibx)*term
        endif
c 2. beta  k -> l
        itmp=ibfrm(korb,ibx)
        kbx=ibto(lorb,itmp)
        if(kbx.ne.0)then
          phase=phbto(korb,itmp)*phbto(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1-v1(jax,kbx)*term
          res2=res2-v2(jax,kbx)*term
        endif
c 2. beta  l -> k
        itmp=ibfrm(lorb,ibx)
        lbx=ibto(korb,itmp)
        if(lbx.ne.0)then
          phase=phbto(lorb,itmp)*phbto(korb,itmp)
          term=tcof*phase*cfrom(iax,lbx)
          res1=res1+v1(jax,ibx)*term
          res2=res2+v2(jax,ibx)*term
        endif
1300    continue
      endif
1100  continue

      if(absym)then
        res1=two*res1
        res2=two*res2
      else
c  2) Beta excitation
      do 2100 ib=1,n1b
      ibxtmp=i1bet(ib,iorb)
      jbx=ibto(jorb,ibxtmp)
      if(jbx.ne.0)then
        ibx=ibto(iorb,ibxtmp)
        tcof=phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
c  I -> J
        do 2200 ixb=ixbpr(ibx),ixbpr(ibx+1)-1
        iax=ibpr(ixb)
c 2. beta  k -> l
        itmp=ibfrm(korb,jbx)
        kbx=ibto(lorb,itmp)
        if(kbx.ne.0)then
          phase=phbto(korb,itmp)*phbto(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1+v1(iax,kbx)*term
          res2=res2+v2(iax,kbx)*term
        endif
c 2. beta  l -> k
        itmp=ibfrm(lorb,ibx)
        lbx=ibto(korb,itmp)
        if(lbx.ne.0)then
          phase=phbto(lorb,itmp)*phbto(korb,itmp)
          term=tcof*phase*cfrom(iax,lbx)
          res1=res1-v1(iax,jbx)*term
          res2=res2-v2(iax,jbx)*term
        endif
c 2. alpha k -> l
        itmp=iafrm(korb,iax)
        kax=iato(lorb,itmp)
        if(kax.ne.0)then
          phase=phato(korb,itmp)*phato(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1+v1(kax,jbx)*term
          res2=res2+v2(kax,jbx)*term
        endif
c 2. alpha l -> k
        itmp=iafrm(lorb,iax)
        lax=iato(korb,itmp)
        if(lax.ne.0)then
          phase=phato(lorb,itmp)*phato(korb,itmp)
          term=tcof*phase*cfrom(lax,ibx)
          res1=res1-v1(iax,jbx)*term
          res2=res2-v2(iax,jbx)*term
        endif
2200    continue
c  J -> I
        do 2300 ixb=ixbpr(jbx),ixbpr(jbx+1)-1
        iax=ibpr(ixb)
c 2. beta  k -> l
        itmp=ibfrm(korb,jbx)
        kbx=ibto(lorb,itmp)
        if(kbx.ne.0)then
          phase=phbto(korb,itmp)*phbto(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1-v1(iax,kbx)*term
          res2=res2-v2(iax,kbx)*term
        endif
c 2. beta  l -> k
        itmp=ibfrm(lorb,ibx)
        lbx=ibto(korb,itmp)
        if(lbx.ne.0)then
          phase=phbto(lorb,itmp)*phbto(korb,itmp)
          term=tcof*phase*cfrom(iax,lbx)
          res1=res1+v1(iax,jbx)*term
          res2=res2+v2(iax,jbx)*term
        endif
c 2. alpha k -> l
        itmp=iafrm(korb,iax)
        kax=iato(lorb,itmp)
        if(kax.ne.0)then
          phase=phato(korb,itmp)*phato(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1-v1(kax,jbx)*term
          res2=res2-v2(kax,jbx)*term
        endif
c 2. alpha l -> k
        itmp=iafrm(lorb,iax)
        lax=iato(korb,itmp)
        if(lax.ne.0)then
          phase=phato(lorb,itmp)*phato(korb,itmp)
          term=tcof*phase*cfrom(lax,ibx)
          res1=res1+v1(iax,jbx)*term
          res2=res2+v2(iax,jbx)*term
        endif
2300    continue
      endif
2100  continue
      endif
      elseif((.not.commut).and..not.sc)then
c  1) Alpha excitation
      do 3100 ia=1,n1a
      iaxtmp=i1alf(ia,iorb)
      jax=iato(jorb,iaxtmp)
      if(jax.ne.0)then
        iax=iato(iorb,iaxtmp)
        tcof=phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
c  I -> J
        do 3200 ixa=ixapr(iax),ixapr(iax+1)-1
        ibx=iapr(ixa)
c 2. alpha k -> l
        itmp=iafrm(korb,jax)
        kax=iato(lorb,itmp)
        if(kax.ne.0)then
          phase=phato(korb,itmp)*phato(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1+v1(kax,ibx)*term
          res2=res2+v2(kax,ibx)*term
        endif
c 2. beta  k -> l
        itmp=ibfrm(korb,ibx)
        kbx=ibto(lorb,itmp)
        if(kbx.ne.0)then
          phase=phbto(korb,itmp)*phbto(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1+v1(jax,kbx)*term
          res2=res2+v2(jax,kbx)*term
        endif
3200    continue
      endif
3100  continue

      if(absym)then
        res1=two*res1
        res2=two*res2
      else
c  2) Beta excitation
      do 4100 ib=1,n1b
      ibxtmp=i1bet(ib,iorb)
      jbx=ibto(jorb,ibxtmp)
      if(jbx.ne.0)then
        ibx=ibto(iorb,ibxtmp)
        tcof=phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
c  I -> J
        do 4200 ixb=ixbpr(ibx),ixbpr(ibx+1)-1
        iax=ibpr(ixb)
c 2. beta  k -> l
        itmp=ibfrm(korb,jbx)
        kbx=ibto(lorb,itmp)
        if(kbx.ne.0)then
          phase=phbto(korb,itmp)*phbto(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1+v1(iax,kbx)*term
          res2=res2+v2(iax,kbx)*term
        endif
c 2. alpha k -> l
        itmp=iafrm(korb,iax)
        kax=iato(lorb,itmp)
        if(kax.ne.0)then
          phase=phato(korb,itmp)*phato(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1+v1(kax,jbx)*term
          res2=res2+v2(kax,jbx)*term
        endif
4200    continue
      endif
4100  continue
      endif
      elseif(commut.and.sc)then
c  1) Alpha excitation
      do 1103 ia=1,n1a
      iaxtmp=i1alf(ia,iorb)
      jax=iato(jorb,iaxtmp)
      if(jax.ne.0)then
        iax=iato(iorb,iaxtmp)
        tcof=phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
c  I -> J
        ibx=ndb-iax+1
c 2. alpha k -> l
        itmp=iafrm(korb,jax)
        kax=iato(lorb,itmp)
        if(kax.ne.0)then
          phase=phato(korb,itmp)*phato(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1+v1(kax,ibx)*term
          res2=res2+v2(kax,ibx)*term
        endif
c 2. alpha l -> k
        itmp=iafrm(lorb,iax)
        lax=iato(korb,itmp)
        if(lax.ne.0)then
          phase=phato(lorb,itmp)*phato(korb,itmp)
          term=tcof*phase*cfrom(lax,ibx)
          res1=res1-v1(jax,ibx)*term
          res2=res2-v2(jax,ibx)*term
        endif
c 2. beta  k -> l
        itmp=ibfrm(korb,ibx)
        kbx=ibto(lorb,itmp)
        if(kbx.ne.0)then
          phase=phbto(korb,itmp)*phbto(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1+v1(jax,kbx)*term
          res2=res2+v2(jax,kbx)*term
        endif
c 2. beta  l -> k
        itmp=ibfrm(lorb,ibx)
        lbx=ibto(korb,itmp)
        if(lbx.ne.0)then
          phase=phbto(lorb,itmp)*phbto(korb,itmp)
          term=tcof*phase*cfrom(iax,lbx)
          res1=res1-v1(jax,ibx)*term
          res2=res2-v2(jax,ibx)*term
        endif
c  J -> I
        ibx=ndb-jax+1
c 2. alpha k -> l
        itmp=iafrm(korb,jax)
        kax=iato(lorb,itmp)
        if(kax.ne.0)then
          phase=phato(korb,itmp)*phato(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1-v1(kax,ibx)*term
          res2=res2-v2(kax,ibx)*term
        endif
c 2. alpha l -> k
        itmp=iafrm(lorb,iax)
        lax=iato(korb,itmp)
        if(lax.ne.0)then
          phase=phato(lorb,itmp)*phato(korb,itmp)
          term=tcof*phase*cfrom(lax,ibx)
          res1=res1+v1(jax,ibx)*term
          res2=res2+v2(jax,ibx)*term
        endif
c 2. beta  k -> l
        itmp=ibfrm(korb,ibx)
        kbx=ibto(lorb,itmp)
        if(kbx.ne.0)then
          phase=phbto(korb,itmp)*phbto(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1-v1(jax,kbx)*term
          res2=res2-v2(jax,kbx)*term
        endif
c 2. beta  l -> k
        itmp=ibfrm(lorb,ibx)
        lbx=ibto(korb,itmp)
        if(lbx.ne.0)then
          phase=phbto(lorb,itmp)*phbto(korb,itmp)
          term=tcof*phase*cfrom(iax,lbx)
          res1=res1+v1(jax,ibx)*term
          res2=res2+v2(jax,ibx)*term
        endif
      endif
1103  continue

      if(absym)then
        res1=two*res1
        res2=two*res2
      else
c  2) Beta excitation
      do 2103 ib=1,n1b
      ibxtmp=i1bet(ib,iorb)
      jbx=ibto(jorb,ibxtmp)
      if(jbx.ne.0)then
        ibx=ibto(iorb,ibxtmp)
        tcof=phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
c  I -> J
        iax=nda-ibx+1
c 2. beta  k -> l
        itmp=ibfrm(korb,jbx)
        kbx=ibto(lorb,itmp)
        if(kbx.ne.0)then
          phase=phbto(korb,itmp)*phbto(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1+v1(iax,kbx)*term
          res2=res2+v2(iax,kbx)*term
        endif
c 2. beta  l -> k
        itmp=ibfrm(lorb,ibx)
        lbx=ibto(korb,itmp)
        if(lbx.ne.0)then
          phase=phbto(lorb,itmp)*phbto(korb,itmp)
          term=tcof*phase*cfrom(iax,lbx)
          res1=res1-v1(iax,jbx)*term
          res2=res2-v2(iax,jbx)*term
        endif
c 2. alpha k -> l
        itmp=iafrm(korb,iax)
        kax=iato(lorb,itmp)
        if(kax.ne.0)then
          phase=phato(korb,itmp)*phato(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1+v1(kax,jbx)*term
          res2=res2+v2(kax,jbx)*term
        endif
c 2. alpha l -> k
        itmp=iafrm(lorb,iax)
        lax=iato(korb,itmp)
        if(lax.ne.0)then
          phase=phato(lorb,itmp)*phato(korb,itmp)
          term=tcof*phase*cfrom(lax,ibx)
          res1=res1-v1(iax,jbx)*term
          res2=res2-v2(iax,jbx)*term
        endif
c  J -> I
        iax=nda-jbx+1
c 2. beta  k -> l
        itmp=ibfrm(korb,jbx)
        kbx=ibto(lorb,itmp)
        if(kbx.ne.0)then
          phase=phbto(korb,itmp)*phbto(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1-v1(iax,kbx)*term
          res2=res2-v2(iax,kbx)*term
        endif
c 2. beta  l -> k
        itmp=ibfrm(lorb,ibx)
        lbx=ibto(korb,itmp)
        if(lbx.ne.0)then
          phase=phbto(lorb,itmp)*phbto(korb,itmp)
          term=tcof*phase*cfrom(iax,lbx)
          res1=res1+v1(iax,jbx)*term
          res2=res2+v2(iax,jbx)*term
        endif
c 2. alpha k -> l
        itmp=iafrm(korb,iax)
        kax=iato(lorb,itmp)
        if(kax.ne.0)then
          phase=phato(korb,itmp)*phato(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1-v1(kax,jbx)*term
          res2=res2-v2(kax,jbx)*term
        endif
c 2. alpha l -> k
        itmp=iafrm(lorb,iax)
        lax=iato(korb,itmp)
        if(lax.ne.0)then
          phase=phato(lorb,itmp)*phato(korb,itmp)
          term=tcof*phase*cfrom(lax,ibx)
          res1=res1+v1(iax,jbx)*term
          res2=res2+v2(iax,jbx)*term
        endif
      endif
2103  continue
      endif
      elseif((.not.commut).and.sc)then
c  1) Alpha excitation
      do 3103 ia=1,n1a
      iaxtmp=i1alf(ia,iorb)
      jax=iato(jorb,iaxtmp)
      if(jax.ne.0)then
        iax=iato(iorb,iaxtmp)
        tcof=phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
c  I -> J
        ibx=ndb-iax+1
c 2. alpha k -> l
        itmp=iafrm(korb,jax)
        kax=iato(lorb,itmp)
        if(kax.ne.0)then
          phase=phato(korb,itmp)*phato(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1+v1(kax,ibx)*term
          res2=res2+v2(kax,ibx)*term
        endif
c 2. beta  k -> l
        itmp=ibfrm(korb,ibx)
        kbx=ibto(lorb,itmp)
        if(kbx.ne.0)then
          phase=phbto(korb,itmp)*phbto(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1+v1(jax,kbx)*term
          res2=res2+v2(jax,kbx)*term
        endif
      endif
3103  continue

      if(absym)then
        res1=two*res1
        res2=two*res2
      else
c  2) Beta excitation
      do 4103 ib=1,n1b
      ibxtmp=i1bet(ib,iorb)
      jbx=ibto(jorb,ibxtmp)
      if(jbx.ne.0)then
        ibx=ibto(iorb,ibxtmp)
        tcof=phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
c  I -> J
        iax=nda-ibx+1
c 2. beta  k -> l
        itmp=ibfrm(korb,jbx)
        kbx=ibto(lorb,itmp)
        if(kbx.ne.0)then
          phase=phbto(korb,itmp)*phbto(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1+v1(iax,kbx)*term
          res2=res2+v2(iax,kbx)*term
        endif
c 2. alpha k -> l
        itmp=iafrm(korb,iax)
        kax=iato(lorb,itmp)
        if(kax.ne.0)then
          phase=phato(korb,itmp)*phato(lorb,itmp)
          term=tcof*phase*cfrom(iax,ibx)
          res1=res1+v1(kax,jbx)*term
          res2=res2+v2(kax,jbx)*term
        endif
      endif
4103  continue
      endif
      endif
      hessorb(iprm1,iprm2)=hessorb(iprm1,iprm2)+oaa2*res1+aa1*res2
      hessorb(iprm2,iprm1)=hessorb(iprm1,iprm2)
9379  continue
      return
      end
