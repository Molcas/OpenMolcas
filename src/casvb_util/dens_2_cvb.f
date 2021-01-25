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
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
c  *********************************************************************
c  *                                                                   *
c  *  DENS   := calculate two-electron (transition) densities          *
c  *                                                                   *
c  *********************************************************************
      subroutine dens_2_cvb(v1,v2,cfrom,d2mata,d2matb,
     > i1alf,i1bet,iafrm,ibfrm,iato,ibto,phato,phbto,
     > d1mata,d1matb,
     > nda,ndb,n1a,n1b,nam1,nbm1,norb,commut,sc,absym)
c  Calculates V1 EijEkl CFROM and V2 EijEkl CFROM
      implicit real*8 (a-h,o-z)
      logical commut,sc,absym
      dimension v1(nda,ndb),v2(nda,ndb),cfrom(nda,ndb)
      dimension d2mata(norb*norb,norb*norb),d2matb(norb*norb,norb*norb)
      dimension i1alf(n1a,norb),i1bet(n1b,norb)
      dimension iafrm(norb,nda),ibfrm(norb,ndb)
      dimension iato(norb,0:nam1),ibto(norb,0:nbm1)
      dimension phato(norb,nam1),phbto(norb,nbm1)
      dimension d1mata(norb,norb),d1matb(norb,norb)
      dimension res1(20*20),res2(20*20)
      save thresh,two
      data thresh/1d-10/,two/2d0/

      call dev1b2_cvb(v1,cfrom,d1mata,
     >  i1alf,i1bet,iato,ibto,phato,phbto,
     >  norb*norb,nda,ndb,n1a,n1b,nam1,nbm1,norb,commut,sc,.true.,absym)
      call dev1b2_cvb(v1,v2,d1matb,
     >  i1alf,i1bet,iato,ibto,phato,phbto,
     >  norb*norb,nda,ndb,n1a,n1b,nam1,nbm1,norb,commut,sc,.true.,absym)

      do 9379 ip1=1,norb*norb
      iorb=(ip1-1)/norb+1
      jorb=ip1-(iorb-1)*norb
      call fzero(res1,ip1)
      call fzero(res2,ip1)

c  1) Alpha excitation
      do 3100 ia=1,n1a
      iaxtmp=i1alf(ia,iorb)
      jax=iato(jorb,iaxtmp)
      if(jax.ne.0)then
        iax=iato(iorb,iaxtmp)
        tcof=phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
        do 93791 ip2=1,ip1
        korb=(ip2-1)/norb+1
        lorb=ip2-(korb-1)*norb
c 2. alpha k -> l
        itmp=iafrm(korb,jax)
        kax=iato(lorb,itmp)
        if(kax.ne.0)then
c  I -> J
          do 3200 ibx=1,ndb
          if(abs(cfrom(iax,ibx)).gt.thresh)then
            phase=phato(korb,itmp)*phato(lorb,itmp)
            term=tcof*phase*cfrom(iax,ibx)
            res1(ip2)=res1(ip2)+v1(kax,ibx)*term
            res2(ip2)=res2(ip2)+v2(kax,ibx)*term
          endif
3200      continue
        endif
93791   continue
c 2. beta  k -> l
        do 3201 ibx=1,ndb
        if(abs(cfrom(iax,ibx)).gt.thresh)then
          do 93792 ip2=1,ip1
          korb=(ip2-1)/norb+1
          lorb=ip2-(korb-1)*norb
          itmp=ibfrm(korb,ibx)
          kbx=ibto(lorb,itmp)
          if(kbx.ne.0)then
            phase=phbto(korb,itmp)*phbto(lorb,itmp)
            term=tcof*phase*cfrom(iax,ibx)
            res1(ip2)=res1(ip2)+v1(jax,kbx)*term
            res2(ip2)=res2(ip2)+v2(jax,kbx)*term
          endif
93792     continue
        endif
3201    continue
      endif
3100  continue

      if(.not.absym)then
c  2) Beta excitation
      do 4100 ib=1,n1b
      ibxtmp=i1bet(ib,iorb)
      jbx=ibto(jorb,ibxtmp)
      if(jbx.ne.0)then
        ibx=ibto(iorb,ibxtmp)
        tcof=phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
        do 93793 ip2=1,ip1
        korb=(ip2-1)/norb+1
        lorb=ip2-(korb-1)*norb
c 2. beta  k -> l
        itmp=ibfrm(korb,jbx)
        kbx=ibto(lorb,itmp)
        if(kbx.ne.0)then
c  I -> J
          do 4200 iax=1,nda
          if(abs(cfrom(iax,ibx)).gt.thresh)then
            phase=phbto(korb,itmp)*phbto(lorb,itmp)
            term=tcof*phase*cfrom(iax,ibx)
            res1(ip2)=res1(ip2)+v1(iax,kbx)*term
            res2(ip2)=res2(ip2)+v2(iax,kbx)*term
          endif
4200      continue
        endif
93793   continue
c 2. alpha k -> l
        do 4201 iax=1,nda
        if(abs(cfrom(iax,ibx)).gt.thresh)then
          do 93794 ip2=1,ip1
          korb=(ip2-1)/norb+1
          lorb=ip2-(korb-1)*norb
          itmp=iafrm(korb,iax)
          kax=iato(lorb,itmp)
          if(kax.ne.0)then
            phase=phato(korb,itmp)*phato(lorb,itmp)
            term=tcof*phase*cfrom(iax,ibx)
            res1(ip2)=res1(ip2)+v1(kax,jbx)*term
            res2(ip2)=res2(ip2)+v2(kax,jbx)*term
          endif
93794     continue
        endif
4201    continue
      endif
4100  continue
      else
        do 93795 ip2=1,ip1
        res1(ip2)=two*res1(ip2)
        res2(ip2)=two*res2(ip2)
93795   continue
      endif

      do 93796 ip2=1,ip1
      korb=(ip2-1)/norb+1
      lorb=ip2-(korb-1)*norb
      d2mata(ip1,ip2)=res1(ip2)
      d2matb(ip1,ip2)=res2(ip2)
c  D(lk|ji) = E_lk E_ji - \delta_kj E_li
      if(korb.eq.jorb)then
        d2mata(ip1,ip2)=d2mata(ip1,ip2)-d1mata(lorb,iorb)
        d2matb(ip1,ip2)=d2matb(ip1,ip2)-d1matb(lorb,iorb)
      endif
      d2mata(ip2,ip1)=d2mata(ip1,ip2)
      d2matb(ip2,ip1)=d2matb(ip1,ip2)
93796 continue
9379  continue
      return
      end
