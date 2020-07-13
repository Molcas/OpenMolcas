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
      subroutine oneexc2_cvb(cfrom,cto,vij,
     > i1alf,i1bet,iato,ibto,phato,phbto,
     > iapr,ixapr,ibpr,ixbpr,npvb,
     > nda,ndb,n1a,n1b,nam1,nbm1,norb,commut,sc,absym,diag,idens,
     > iPvb)
c  Calculates Cto = Pvb Eij Cfrom
      implicit real*8 (a-h,o-z)
      logical commut,sc,absym,diag
      dimension cfrom(nda,ndb),cto(nda,ndb)
      dimension vij(*)
      dimension i1alf(n1a,norb),i1bet(n1b,norb)
      dimension iato(norb,0:nam1),ibto(norb,0:nbm1)
      dimension phato(norb,nam1),phbto(norb,nbm1)
      dimension iapr(npvb),ixapr(nda+1),ibpr(npvb),ixbpr(ndb+1)
      save two,half,thresh
      data two/2d0/,half/.5d0/,thresh/1.d-10/

      if(diag)then
        nvij=norb*norb
      elseif(idens.eq.1)then
        nvij=norb*(norb-1)
      endif
      if(absym)then
        if(idens.eq.0)then
          call dscal_(nda*ndb,half,cto,1)
        else
          call dscal_(nvij,half,vij,1)
        endif
      endif
      iprm=0
      do 1100 jorb=1,norb
      do 1200 iorb=1,norb
      if(iorb.eq.jorb.and..not.diag)goto 1200
      iprm=iprm+1
      if(idens.eq.0.and.abs(vij(iprm)).lt.thresh)goto 1200
      if(.not.sc)then
c  a) Alpha excitation
      do 2100 ia=1,n1a
      iaxtmp=i1alf(ia,iorb)
      jax=iato(jorb,iaxtmp)
      if(jax.ne.0)then
        iax=iato(iorb,iaxtmp)
        if(idens.eq.0)then
          tcof=vij(iprm)*phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
          if(iPvb.eq.0)then
            call daxpy_(ndb,tcof,cfrom(jax,1),nda,cto(iax,1),nda)
          elseif(iPvb.eq.1)then
            do 2200 ixa=ixapr(jax),ixapr(jax+1)-1
            ibx=iapr(ixa)
            cto(iax,ibx)=cto(iax,ibx)+tcof*cfrom(jax,ibx)
2200        continue
          elseif(iPvb.eq.2)then
            do 2300 ixa=ixapr(iax),ixapr(iax+1)-1
            ibx=iapr(ixa)
            cto(iax,ibx)=cto(iax,ibx)+tcof*cfrom(jax,ibx)
2300        continue
          endif
        else
          tcof=phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
          if(iPvb.eq.0)then
            vij(iprm)=vij(iprm)+tcof*
     >        ddot_(ndb,cto(iax,1),nda,cfrom(jax,1),nda)
          elseif(iPvb.eq.1)then
            do 2400 ixa=ixapr(jax),ixapr(jax+1)-1
            ibx=iapr(ixa)
            vij(iprm)=vij(iprm)+tcof*cto(iax,ibx)*cfrom(jax,ibx)
2400        continue
          elseif(iPvb.eq.2)then
            do 2500 ixa=ixapr(iax),ixapr(iax+1)-1
            ibx=iapr(ixa)
            vij(iprm)=vij(iprm)+tcof*cto(iax,ibx)*cfrom(jax,ibx)
2500        continue
          endif
        endif
      endif
2100  continue

      if(.not.absym)then
c  c) Beta excitation
        do 3100 ib=1,n1b
        ibxtmp=i1bet(ib,iorb)
        jbx=ibto(jorb,ibxtmp)
        if(jbx.ne.0)then
          ibx=ibto(iorb,ibxtmp)
          if(idens.eq.0)then
            tcof=vij(iprm)*phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
            if(iPvb.eq.0)then
              call daxpy_(nda,tcof,cfrom(1,jbx),1,cto(1,ibx),1)
            elseif(iPvb.eq.1)then
              do 3200 ixb=ixbpr(jbx),ixbpr(jbx+1)-1
              iax=ibpr(ixb)
              cto(iax,ibx)=cto(iax,ibx)+tcof*cfrom(iax,jbx)
3200          continue
            elseif(iPvb.eq.2)then
              do 3300 ixb=ixbpr(ibx),ixbpr(ibx+1)-1
              iax=ibpr(ixb)
              cto(iax,ibx)=cto(iax,ibx)+tcof*cfrom(iax,jbx)
3300          continue
            endif
          else
            tcof=phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
            if(iPvb.eq.0)then
              vij(iprm)=vij(iprm)+tcof*
     >          ddot_(nda,cto(1,ibx),1,cfrom(1,jbx),1)
            elseif(iPvb.eq.1)then
              do 3400 ixb=ixbpr(jbx),ixbpr(jbx+1)-1
              iax=ibpr(ixb)
              vij(iprm)=vij(iprm)+tcof*cto(iax,ibx)*cfrom(iax,jbx)
3400          continue
            elseif(iPvb.eq.2)then
              do 3500 ixb=ixbpr(ibx),ixbpr(ibx+1)-1
              iax=ibpr(ixb)
              vij(iprm)=vij(iprm)+tcof*cto(iax,ibx)*cfrom(iax,jbx)
3500          continue
            endif
          endif
        endif
3100    continue
      endif
      elseif(sc)then
c  a) Alpha excitation
      do 4100 ia=1,n1a
      iaxtmp=i1alf(ia,iorb)
      jax=iato(jorb,iaxtmp)
      if(jax.ne.0)then
        iax=iato(iorb,iaxtmp)
        if(idens.eq.0)then
          tcof=vij(iprm)*phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
          if(iPvb.eq.0)then
            call daxpy_(ndb,tcof,cfrom(jax,1),nda,cto(iax,1),nda)
          elseif(iPvb.eq.1)then
            ibx=ndb-jax+1
            cto(iax,ibx)=cto(iax,ibx)+tcof*cfrom(jax,ibx)
          elseif(iPvb.eq.2)then
            ibx=ndb-iax+1
            cto(iax,ibx)=cto(iax,ibx)+tcof*cfrom(jax,ibx)
          endif
        else
          tcof=phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
          if(iPvb.eq.0)then
            vij(iprm)=vij(iprm)+tcof*
     >        ddot_(ndb,cto(iax,1),nda,cfrom(jax,1),nda)
          elseif(iPvb.eq.1)then
            ibx=ndb-jax+1
            vij(iprm)=vij(iprm)+tcof*cto(iax,ibx)*cfrom(jax,ibx)
          elseif(iPvb.eq.2)then
            ibx=ndb-iax+1
            vij(iprm)=vij(iprm)+tcof*cto(iax,ibx)*cfrom(jax,ibx)
          endif
        endif
      endif
4100  continue

      if(.not.absym)then
c  c) Beta excitation
        do 5100 ib=1,n1b
        ibxtmp=i1bet(ib,iorb)
        jbx=ibto(jorb,ibxtmp)
        if(jbx.ne.0)then
          ibx=ibto(iorb,ibxtmp)
          if(idens.eq.0)then
            tcof=vij(iprm)*phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
            if(iPvb.eq.0)then
              call daxpy_(nda,tcof,cfrom(1,jbx),1,cto(1,ibx),1)
            elseif(iPvb.eq.1)then
              iax=nda-jbx+1
              cto(iax,ibx)=cto(iax,ibx)+tcof*cfrom(iax,jbx)
            elseif(iPvb.eq.2)then
              iax=nda-ibx+1
              cto(iax,ibx)=cto(iax,ibx)+tcof*cfrom(iax,jbx)
            endif
          else
            tcof=phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
            if(iPvb.eq.0)then
              vij(iprm)=vij(iprm)+tcof*
     >          ddot_(nda,cto(1,ibx),1,cfrom(1,jbx),1)
            elseif(iPvb.eq.1)then
              iax=nda-jbx+1
              vij(iprm)=vij(iprm)+tcof*cto(iax,ibx)*cfrom(iax,jbx)
            elseif(iPvb.eq.2)then
              iax=nda-ibx+1
              vij(iprm)=vij(iprm)+tcof*cto(iax,ibx)*cfrom(iax,jbx)
            endif
          endif
        endif
5100    continue
      endif
      endif
1200  continue
1100  continue
      if(absym)then
        if(idens.eq.0)then
          do 6000 ia=1,nda
          do 6001 ib=ia,nda
          cto(ia,ib)=cto(ia,ib)+cto(ib,ia)
          cto(ib,ia)=cto(ia,ib)
6001      continue
6000      continue
        else
          call dscal_(nvij,two,vij,1)
        endif
      endif
      return
c Avoid unused argument warnings
      if (.false.) call Unused_logical(commut)
      end
