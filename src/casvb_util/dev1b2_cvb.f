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
      subroutine dev1b2_cvb(cfrom,cto,dmat,
     > i1alf,i1bet,iato,ibto,phato,phbto,
     > iparmx,nda,ndb,n1a,n1b,nam1,nbm1,norb,commut,sc,absym,diag)
c  Calculates all CTO Eij CFROM
      implicit real*8 (a-h,o-z)
      logical commut,sc,absym,diag
      dimension cfrom(nda,ndb),cto(nda,ndb),dmat(iparmx)
      dimension i1alf(n1a,norb),i1bet(n1b,norb)
      dimension iato(norb,0:nam1),ibto(norb,0:nbm1)
      dimension phato(norb,nam1),phbto(norb,nbm1)
      save two
      data two/2d0/

      call fzero(dmat,iparmx)

      iparm=0
      do 1100 iorb=1,norb
      do 1200 jorb=1,norb
      if(jorb.eq.iorb.and..not.diag)goto 1200
      iparm=iparm+1
      if(iparm.gt.iparmx)return

c  a) Alpha excitation
      do 2100 ia=1,n1a
      iaxtmp=i1alf(ia,iorb)
      jax=iato(jorb,iaxtmp)
      if(jax.ne.0)then
        iax=iato(iorb,iaxtmp)
        tcof=phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
        dmat(iparm)=dmat(iparm)+tcof*
     >    ddot_(ndb,cto(jax,1),nda,cfrom(iax,1),nda)
      endif
2100  continue

      if(.not.absym)then
c  c) Beta excitation
        do 3100 ib=1,n1b
        ibxtmp=i1bet(ib,iorb)
        jbx=ibto(jorb,ibxtmp)
        if(jbx.ne.0)then
          ibx=ibto(iorb,ibxtmp)
          tcof=phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
          dmat(iparm)=dmat(iparm)+tcof*
     >      ddot_(nda,cto(1,jbx),1,cfrom(1,ibx),1)
        endif
3100    continue
      else
        dmat(iparm)=two*dmat(iparm)
      endif
1200  continue
1100  continue
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_logical(commut)
        call Unused_logical(sc)
      end if
      end
