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
      subroutine change1_cvb()
      implicit real*8 (a-h,o-z)
      logical changed
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "frag_cvb.fh"
#include "rls_cvb.fh"

c  Arrays for determinant handling (E_ij) and definition
c  of VB wavefunction :
      changed=.false.
      if(chpcmp_cvb(norb))changed=.true.
      if(chpcmp_cvb(nalf))changed=.true.
      if(chpcmp_cvb(nbet))changed=.true.
      if(chpcmp_cvb(nel))changed=.true.

      if(changed)call touch_cvb('CASPRINT')

      if(chpcmp_cvb(nconf))changed=.true.

      if(.not.changed)call cnfchk_cvb()
      nvb=nvb_cvb(kbasis)

      if(chpcmp_cvb(ndetvb))changed=.true.
      if(chpcmp_cvb(mxion))changed=.true.
      if(chpcmp_cvb(mnion))changed=.true.
      if(changed)call touch_cvb('MEM1')
      return

      entry chop1_cvb()
      if(release(1))call mfreei_cvb(ll(1))
      release(1)=.true.
      release(2)=.false.

c  Dimensions
      call icomb_cvb(norb,nalf,nda)
      call icomb_cvb(norb,nbet,ndb)
      do i=1,nfrag
      call icomb_cvb(norb,nalf_fr(i,1),nda_fr(i,1))
      call icomb_cvb(norb,nbet_fr(i,1),ndb_fr(i,1))
      enddo
      call icomb_cvb(norb-1,nalf-1,n1a)
      call icomb_cvb(norb-1,nbet-1,n1b)
      call icomb_cvb(norb,nalf-1,nam1)
      call icomb_cvb(norb,nbet-1,nbm1)
      ndet = nda*ndb
c  Symmetry of determinant strings :
      call getnci_cvb(ncivb,nel,nalf-nbet,0)
c  Identical indexing arrays may share memory :
      ll(1) = mstacki_cvb(norb*n1a)
      ll(2) = ll(1)
      if(.not.absym(4))ll(2) = mstacki_cvb(norb*n1b)
      ll(3) = mstacki_cvb(norb*nda)
      ll(4) = ll(3)
      if(.not.absym(4))ll(4) = mstacki_cvb(norb*ndb)
      ll(5) = mstacki_cvb(norb*(nam1+1))
      ll(6) = ll(5)
      if(.not.absym(4))ll(6) = mstacki_cvb(norb*(nbm1+1))

c  7 & 8 (PHAFRM & PHBFRM) taken out

      ll(9) = mstackr_cvb(norb*nam1)
      ll(10)= ll(9)
      if(.not.absym(4))ll(10)= mstackr_cvb(norb*nbm1)

c  Determinant dimensioning for VB wavefunction :
      ndavb=0
      ndbvb=0
      naprodvb=1
      nbprodvb=1
      nvbprod=1
      npvb=1
      do ifrag=1,nfrag
      ndavb=ndavb+nda_fr(1,ifrag)+1
      ndbvb=ndbvb+ndb_fr(1,ifrag)+1
      naprodvb=naprodvb*nda_fr(1,ifrag)
      nbprodvb=nbprodvb*ndb_fr(1,ifrag)
      nvbprod=nvbprod*ndetvb_fr(ifrag)
      npvb=npvb*ndetvb_fr(ifrag)
      enddo
      if(nfrag.le.1)then
        naprodvb=0
        nbprodvb=0
        nvbprod=0
      endif

      ll(11)= mstacki_cvb(npvb)
      ll(12)= mstacki_cvb(nda+1)
      ll(13)= mstacki_cvb(npvb)
      ll(14)= mstacki_cvb(ndb+1)
      ll(15)= mstacki_cvb(nconf*noe)
c  16 obsolete (former ioncty)
      ll(17)= mstacki_cvb(ndetvb)
c  Use 7 & 8 for IA12IND / IB12IND
      ll(7)= mstacki_cvb(naprodvb)
      ll(8)= mstacki_cvb(nbprodvb)

      ll(20)= mstacki_cvb(ndetvb)
      ll(21)= mstacki_cvb(ndavb)
      ll(22)= mstacki_cvb(ndetvb)
      ll(23)= mstacki_cvb(ndbvb)
      return
      end
