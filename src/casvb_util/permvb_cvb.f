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
c  *  PERMVB := Permute orbitals in V1 according to IPERM.             *
c  *  PERMCI := Permute orbitals in V1 according to IPERM.             *
c  *                                                                   *
c  *  V1 is either full CI vector (NDET), or just                      *
c  *  VB determinants (NDETVB).                                        *
c  *                                                                   *
c  *********************************************************************
      subroutine permvb_cvb(v1,iperm)
c  Permutes orbitals in V1 according to IPERM.
      implicit real*8 (a-h,o-z)
      logical vb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "WrkSpc.fh"
      dimension iperm(norb),v1(*)

      vb=.true.
      k1 = mstacki_cvb((norb+1)*(nalf+1))
      k2 = mstacki_cvb((norb+1)*(nbet+1))
      k5 = mstacki_cvb(norb+1)
      k6 = mstacki_cvb(norb+1)
      k7 = mstacki_cvb(norb+1)
      k8 = mstacki_cvb(norb)
      k9 = mstacki_cvb(norb)
      k10= mstacki_cvb(norb)
      k11= mstacki_cvb(norb)
      k12= mstacki_cvb(norb)
      k13= mstacki_cvb(nda)
      k14= mstackr_cvb(nda)
      k15= mstacki_cvb(ndb)
      k16= mstackr_cvb(ndb)
      if(vb)then
        k17= mstackr_cvb(ndetvb)
      elseif(mavailr_cvb().ge.ndet)then
        ialg=1
        k17= mstackr_cvb(ndet)
      else
        ialg=2
        k17= mstackr_cvb(nda)
      endif
      if(vb)then
        call permvb2_cvb(v1,iperm,vb,iwork(ll(11)),iwork(ll(12)),
     >    iwork(k1),iwork(k2),iwork(k5),iwork(k6),
     >    iwork(k7),iwork(k8),iwork(k9),iwork(k10),iwork(k11),
     >    iwork(k12),
     >    iwork(k13),work(k14),iwork(k15),work(k16),work(k17),ialg)
      else
        iv1=nint(v1(1))
cDLC iwork(maddr_r2i_cvb(iaddr_ci(iv1))) --> work(iaddr_ci(iv1))
        call permvb2_cvb(work(iaddr_ci(iv1)),
     >    iperm,vb,iwork(ll(11)),iwork(ll(12)),
     >    iwork(k1),iwork(k2),iwork(k5),iwork(k6),
     >    iwork(k7),iwork(k8),iwork(k9),iwork(k10),iwork(k11),
     >    iwork(k12),
     >    iwork(k13),work(k14),iwork(k15),work(k16),work(k17),ialg)
      endif
      call mfreei_cvb(k1)
      return
      end
      subroutine permci_cvb(v1,iperm)
      implicit real*8 (a-h,o-z)
      logical vb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "WrkSpc.fh"
      dimension iperm(norb),v1(*)
      vb=.false.
      k1 = mstacki_cvb((norb+1)*(nalf+1))
      k2 = mstacki_cvb((norb+1)*(nbet+1))
      k5 = mstacki_cvb(norb+1)
      k6 = mstacki_cvb(norb+1)
      k7 = mstacki_cvb(norb+1)
      k8 = mstacki_cvb(norb)
      k9 = mstacki_cvb(norb)
      k10= mstacki_cvb(norb)
      k11= mstacki_cvb(norb)
      k12= mstacki_cvb(norb)
      k13= mstacki_cvb(nda)
      k14= mstackr_cvb(nda)
      k15= mstacki_cvb(ndb)
      k16= mstackr_cvb(ndb)
      if(vb)then
        k17= mstackr_cvb(ndetvb)
      elseif(mavailr_cvb().ge.ndet)then
        ialg=1
        k17= mstackr_cvb(ndet)
      else
        ialg=2
        k17= mstackr_cvb(nda)
      endif
      if(vb)then
        call permvb2_cvb(v1,iperm,vb,iwork(ll(11)),iwork(ll(12)),
     >    iwork(k1),iwork(k2),iwork(k5),iwork(k6),
     >    iwork(k7),iwork(k8),iwork(k9),iwork(k10),iwork(k11),
     >    iwork(k12),
     >    iwork(k13),work(k14),iwork(k15),work(k16),work(k17),ialg)
      else
        iv1=nint(v1(1))
cDLC iwork(maddr_r2i_cvb(iaddr_ci(iv1))) --> work(iaddr_ci(iv1))
        call permvb2_cvb(work(iaddr_ci(iv1)),
     >    iperm,vb,iwork(ll(11)),iwork(ll(12)),
     >    iwork(k1),iwork(k2),iwork(k5),iwork(k6),
     >    iwork(k7),iwork(k8),iwork(k9),iwork(k10),iwork(k11),
     >    iwork(k12),
     >    iwork(k13),work(k14),iwork(k15),work(k16),work(k17),ialg)
      endif
      call mfreei_cvb(k1)
      return
      end
