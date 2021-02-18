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

#include "malloc_cvb.fh"
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
        call permvb2_cvb(v1,iperm,vb,iw(ll(11)),iw(ll(12)),
     >    iw(k1),iw(k2),iw(k5),iw(k6),
     >    iw(k7),iw(k8),iw(k9),iw(k10),iw(k11),iw(k12),
     >    iw(k13),w(k14),iw(k15),w(k16),w(k17),ialg)
      else
        iv1=nint(v1(1))
cDLC iw(maddr_r2i_cvb(iaddr_ci(iv1))) --> w(iaddr_ci(iv1))
        call permvb2_cvb(w(iaddr_ci(iv1)),
     >    iperm,vb,iw(ll(11)),iw(ll(12)),
     >    iw(k1),iw(k2),iw(k5),iw(k6),
     >    iw(k7),iw(k8),iw(k9),iw(k10),iw(k11),iw(k12),
     >    iw(k13),w(k14),iw(k15),w(k16),w(k17),ialg)
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

#include "malloc_cvb.fh"
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
        call permvb2_cvb(v1,iperm,vb,iw(ll(11)),iw(ll(12)),
     >    iw(k1),iw(k2),iw(k5),iw(k6),
     >    iw(k7),iw(k8),iw(k9),iw(k10),iw(k11),iw(k12),
     >    iw(k13),w(k14),iw(k15),w(k16),w(k17),ialg)
      else
        iv1=nint(v1(1))
cDLC iw(maddr_r2i_cvb(iaddr_ci(iv1))) --> w(iaddr_ci(iv1))
        call permvb2_cvb(w(iaddr_ci(iv1)),
     >    iperm,vb,iw(ll(11)),iw(ll(12)),
     >    iw(k1),iw(k2),iw(k5),iw(k6),
     >    iw(k7),iw(k8),iw(k9),iw(k10),iw(k11),iw(k12),
     >    iw(k13),w(k14),iw(k15),w(k16),w(k17),ialg)
      endif
      call mfreei_cvb(k1)
      return
      end
