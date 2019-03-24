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
      subroutine gaussj_cvb(orbs,igjorb)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension orbs(norb,norb),igjorb(*)
      dimension ogjorb(norb,norb)
c  *********************************************************************
c  *                                                                   *
c  *  GAUSSJ    := Define sequence of simple updates from orb transf.  *
c  *                                                                   *
c  *********************************************************************

      k1 = mstackr_cvb(norb*norb)
      k2 = mstacki_cvb(norb)
      k3 = mstacki_cvb(norb)
      k4 = mstacki_cvb(norb)
      call fmove_cvb(orbs,w(k1),norb*norb)
      ioff=idbl_cvb(norb*norb)
!     IFG: Note that originally ogjorb was igjorb.
!          There's probably a bug somewhere, but for now we just
!          make the compiler happy with the transfer statement
      call gaussj2_cvb(w(k1),iw(k2),iw(k3),iw(k4),
     >  igjorb(1+ioff),igjorb(1+norb+ioff),ogjorb,norb)
      igjorb(1:norb*norb)=transfer(ogjorb,igjorb(1),norb*norb)
      call imove_cvb(igjorb(1+ioff),iw(k2),norb)
      do 100 i=1,norb
100   igjorb(iw(i+k2-1)+ioff)=i
      call mfreer_cvb(k1)
      return
      end
