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
      subroutine igaussj_cvb(orbs,igjorb)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension orbs(norb,norb),igjorb(*)
c  *********************************************************************
c  *                                                                   *
c  *  GAUSSJ    := Define sequence of simple updates from orb transf.  *
c  *                                                                   *
c  *********************************************************************

      call igaussj_cvb_internal(igjorb)
*
*     This is to allow type punning without an explicit interface
      contains
      subroutine igaussj_cvb_internal(igjorb)
      use iso_c_binding
      integer, target :: igjorb(*)
      real*8, pointer :: gjorb(:)
      k1 = mstackr_cvb(norb*norb)
      k2 = mstacki_cvb(norb)
      k3 = mstacki_cvb(norb)
      k4 = mstacki_cvb(norb)
      call fmove_cvb(orbs,w(k1),norb*norb)
      ioff=idbl_cvb(norb*norb)
      call c_f_pointer(c_loc(igjorb(1)),gjorb,[norb*norb])
      call gaussj2_cvb(w(k1),iw(k2),iw(k3),iw(k4),
     >  igjorb(1+ioff),igjorb(1+norb+ioff),gjorb,norb)
      nullify(gjorb)
      call imove_cvb(igjorb(1+ioff),iw(k2),norb)
      do 100 i=1,norb
      igjorb(iw(i+k2-1)+ioff)=i
100   continue
      call mfreer_cvb(k1)
      return
      end subroutine igaussj_cvb_internal
*
      end
*
      subroutine gaussj_cvb(orbs,gjorb)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
      dimension orbs(norb,norb),gjorb(*)
      call gaussj_cvb_internal(gjorb)
*
*     This is to allow type punning without an explicit interface
      contains
      subroutine gaussj_cvb_internal(gjorb)
      use iso_c_binding
      real*8, target :: gjorb(*)
      integer, pointer :: igjorb(:)
      call c_f_pointer(c_loc(gjorb(1)),igjorb,[1])
      call igaussj_cvb(orbs,igjorb)
      nullify(igjorb)
      end subroutine gaussj_cvb_internal
*
      end
