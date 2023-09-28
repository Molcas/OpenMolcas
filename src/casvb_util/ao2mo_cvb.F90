!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine ao2mo_cvb(orbsao,orbs,norb1)

use casvb_global, only: nbas_mo
use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
integer(kind=iwp) :: norb1
real(kind=wp) :: orbsao(nbas_mo,norb1), orbs(norb,norb1)
#include "WrkSpc.fh"
integer(kind=iwp) :: i1
integer(kind=iwp), external :: mstackr_cvb

if (norb1 == 0) return
i1 = mstackr_cvb(nbas_mo*norb)
call getmo_cvb(work(i1),3)
call mxattb_cvb(work(i1),orbsao,norb,nbas_mo,norb1,orbs)
call mfreer_cvb(i1)

return

end subroutine ao2mo_cvb
