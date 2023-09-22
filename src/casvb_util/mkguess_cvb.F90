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

subroutine mkguess_cvb()

use casvb_global, only: nbas_mo
use Definitions, only: iwp

implicit none
#include "main_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iorbsao, irdorbs
integer(kind=iwp), external :: mstacki_cvb, mstackr_cvb

irdorbs = mstacki_cvb(norb)
iorbsao = mstackr_cvb(nbas_mo*norb)
call mkguess2_cvb(work(lv(1)),work(lv(2)),iwork(irdorbs),work(iorbsao))
call mfreei_cvb(irdorbs)

return

end subroutine mkguess_cvb
