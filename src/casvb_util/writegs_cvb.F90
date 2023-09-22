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

subroutine writegs_cvb()

use Definitions, only: iwp

implicit none
#include "main_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i1
integer(kind=iwp), external :: mstacki_cvb

i1 = mstacki_cvb(ndetvb)
call writegs2_cvb(work(lv(1)),work(lv(2)),work(lw(9)),iwork(ll(11)),iwork(ll(12)),iwork(i1))
call mfreei_cvb(i1)

return

end subroutine writegs_cvb
