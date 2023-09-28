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

subroutine asonc12s_cvb( &
#                       define _CALLING_
#                       include "ddasonc_interface.fh"
                       )
! Applies S and H on c vector(s).

use Definitions, only: wp, iwp

implicit none
#include "ddasonc_interface.fh"
#include "main_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i1
integer(kind=iwp), external :: mstackr_cvb

#include "macros.fh"
unused_var(axc)

i1 = mstackr_cvb(nvb+nprorb)
call asonc12s2_cvb(c,sxc,nvec,nprm,work(lc(3)),work(lc(2)),work(lv(1)),work(lw(4)),work(lw(5)),work(lw(6)),work(lw(9)), &
                   work(lv(2)),work(i1))
call mfreer_cvb(i1)

return

end subroutine asonc12s_cvb
