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
!*************************************************************
!** Routines for imposing constraints on VB wfn. parameters **
!*************************************************************
!*********************
!** Set-up routines **
!*********************

subroutine construc_cvb(tconstr,ipermzeta)

use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: tconstr(nvb,nvb)
integer(kind=iwp) :: ipermzeta(norb,nzeta)
#include "WrkSpc.fh"
integer(kind=iwp) :: i1, i2, i3
integer(kind=iwp), external :: mstackr_cvb

i1 = mstackr_cvb(norb*norb)
i2 = mstackr_cvb(norb*norb)
i3 = mstackr_cvb(norb*norb)
call setipermzeta_cvb(ipermzeta,work(lv(1)),work(ls(1)),iwork(ls(13)),work(i1),work(i2),work(i3))
call mfreer_cvb(i1)
if (iconstruc == 2) call construc2_cvb(tconstr)

return

end subroutine construc_cvb
