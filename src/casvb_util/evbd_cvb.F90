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

subroutine evbd_cvb(orbs,cvb,fx,ioptc,iter)

use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: orbs(norb,norb), cvb(nvb), fx
integer(kind=iwp) :: ioptc, iter
#include "optze_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i1, i2, i3, i4, i5, i6, i7, mstackr_cvb1
integer(kind=iwp), external :: mavailr_cvb, mstackr_cvb

mstackr_cvb1 = mavailr_cvb()
maxdav = min((mstackr_cvb1-ndet-nvb)/(3*nvb),nvb)
i1 = mstackr_cvb(nvb*maxdav)
i2 = mstackr_cvb(nvb*maxdav)
i3 = mstackr_cvb(nvb*maxdav)
i4 = mstackr_cvb(nvb)
i5 = mstackr_cvb(maxdav*maxdav)
i6 = mstackr_cvb(maxdav)
i7 = mstackr_cvb(maxdav)
call evbd2_cvb(orbs,cvb,fx,ioptc,iter,work(lw(4)),work(lw(5)),work(lw(6)),work(i1),work(i2),work(i3),work(i4),work(i5),work(i6), &
               work(i7))
call mfreer_cvb(i1)

return

end subroutine evbd_cvb
