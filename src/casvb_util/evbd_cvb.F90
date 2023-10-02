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

!IFG trivial
subroutine evbd_cvb(orbs,cvb,fx,ioptc,iter)

use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: orbs(norb,norb), cvb(nvb), fx
integer(kind=iwp) :: ioptc, iter
#include "optze_cvb.fh"
#include "WrkSpc.fh"

call evbd2_cvb(orbs,cvb,fx,ioptc,iter,work(lw(4)),work(lw(5)),work(lw(6)))

return

end subroutine evbd_cvb
