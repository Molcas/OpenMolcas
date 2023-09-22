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
subroutine o7b_cvb(nparm1,dxnrm,grdnrm,close2conv)

use casvb_global, only: ix
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nparm1
real(kind=wp) :: dxnrm, grdnrm
logical(kind=iwp) :: close2conv
#include "WrkSpc.fh"

call o7b2_cvb(nparm1,work(ix(1)),dxnrm,grdnrm,close2conv)

return

end subroutine o7b_cvb
