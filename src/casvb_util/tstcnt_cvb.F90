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
! Contents of CI vectors:

!IFG trivial
function tstcnt_cvb(civec,idep1)

use casvb_global, only: icnt_ci, ndet
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp) :: tstcnt_cvb
real(kind=wp) :: civec(0:ndet)
integer(kind=iwp) :: idep1
integer(kind=iwp) :: ident_ci

ident_ci = nint(civec(0))
tstcnt_cvb = icnt_ci(ident_ci) == idep1

return

end function tstcnt_cvb
