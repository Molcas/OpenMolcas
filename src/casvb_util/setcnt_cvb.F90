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
subroutine setcnt_cvb(xident_ci,idep1)

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: xident_ci(1)
integer(kind=iwp) :: idep1
#include "main_cvb.fh"
integer(kind=iwp) :: ident_ci

ident_ci = nint(xident_ci(1))
call setcnt2_cvb(ident_ci,idep1)

return

end subroutine setcnt_cvb
