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
function igetcnt2_cvb(ident_ci)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: igetcnt2_cvb
integer(kind=iwp) :: ident_ci
#include "main_cvb.fh"

igetcnt2_cvb = icnt_ci(ident_ci)

return

end function igetcnt2_cvb
