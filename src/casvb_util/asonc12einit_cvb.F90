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
subroutine asonc12einit_cvb(ippinp)

use casvb_global, only: ipp12e, iter12e
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ippinp
#include "main_cvb.fh"

iter12e = 0
ipp12e = ippinp

return

end subroutine asonc12einit_cvb
