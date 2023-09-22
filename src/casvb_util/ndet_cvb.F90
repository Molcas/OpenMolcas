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
function ndet_cvb(nel1,nalf1)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ndet_cvb
integer(kind=iwp) :: nel1, nalf1
integer(kind=iwp) :: nd

call icomb_cvb(nel1,nalf1,nd)
ndet_cvb = nd

return

end function ndet_cvb
