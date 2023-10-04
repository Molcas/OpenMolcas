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

!*********************************************************************
!*                                                                   *
!*  PROJ      := Project CASSCF vector onto irrep(s).                *
!*                                                                   *
!*********************************************************************
subroutine proj_cvb(civec)

use casvb_global, only: civbvec
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: civec(*)
#include "main_cvb.fh"
integer(kind=iwp) :: icivec
real(kind=wp) :: dum(mxirrep)

if (projsym) then
  icivec = nint(civec(1))
  call psym1_cvb(civbvec(:,icivec),civbvec(:,icivec),dum,1)
end if

return

end subroutine proj_cvb
