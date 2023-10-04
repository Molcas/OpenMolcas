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

!***********************************************************************
!*                                                                     *
!*  SYMWEIGHT := CASSCF scalar product divided into irreps.            *
!*                                                                     *
!***********************************************************************
subroutine symweight_cvb(civec1,civec2,osym)

use casvb_global, only: civbvec
use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: civec1(*), civec2(*), osym(mxirrep)
integer(kind=iwp) :: icivec1, icivec2

icivec1 = nint(civec1(1))
icivec2 = nint(civec2(1))
call psym1_cvb(civbvec(:,icivec1),civbvec(:,icivec2),osym,2)

return

end subroutine symweight_cvb
