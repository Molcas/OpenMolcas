!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine genprexyz4(preXZ)

use Constants, only: Two
use Definitions, only: wp

implicit none
#include "para.fh"
real(kind=wp) :: preXZ(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)
real(kind=wp), parameter :: roottwo = sqrt(Two)

!bs ####################################################################
!bs   prefactors preXZ und preY include the factors 1/root(2)
!bs   for the +/- linear combinations of spherical harmonics
!bs ####################################################################
preXZ(:,0,:,:) = preXZ(:,0,:,:)*roottwo

return

end subroutine genprexyz4
