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

subroutine genprexyz6(preY,preXZ)

use Definitions, only: wp

implicit none
#include "para.fh"
real(kind=wp) :: preY(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax), preXZ(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)

!bs ####################################################################
!bs   prefactors preXZ und preY include the factors 1/root(2)
!bs   for the +/- linear combinations of spherical harmonics
!bs ####################################################################
preY(:,:,:,:) = preXZ(:,:,:,:)

return

end subroutine genprexyz6
