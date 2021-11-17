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
! Copyright (C) 1996, Niclas Forsberg                                  *
!***********************************************************************

subroutine TransEnergy(G01,x_anharm1,harmfreq1,level1,G02,x_anharm2,harmfreq2,level2,energy,nDim)
!  Purpose:
!    Calculate transition energy between two states.
!
!  Input:
!    G01,G02    : Real - (G02-G01) = energy difference between the two states.
!    x_anharm1  : Real two dimensional array - anharmonicity constants for ground state.
!    x_anharm2  : Real two dimensional array - anharmonicity constants for excited state.
!    harmfreq1  : Real array - harmonic frequencies for ground state.
!    harmfreq2  : Real array - harmonic frequencies for excited state.
!    level1     : Integer array - quanta for ground state.
!    level2     : Integer array - quanta for excited state.
!
!  Output:
!    energy     : Real variable - energy for transition between level1 and level2.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDim, level1(nDim), level2(nDim)
real(kind=wp), intent(in) :: G01, x_anharm1(nDim,nDim), harmfreq1(nDim), G02, x_anharm2(nDim,nDim), harmfreq2(nDim)
real(kind=wp), intent(out) :: energy
integer(kind=iwp) :: i, j
real(kind=wp) :: G0, G1

! Calculate energy for level 1.
G0 = G01
do i=1,nDim
  G0 = G0+harmfreq1(i)*(level1(i)+Half)
  do j=i,nDim
    G0 = G0+x_anharm1(i,j)*(level1(i)+Half)*(level1(j)+Half)
  end do
end do

! Calculate energy for level 2.
G1 = G02
do i=1,nDim
  G1 = G1+harmfreq2(i)*(level2(i)+Half)
  do j=i,nDim
    G1 = G1+x_anharm2(i,j)*(level2(i)+Half)*(level2(j)+Half)
  end do
end do

! Calculate energy difference.
energy = G1-G0

end subroutine TransEnergy
