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
!    G01,G02    : Real*8 - (G02-G01) = energy difference between
!                 the two states.
!    x_anharm1  : Real*8 two dimensional array - anharmonicity
!                 constants for ground state.
!    x_anharm2  : Real*8 two dimensional array - anharmonicity
!                 constants for excited state.
!    harmfreq1  : Real*8 array - harmonic frequencies for
!                 ground state.
!    harmfreq2  : Real*8 array - harmonic frequencies for
!                 excited state.
!    level1     : Integer array - quanta for ground state.
!    level2     : Integer array - quanta for excited state.
!
!  Output:
!    energy     : Real*8 variable - energy for transition
!                 between level1 and level2.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use Constants, only: Half

implicit real*8(a-h,o-z)
#include "Constants_mula.fh"
real*8 sum, G0, G1, G01, G02
real*8 energy
real*8 x_anharm1(nDim,nDim), x_anharm2(nDim,nDim)
real*8 harmfreq1(nDim), harmfreq2(nDim)
integer level1(nDim), level2(nDim)

! Calculate energy for level 1.
sum = G01
do i=1,nDim
  sum = sum+harmfreq1(i)*(level1(i)+Half)
  do j=i,nDim
    sum = sum+x_anharm1(i,j)*(level1(i)+Half)*(level1(j)+Half)
  end do
end do
G0 = sum

! Calculate energy for level 2.
sum = G02
do i=1,nDim
  sum = sum+harmfreq2(i)*(level2(i)+Half)
  do j=i,nDim
    sum = sum+x_anharm2(i,j)*(level2(i)+Half)*(level2(j)+Half)
  end do
end do
G1 = sum

! Calculate energy difference.
energy = G1-G0

end subroutine TransEnergy
