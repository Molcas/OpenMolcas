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

subroutine AnharmonicFreq(x_anharm,harmfreq,anharmfreq,nOsc)
!  Purpose:
!    Calculate the anharmonic frequencies.
!
!  Input:
!    x_anharm   : Real*8 two dimensional array - anharmonicity
!                 constants.
!    harmfreq   : Real*8 array - harmonic frequencies.
!
!  Output:
!    anharmfreq : Real*8 array - anharmonic frequencies.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero

implicit real*8(a-h,o-z)
#include "Constants_mula.fh"
real*8 G1, G2
real*8 x_anharm(nosc,nosc)
real*8 harmfreq(nosc)
real*8 anharmfreq(nosc)
integer, allocatable :: level1(:), level2(:)

!nDim = nOsc
call mma_allocate(level1,nOsc,label='level1')
call mma_allocate(level2,nOsc,label='level2')

G1 = Zero
G2 = G1
level1(:) = 0
do istate=1,nOsc
  level2(:) = 0
  level2(istate) = 1
  !l_harm = nOsc
  call TransEnergy(G1,x_anharm,harmfreq,level1,G2,x_anharm,harmfreq,level2,anharmfreq(istate),nOsc)
end do

call mma_deallocate(level1)
call mma_deallocate(level2)

end subroutine AnharmonicFreq
