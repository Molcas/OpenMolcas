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
!    x_anharm   : Real two dimensional array - anharmonicity constants.
!    harmfreq   : Real array - harmonic frequencies.
!
!  Output:
!    anharmfreq : Real array - anharmonic frequencies.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nOsc
real(kind=wp), intent(in) :: x_anharm(nOsc,nOsc), harmfreq(nOsc)
real(kind=wp), intent(out) :: anharmfreq(nOsc)
integer(kind=iwp) :: istate
real(kind=wp) :: G1, G2
integer(kind=iwp), allocatable :: level1(:), level2(:)

call mma_allocate(level1,nOsc,label='level1')
call mma_allocate(level2,nOsc,label='level2')

G1 = Zero
G2 = G1
level1(:) = 0
do istate=1,nOsc
  level2(:) = 0
  level2(istate) = 1
  call TransEnergy(G1,x_anharm,harmfreq,level1,G2,x_anharm,harmfreq,level2,anharmfreq(istate),nOsc)
end do

call mma_deallocate(level1)
call mma_deallocate(level2)

end subroutine AnharmonicFreq
