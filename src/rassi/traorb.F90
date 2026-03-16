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

subroutine TRAORB(NSYM,NOSH,NBASF,NCXA,CXA,NCMO,CMO)
! TRANSFORM ORBITAL COEFFICIENTS CMO BY MULTIPLYING WITH
! TRANSFORMATION MATRIX CXA.

use definitions, only: iwp, wp
use constants, only: Zero, One
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer(kind=iwp), intent(in) :: NSYM, NCXA, NCMO
integer(kind=iwp), intent(in) :: NOSH(NSYM), NBASF(NSYM)
real(kind=wp), intent(in) :: CXA(NCXA)
real(kind=wp), intent(inout) :: CMO(NCMO)
real(kind=wp), allocatable :: CNew(:)
integer(kind=iwp) ISTA1, ISTA2, IS, NO, NB

call mma_allocate(CNEW,NCMO,Label='CNEW')
ISTA1 = 1
ISTA2 = 1
do IS=1,NSYM
  NO = NOSH(IS)
  if (NO == 0) cycle
  NB = NBASF(IS)
  if (NB /= 0) then
    call DGEMM_('N','N',NB,NO,NO,One,CMO(ISTA1),NB,CXA(ISTA2),NO,Zero,CNEW(ISTA1),NB)
    ISTA1 = ISTA1+NO*NB
  end if
  ISTA2 = ISTA2+NO**2
end do
call DCOPY_(NCMO,CNEW,1,CMO,1)
call mma_deallocate(CNEW)

end subroutine TRAORB
