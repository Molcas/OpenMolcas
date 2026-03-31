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

subroutine HSHINI(NSIZE,ITAB,INULL)

use rassi_aux, only: NHASH
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NSIZE, INULL
integer(kind=iwp), intent(out) :: ITAB(NSIZE,2)
integer(kind=iwp) :: I, IFREE

if (NSIZE < NHASH) then
  write(u6,*) ' HSHINI: Table size must be at least as'
  write(u6,*) '         big as NHASH, presently =',NHASH
  call ABEND()
end if
do I=1,NHASH
  ITAB(I,1) = INULL
  ITAB(I,2) = INULL
end do
IFREE = NHASH+1
do I=IFREE,NSIZE-1
  ITAB(I,1) = I+1
  ITAB(I,2) = INULL
end do
ITAB(NSIZE,1) = INULL
ITAB(NSIZE,2) = IFREE

end subroutine HSHINI
