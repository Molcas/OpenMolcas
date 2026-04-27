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

subroutine HSHPUT(KEYDIM,NCOMP,ITEM,NSIZE,ITAB,ITEMID)

use rassi_aux, only: MULT, NHASH
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: KEYDIM, NCOMP, ITEM(NCOMP,*), NSIZE, ITEMID
integer(kind=iwp), intent(inout) :: ITAB(NSIZE,2)
integer(kind=iwp) :: I, IFREE, IND, INULL, LOOKAT, NEXT

if (NSIZE < NHASH) then
  write(u6,*) ' HSHPUT: Table size must be at least as'
  write(u6,*) '         big as NHASH, presently =',NHASH
  call ABEND()
end if
INULL = ITAB(NSIZE,1)
IFREE = ITAB(NSIZE,2)
if (ITAB(IFREE,1) == INULL) then
  write(u6,*) ' HSHPUT: Table is already full.'
  write(u6,*) ' Size NSIZE is too small, NSIZE =',NSIZE
  call ABEND()
end if
! Evaluate hash index:
IND = mod(ITEM(1,ITEMID),NHASH)
do I=2,KEYDIM
  IND = mod(ITEM(I,ITEMID)+MULT*IND,NHASH)
end do
IND = IND+1
! IND is a hashed index in interval 1..NHASH < NSIZE
! Find the last item with this key:

LOOKAT = IND
do
  ! Are there already items with that hash signature?
  if (ITAB(LOOKAT,1) == INULL) exit
  LOOKAT = ITAB(LOOKAT,1)
end do

! No more items with the same signature.
! Put the new item in the table at a free location.
ITAB(LOOKAT,1) = IFREE
ITAB(LOOKAT,2) = ITEMID
NEXT = ITAB(IFREE,1)
ITAB(IFREE,1) = INULL
ITAB(NSIZE,2) = NEXT

end subroutine HSHPUT
