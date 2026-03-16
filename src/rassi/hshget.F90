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

subroutine HSHGET(KEY,KEYDIM,NCOMP,ITEM,NSIZE,ITAB,ITEMID)

use definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: KEYDIM, NCOMP, NSIZE
integer(kind=iwp), intent(in) :: KEY(KEYDIM), ITEM(NCOMP,*)
integer(kind=iwp), intent(in) :: ITAB(NSIZE,2)
integer(kind=iwp), intent(out) :: ITEMID
! These parameters determine the hash function:
integer(kind=iwp), parameter :: MULT = 37, NHASH = 997
integer(kind=iwp) NULL, IND, I, LOOKAT
logical(kind=iwp) FAILED_ID

if (NSIZE < NHASH) then
  write(u6,*) ' HSHGET: Table size must be at least as'
  write(u6,*) '         big as NHASH, presently =',NHASH
  call ABEND()
end if
NULL = ITAB(NSIZE,1)
! Evaluate hash index:
IND = mod(KEY(1),NHASH)
do I=2,KEYDIM
  IND = mod(KEY(I)+MULT*IND,NHASH)
end do
IND = IND+1
! IND is a hashed index in interval 1..NHASH < NSIZE
! Find the item with this key:

LOOKAT = IND

Outer: do

  ! Are there (more) items with that hash signature?
  if (ITAB(LOOKAT,1) == NULL) then
    ! Here, if we have failed to find such an item.
    ITEMID = 0
    return
  end if

  ! Try to identify an item which has the given key:
  ITEMID = ITAB(LOOKAT,2)
  FAILED_ID = .false.
  do I=1,KEYDIM
    if (ITEM(I,ITEMID) /= KEY(I)) then
      FAILED_ID = .true.
      exit
    end if
  end do

  if (FAILED_ID) then
    ! Here, if we have not yet identified the item.
    LOOKAT = ITAB(LOOKAT,1)
  else
    ! Here, if we have identified the item.
    exit
  end if

end do Outer

end subroutine HSHGET
