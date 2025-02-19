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
! Copyright (C) 2001, Jeppe Olsen                                      *
!***********************************************************************

subroutine CMP_IVEC_ILIST(IVEC,ILIST,LLIST,NLIST,INUM)
! An integer IVEC of LLIST entries are given.
! compare with list of vectors in ILIST and find first
! vector in LLIST that is identical to IVEC.
!
! If INUM = 0, the list was not found
!
! Jeppe Olsen, December 2001

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LLIST, IVEC(LLIST), NLIST, ILIST(LLIST,NLIST)
integer(kind=iwp), intent(out) :: INUM
integer(kind=iwp) :: IELMNT, JLIST, NTEST

INUM = 0
do JLIST=1,NLIST
  do IELMNT=1,LLIST
    if (IVEC(IELMNT) /= ILIST(IELMNT,JLIST)) exit
  end do
  if (IELMNT > LLIST) then
    INUM = JLIST
    exit
  end if
end do

NTEST = 0
if (NTEST >= 100) then
  write(u6,*) ' Input list :'
  call IWRTMA(IVEC,1,LLIST,1,LLIST)
  write(u6,*) ' Address of list : ',INUM
end if

end subroutine CMP_IVEC_ILIST
