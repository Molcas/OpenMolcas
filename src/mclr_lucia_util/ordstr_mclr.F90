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
! Copyright (C) 1984,1989-1993, Jeppe Olsen                            *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine ORDSTR_MCLR(IINST,IOUTST,NELMNT,ISGN)
! ORDER A STRING OF INTEGERS TO ASCENDING ORDER
!
! IINST  : INPUT STRING IS IINST
! IOUTST : OUTPUT STRING IS IOUTST
! NELMNT : NUMBER OF INTEGERS IN STRING
! ISGN   : SIGN OF PERMUTATION : + 1 : EVEN PERMUTATIONN
!                                - 1 : ODD  PERMUTATION
!
! THIS CODE CONTAINS THE OLD ORDER CODE OF JOE GOLAB
! (HE IS HEREBY AKNOWLEDGED, AND I AM EXCUSED)
!
! IMPLEMENTED MORE TRANSPARENT BUBBLE SORTING INSTEAD
!               JR NOV 2006

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: NELMNT
integer(kind=iwp), intent(inout) :: IINST(NELMNT)
integer(kind=iwp), intent(out) :: IOUTST(NELMNT), ISGN
integer(kind=iwp) :: I, iPass, iTemp

ISGN = 1
if (NELMNT == 0) return

iTEMP = 0

iPass = 1
do while (iPass /= 0)
  iPass = 0
  do I=1,NELMNT-1
    if (IINST(I) > IINST(I+1)) then
      iTEMP = IINST(I)
      IINST(I) = IINST(I+1)
      IINST(I+1) = iTEMP
      ISGN = -1*ISGN
      iPass = 1
    end if
  end do
end do

IOUTST(:) = IINST(:)

!ISGN = 1

!JOE = 1
!I = JOE
!do
!  if (I == NELMNT) exit
!  if (IOUTST(I) > IOUTST(I+1)) then
!    JOE = I+1
!    do
!      iSWAP = IOUTST(I)
!      ISGN = -ISGN
!      IOUTST(I) = IOUTST(I+1)
!      IOUTST(I+1) = iSWAP
!      if (I == 1) exit
!      I = I-1
!      if (IOUTST(I) <= IOUTST(I+1)) exit
!    end do
!    I = JOE
!  else
!    I = I+1
!  end if
!end do
!
! END ORDER
!
!#ifdef _DEBUGPRINT_
!write(u6,*) ' INPUT STRING ORDERED STRING ISGN ',NELMNT
!call IWRTMA(IINST,1,NELMNT,1,NELMNT)
!call IWRTMA(IOUTST,1,NELMNT,1,NELMNT)
!write(u6,*) ' ISGN : ',ISGN
!#endif

end subroutine ORDSTR_MCLR
