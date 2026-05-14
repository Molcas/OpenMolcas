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

!#define _DEBUGPRINT_
subroutine ORDINT(IINST,IOUTST,NELMNT,INO)
! ORDER A STRING OF INTEGERS TO ASCENDING ORDER
!
! IINST : INPUT STRING
! IOUTST : OUTPUT STRING
! NELMNT : NUMBER OF INTEGERS
! INO : Mapping array from new to old order
!
! THIS CODE CONTAINS THE OLD ORDER CODE OF JOE GOLAB
! (HE IS HEREBY ACKNOWLEDGED, AND I AM EXCUSED)

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NELMNT, IINST(NELMNT)
integer(kind=iwp), intent(out) :: IOUTST(NELMNT), INO(NELMNT)
integer(kind=iwp) :: I, ISWAP, JOE

if (NELMNT /= 0) then
  IOUTST(:) = IINST(:)
  INO(:) = [(I,I=1,NELMNT)]

  ! BEGIN TO ORDER

  JOE = 1
  I = JOE
  do
    if (I == NELMNT) exit
    if (IOUTST(I) <= IOUTST(I+1)) then
      I = I+1
      cycle
    end if
    JOE = I+1
    do
      ISWAP = IOUTST(I)
      IOUTST(I) = IOUTST(I+1)
      IOUTST(I+1) = ISWAP
      ISWAP = INO(I)
      INO(I) = INO(I+1)
      INO(I+1) = ISWAP
      if (I == 1) exit
      I = I-1
      if (IOUTST(I) <= IOUTST(I+1)) exit
    end do
    I = JOE
  end do

  ! END ORDER
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' Result from ORDINT'
write(u6,*)
write(u6,*) ' Input string'
call IWRTMA(IINST,1,NELMNT,1,NELMNT)
write(u6,*) ' Ordered string'
call IWRTMA(IOUTST,1,NELMNT,1,NELMNT)
write(u6,*) ' New to old order'
call IWRTMA(INO,1,NELMNT,1,NELMNT)
#endif

end subroutine ORDINT
