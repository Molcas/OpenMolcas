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

subroutine ORDINT(IINST,IOUTST,NELMNT,INO,IPRNT)
! ORDER A STRING OF INTEGERS TO ASCENDING ORDER
!
! IINST : INPUT STRING
! IOUTST : OUTPUT STRING
! NELMNT : NUMBER OF INTEGERS
! INO : Mapping array from new to old order
!
! THIS CODE CONTAINS THE OLD ORDER CODE OF JOE GOLAB
! (HE IS HEREBY ACKNOWLEDGED, AND I AM EXCUSED)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: NELMNT, IINST(NELMNT), IOUTST(NELMNT), INO(NELMNT), IPRNT
integer(kind=iwp) :: I, ISWAP, JOE, NTEST

if (NELMNT == 0) goto 1001
call ICOPVE(IINST,IOUTST,NELMNT)
do I=1,NELMNT
  INO(I) = I
end do

! BEGIN TO ORDER

JOE = 1
10 I = JOE
20 continue
if (I == NELMNT) GO TO 50
if (IOUTST(I) <= IOUTST(I+1)) GO TO 40
JOE = I+1
30 ISWAP = IOUTST(I)
IOUTST(I) = IOUTST(I+1)
IOUTST(I+1) = ISWAP
ISWAP = INO(I)
INO(I) = INO(I+1)
INO(I+1) = ISWAP
if (I == 1) GO TO 10
I = I-1
if (IOUTST(I) > IOUTST(I+1)) GO TO 30
GO TO 10
40 I = I+1
GO TO 20

! END ORDER

50 continue

1001 continue
NTEST = 0
NTEST = max(NTEST,IPRNT)
if (NTEST >= 200) then
  write(u6,*) ' Result from ORDINT'
  write(u6,*)
  write(u6,*) ' Input string'
  call IWRTMA(IINST,1,NELMNT,1,NELMNT)
  write(u6,*) ' Ordered string'
  call IWRTMA(IOUTST,1,NELMNT,1,NELMNT)
  write(u6,*) ' New to old order'
  call IWRTMA(INO,1,NELMNT,1,NELMNT)
end if

end subroutine ORDINT
