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
! Copyright (C) Joseph Golab                                           *
!               1987, Jeppe Olsen                                      *
!               1989, Markus P. Fuelscher                              *
!***********************************************************************

subroutine ORDSTR(IINST,IOUTST,NELMNT,ISGN,IPRINT)
! AUTHOR:        J. OLSEN, UNIV. OF LUND, SWEDEN, APRIL 1987
! MODIFICATIONS: INCLUSION INTO THE RASSCF METHOD
!                M.P. FUELSCHER, UNIV. OF LUND, SWEDEN, MAY 1989
!
! PURPOSE:
!
! ORDER A STRING OF INTEGERS TO ASCENDING ORDER
! IINST : INPUT STRING IS IINST
! IOUTST : OUTPUT STRING IS IOUTST
! NELMNT : NUMBER OF INTEGERS IN STRING
! ISGN : SIGN OF PERMUTATION : +1 : EVEN PERMUTATION
!                              -1 : ODD  PERMUTATION
!
! THIS CODE CONTAINS THE OLD ORDER CODE OF JOE GOLAB
! ( HE IS HEREBY AKNOWLEDGED , AND I AM EXCUSED )

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: NELMNT, IINST(NELMNT), IOUTST(NELMNT), ISGN, IPRINT
integer(kind=iwp) :: I, JOE, SWAP

if (NELMNT == 0) return

call ICOPY(NELMNT,IINST,1,IOUTST,1)
ISGN = 1

! BEGIN TO ORDER

JOE = 1
10 continue
I = JOE
20 continue
if (I == NELMNT) GO TO 50
if (IOUTST(I) <= IOUTST(I+1)) GO TO 40
JOE = I+1
30 continue
SWAP = IOUTST(I)
ISGN = -ISGN
IOUTST(I) = IOUTST(I+1)
IOUTST(I+1) = SWAP
if (I == 1) GO TO 10
I = I-1
if (IOUTST(I) > IOUTST(I+1)) GO TO 30
GO TO 10
40 continue
I = I+1
GO TO 20

! END ORDER

50 continue
if (IPRINT > 30) then
  write(6,*) ' INPUT STRING ORDERED STRING ISGN '
  call IWRTMA(IINST,1,NELMNT,1,NELMNT)
  call IWRTMA(IOUTST,1,NELMNT,1,NELMNT)
  write(6,*) ' ISGN : ',ISGN
end if

return

end subroutine ORDSTR
