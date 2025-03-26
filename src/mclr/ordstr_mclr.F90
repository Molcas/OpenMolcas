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
subroutine ORDSTR_MCLR(IINST,IOUTST,NELMNT,ISIGN)
! ORDER A STRING OF INTEGERS TO ASCENDING ORDER
!
! IINST : INPUT STRING IS IINST
! IOUTST : OUTPUT STRING IS IOUTST
! NELMNT : NUMBER OF INTEGERS IN STRING
! ISIGN :  SIGN OF PERMUTATION : + 1 : EVEN PERMUTATIONN
!                                - 1 : ODD  PERMUTATION
!
! THIS CODE CONTAINS THE OLD ORDER CODE OF JOE GOLAB
! (HE IS HEREBY AKNOWLEDGED, AND I AM EXCUSED)
!
! IMPLEMENTED MORE TRANSPARENT BUBBLE SORTING INSTEAD
!               JR NOV 2006

implicit none
integer NELMNT
integer IINST(NELMNT), IOUTST(NELMNT)
integer ISIGN

integer iTemp, iPass, I

if (NELMNT == 0) return

ISIGN = 1
iTEMP = 0

10 continue
iPass = 0
do I=1,NELMNT-1
  if (IINST(I) > IINST(I+1)) then
    iTEMP = IINST(I)
    IINST(I) = IINST(I+1)
    IINST(I+1) = iTEMP
    ISIGN = -1*ISIGN
    iPass = 1
  end if
end do
if (IPASS /= 0) goto 10

do I=1,NELMNT
  IOUTST(I) = IINST(I)
end do

!call iCOPY(NELMNT,IINST,1,IOUTST,1)
!ISIGN = 1

!JOE = 1
!10 continue
!I = JOE
!20 continue
!if (I == NELMNT) goto 50
!if (IOUTST(I) <= IOUTST(I+1)) goto 40
!JOE = I+1
!30 continue
!iSWAP = IOUTST(I)
!ISIGN = -ISIGN
!IOUTST(I) = IOUTST(I+1)
!IOUTST(I+1) = iSWAP
!if (I == 1) goto 10
!I = I-1
!if (IOUTST(I) > IOUTST(I+1)) goto 30
!goto 10
!40 continue
!I = I+1
!goto 20
!
! END ORDER
!
!50 continue
!#ifdef _DEBUGPRINT_
!write(6,*) ' INPUT STRING ORDERED STRING ISIGN ',NELMNT
!call IWRTMA(IINST,1,NELMNT,1,NELMNT)
!call IWRTMA(IOUTST,1,NELMNT,1,NELMNT)
!write(6,*) ' ISIGN : ',ISIGN
!#endif

end subroutine ORDSTR_MCLR
