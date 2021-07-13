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
! Copyright (C) 1987, Jeppe Olsen                                      *
!               1989, Markus P. Fuelscher                              *
!***********************************************************************
      SUBROUTINE ORDSTR(IINST,IOUTST,NELMNT,ISIGN,IPRINT)
!
!     AUTHOR:        J. OLSEN, UNIV. OF LUND, SWEDEN, APRIL 1987
!     MODIFICATIONS: INCLUSION INTO THE RASSCF METHOD
!                    M.P. FUELSCHER, UNIV. OF LUND, SWEDEN, MAY 1989
!
!     PURPOSE:
!
!     ORDER A STRING OF INTEGERS TO ASCENDING ORDER
!     IINST : INPUT STRING IS IINST
!     IOUTST : OUTPUT STRING IS IOUTST
!     NELMNT : NUMBER OF INTEGERS IN STRING
!     ISIGN :  SIGN OF PERMUTATION : + 1 : EVEN PERMUTATIONN
!                                    - 1 : ODD  PERMUTATION
!
!     SUBROUTINE CALLS:
!
!     ICOPY,IWRTMA
!
!
!     THIS CODE CONTAINS THE OLD ORDER CODE OF JOE GOLAB
!     ( HE IS HEREBY AKNOWLEDGED , AND I AM EXCUSED )
!
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IINST(NELMNT),IOUTST(NELMNT)
      INTEGER SWAP
!
      IF (NELMNT.EQ.0) RETURN
!
      CALL ICOPY(NELMNT,IINST,1,IOUTST,1)
      ISIGN = 1
!
!       BEGIN TO ORDER
!
        JOE = 1
10      I = JOE
20      CONTINUE
        IF(I.EQ.NELMNT) GO TO 50
        IF(IOUTST(I).LE.IOUTST(I+1)) GO TO 40
        JOE = I + 1
30      SWAP = IOUTST(I)
        ISIGN = - ISIGN
        IOUTST(I) = IOUTST(I+1)
        IOUTST(I+1) = SWAP
        IF(I.EQ.1) GO TO 10
        I = I - 1
        IF(IOUTST(I).GT.IOUTST(I+1)) GO TO 30
        GO TO 10
40      I = I + 1
      GO TO 20
!
!     END ORDER
!
50    CONTINUE
      IF( IPRINT.GT.30 ) THEN
        Write(6,*)  ' INPUT STRING ORDERED STRING ISIGN '
        CALL IWRTMA(IINST,1,NELMNT,1,NELMNT)
        CALL IWRTMA(IOUTST,1,NELMNT,1,NELMNT)
        Write(6,*) ' ISIGN : ', ISIGN
      END IF
!
      RETURN
      END
