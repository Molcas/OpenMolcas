************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE SQUARE(A, B,ICB,IRB, NROW)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(*)

* PAM Sep 06: The two special cases here account for
* almost all calls of this code, and written such as
* not to store with non-zero stride--at least some
* rudimentary consideration for cache.
      if(ICB.EQ.1) GOTO 100
      if(IRB.EQ.1) GOTO 200
* General and inefficient code:
      IND=0
      DO 10 IROW=0,NROW-1
      DO 10 ICOL=0,IROW
         IND=IND+1
         B(1+IROW*ICB+ICOL*IRB)=A(IND)
         B(1+ICOL*ICB+IROW*IRB)=A(IND)
10    CONTINUE
      GOTO 900

 100  CONTINUE
      DO IC=0,NROW-1
       DO IR=0,IC
        B(1+IR+IC*IRB)=A(1+IR+(IC*(IC+1))/2)
       END DO
      END DO
      DO IC=0,NROW-2
       DO IR=IC+1,NROW-1
        B(1+IR+IC*IRB)=B(1+IC+IR*IRB)
       END DO
      END DO
      GOTO 900

 200  CONTINUE
      DO IC=0,NROW-1
       DO IR=0,IC
        B(1+IR+IC*ICB)=A(1+IR+(IC*(IC+1))/2)
       END DO
      END DO
      DO IC=0,NROW-2
       DO IR=IC+1,NROW-1
        B(1+IR+IC*ICB)=B(1+IC+IR*ICB)
       END DO
      END DO

 900  CONTINUE
      RETURN

      END
