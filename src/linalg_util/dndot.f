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
      SUBROUTINE DNDOT(N,M,S,INCS,ISW,X,INCXI,INCXO,Y,INCYI,INCYO)
C
C     COMPUTE DOT PRODUCT N TIMES
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 S(1+(N-1)*INCS),
     1       X((1+(M-1)*INCXI)*(1+(N-1)*INCXO)),
     2       Y((1+(M-1)*INCYI)*(1+(N-1)*INCYO))
      IF (ISW.EQ.1) THEN
         DO 10 I = 1, N
            S(1+(I-1)*INCS) = DDOT_(M,X(1+(I-1)*INCXO),INCXI,
     1                               Y(1+(I-1)*INCYO),INCYI)
   10    CONTINUE
      ELSE IF (ISW.EQ.2) THEN
         DO 20 I = 1, N
            S(1+(I-1)*INCS) =-DDOT_(M,X(1+(I-1)*INCXO),INCXI,
     1                               Y(1+(I-1)*INCYO),INCYI)
   20    CONTINUE
      ELSE IF (ISW.EQ.3) THEN
         DO 30 I = 1, N
            S(1+(I-1)*INCS) = S(1+(I-1)*INCS) +
     1                        DDOT_(M,X(1+(I-1)*INCXO),INCXI,
     2                               Y(1+(I-1)*INCYO),INCYI)
   30    CONTINUE
      ELSE IF (ISW.EQ.4) THEN
         DO 40 I = 1, N
            S(1+(I-1)*INCS) = S(1+(I-1)*INCS) -
     1                        DDOT_(M,X(1+(I-1)*INCXO),INCXI,
     2                               Y(1+(I-1)*INCYO),INCYI)
   40    CONTINUE
      ELSE
         Call SysAbendMsg('dndot','ISW IS OUT OF RANGE IN DNDOT',' ')
      END IF
      RETURN
      END
