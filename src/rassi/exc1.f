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
      SUBROUTINE EXC1 (CPQ,NUP,A,B,NCP,ICOUP,VTAB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NUP,*),B(NUP,*),VTAB(*)
      INTEGER   ICOUP(3,NCP)
C CASE: ADD EPQ*A TO B, WHERE Q<P<=MIDLEV AND A AND B ARE SINGLE
C MATRIX BLOCKS.
C ONLY LOWER WALKS AFFECTED => COLUMN OPERATION.
C P>Q, EXCITING OPERATOR  => A IS  LEFTHAND WALK.
C USE VECTOR ROUTINE CALL IF LONG ENOUGH:
      IF(NUP.GT.20) THEN
        DO ICP=1,NCP
          JLFT=ICOUP(1,ICP)
          JRGT=ICOUP(2,ICP)
          X=CPQ*VTAB(ICOUP(3,ICP))
          CALL DAXPY_(NUP,X,A(1,JLFT),1,B(1,JRGT),1)
        END DO
      ELSE
        DO ICP=1,NCP
          JLFT=ICOUP(1,ICP)
          JRGT=ICOUP(2,ICP)
          X=CPQ*VTAB(ICOUP(3,ICP))
          DO I=1,NUP
            B(I,JRGT)=B(I,JRGT)+X*A(I,JLFT)
          END DO
        END DO
      END IF
      RETURN
      END
