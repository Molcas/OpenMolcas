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
      SUBROUTINE DEX1 (CPQ,NUP,A,B,NCP,ICOUP,VTAB)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER   ICOUP(3,NCP)
      DIMENSION A(NUP,*),B(NUP,*),VTAB(*)
C CASE: ADD EPQ*A TO B, WHERE P<Q<=MIDLEV AND A AND B ARE SINGLE
C MATRIX BLOCKS.
C ONLY LOWER WALKS AFFECTED => COLUMN OPERATION.
C P<Q, DEEXCITING OPERATOR  => A IS RIGHTHAND WALK.
C CHOSE DAXPY IF LARGE VECTOR LENGTH:
      IF(NUP.GT.20) THEN
        DO ICP=1,NCP
          JLFT=ICOUP(1,ICP)
          JRGT=ICOUP(2,ICP)
          X=CPQ*VTAB(ICOUP(3,ICP))
          CALL DAXPY_(NUP,X,A(1,JRGT),1,B(1,JLFT),1)
        END DO
      ELSE
        DO ICP=1,NCP
          JLFT=ICOUP(1,ICP)
          JRGT=ICOUP(2,ICP)
          X=CPQ*VTAB(ICOUP(3,ICP))
          DO I=1,NUP
            B(I,JLFT)=B(I,JLFT)+X*A(I,JRGT)
          END DO
        END DO
      END IF
      RETURN
      END
