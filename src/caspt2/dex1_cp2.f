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
      SUBROUTINE DEX1_CP2 (CPQ,NUP,A,B,NCP,ICOUP,VTAB)
      IMPLICIT NONE
      INTEGER NUP,NCP,ICOUP(3,NCP)
      REAL*8 CPQ,A(NUP,*),B(NUP,*),VTAB(*)
      INTEGER ICP,JLFT,JRGT
      REAL*8 X
C CASE: ADD EPQ*A TO B, WHERE P<Q<=MIDLEV AND A AND B ARE SINGLE
C MATRIX BLOCKS.
C ONLY LOWER WALKS AFFECTED => COLUMN OPERATION.
C P<Q, DEEXCITING OPERATOR  => A IS RIGHTHAND WALK.
C CHOSE DAXPY IF LARGE VECTOR LENGTH:
      DO ICP=1,NCP
        JLFT=ICOUP(1,ICP)
        JRGT=ICOUP(2,ICP)
        X=CPQ*VTAB(ICOUP(3,ICP))
        CALL DAXPY_(NUP,X,A(1,JRGT),1,B(1,JLFT),1)
      END DO
      RETURN
      END
