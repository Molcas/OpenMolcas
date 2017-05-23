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
      SUBROUTINE EXC2 (CPQ,NDWN,NUPA,A,NUPB,B,NCP,ICOUP,VTAB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NUPA,NDWN),B(NUPB,NDWN),
     *          VTAB(*)
      INTEGER   ICOUP(3,NCP)
C CASE: ADD EPQ*A TO B, WHERE MIDLEV<Q<P AND A AND B ARE SINGLE
C MATRIX BLOCKS.
C ONLY UPPER WALKS AFFECTED =>    ROW OPERATION.
C Q<P, EXCITING OPERATOR  => A IS  LEFTHAND WALK.
      IF(NDWN.GT.15) THEN
        DO ICP=1,NCP
          ILFT=ICOUP(1,ICP)
          IRGT=ICOUP(2,ICP)
          X=CPQ*VTAB(ICOUP(3,ICP))
          CALL DAXPY_(NDWN,X,A(ILFT,1),NUPA,B(IRGT,1),NUPB)
        END DO
      ELSE
        DO ICP=1,NCP
          ILFT=ICOUP(1,ICP)
          IRGT=ICOUP(2,ICP)
          X=CPQ*VTAB(ICOUP(3,ICP))
          DO J=1,NDWN
            B(IRGT,J)=B(IRGT,J)+X*A(ILFT,J)
          END DO
        END DO
      END IF
      RETURN
      END
