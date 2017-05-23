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
      SUBROUTINE DEX2_CP2 (CPQ,NDWN,NUPA,A,NUPB,B,NCP,ICOUP,VTAB)
      IMPLICIT NONE
      INTEGER NDWN,NUPA,NUPB,NCP,ICOUP(3,NCP)
      REAL*8 CPQ,A(NUPA,NDWN),B(NUPB,NDWN),VTAB(*)
      INTEGER ICP,ILFT,IRGT
      REAL*8 X
C CASE: ADD EPQ*A TO B, WHERE MIDLEV<P<Q AND A AND B ARE SINGLE
C MATRIX BLOCKS.
C ONLY UPPER WALKS AFFECTED =>    ROW OPERATION.
C P<Q, DEEXCITING OPERATOR  => A IS RIGHTHAND WALK.
C CHOSE DAXPY IF LARGE LENGTH.
      DO ICP=1,NCP
        ILFT=ICOUP(1,ICP)
        IRGT=ICOUP(2,ICP)
        X=CPQ*VTAB(ICOUP(3,ICP))
        CALL DAXPY_(NDWN,X,A(IRGT,1),NUPA,B(ILFT,1),NUPB)
      END DO
      RETURN
      END
