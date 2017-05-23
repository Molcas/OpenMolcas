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
      SUBROUTINE VPACK4(IN,IND,NBI,NINT)
      implicit none
      integer IN,IND,NBI,NINT, I
CR8   1
C     INTEGER*2 IND
      DIMENSION IN(*),IND(*)
CREA  2
C     DATA I1,I2,I3 /48,32,16/
C      NBV=4*NBI
CR8   2
      DO 2 I=1,NINT
 2    IND(I)=IN(I)
CREA  4
C     IF(NINT.NE.NBV)THEN
C     DO 1 I=NINT+1,NBV
C1     IN(I)=0
C      ENDIF
CCDC 11
C     CALL Q8SHIFTV (X'08',,IN(1;NBI),,I1,,IND(1;NBI))
C      NN=NBI+1
C     CALL Q8LINKV  (X'10')
C     CALL Q8SHIFTV (X'08',,IN(NN;NBI),,I2,,IN(1;NBI))
C     CALL Q8ORV    (X'02',,IND(1;NBI),,IN(1;NBI),,IND(1;NBI))
C      NN=NBI+NN
C     CALL Q8LINKV  (X'10')
C     CALL Q8SHIFTV (X'08',,IN(NN;NBI),,I3,,IN(1;NBI))
C     CALL Q8ORV    (X'02',,IND(1;NBI),,IN(1;NBI),,IND(1;NBI))
C      NN=NBI+NN
C      CALL Q8ORV    (X'02',,IND(1;NBI),,IN(NN;NBI),,IND(1;NBI))
CREA 11
C      DO 2 I=1,NBI
C2     IND(I)=SHIFT(IN(I),I1)
C      NN=NBI
C      DO 3 I=1,NBI
C3     IND(I)=OR(IND(I),SHIFT(IN(I+NN),I2))
C      NN=NBI*2
C      DO 4 I=1,NBI
C4     IND(I)=OR(IND(I),SHIFT(IN(I+NN),I3))
C      NN=NBI*3
C      DO 5 I=1,NBI
C5     IND(I)=OR(IND(I),      IN(I+NN))
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(NBI)
      END
