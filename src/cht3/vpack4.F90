!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SUBROUTINE VPACK4(IN,IND,NBI,NINT)
      implicit none
      integer IN,IND,NBI,NINT, I
!R8   1
!     INTEGER*2 IND
      DIMENSION IN(*),IND(*)
!REA  2
!     DATA I1,I2,I3 /48,32,16/
!      NBV=4*NBI
!R8   2
      DO 2 I=1,NINT
      IND(I)=IN(I)
 2    CONTINUE
!REA  4
!     IF(NINT.NE.NBV)THEN
!     DO 1 I=NINT+1,NBV
!1     IN(I)=0
!      ENDIF
!CDC 11
!     CALL Q8SHIFTV (X'08',,IN(1;NBI),,I1,,IND(1;NBI))
!      NN=NBI+1
!     CALL Q8LINKV  (X'10')
!     CALL Q8SHIFTV (X'08',,IN(NN;NBI),,I2,,IN(1;NBI))
!     CALL Q8ORV    (X'02',,IND(1;NBI),,IN(1;NBI),,IND(1;NBI))
!      NN=NBI+NN
!     CALL Q8LINKV  (X'10')
!     CALL Q8SHIFTV (X'08',,IN(NN;NBI),,I3,,IN(1;NBI))
!     CALL Q8ORV    (X'02',,IND(1;NBI),,IN(1;NBI),,IND(1;NBI))
!      NN=NBI+NN
!      CALL Q8ORV    (X'02',,IND(1;NBI),,IN(NN;NBI),,IND(1;NBI))
!REA 11
!      DO 2 I=1,NBI
!2     IND(I)=SHIFT(IN(I),I1)
!      NN=NBI
!      DO 3 I=1,NBI
!3     IND(I)=OR(IND(I),SHIFT(IN(I+NN),I2))
!      NN=NBI*2
!      DO 4 I=1,NBI
!4     IND(I)=OR(IND(I),SHIFT(IN(I+NN),I3))
!      NN=NBI*3
!      DO 5 I=1,NBI
!5     IND(I)=OR(IND(I),      IN(I+NN))
      RETURN
! Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(NBI)
      END
