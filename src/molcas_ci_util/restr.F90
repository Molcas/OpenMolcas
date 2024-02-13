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
      SUBROUTINE RESTR(NVERT0,IDRT0,IDOWN0,IVER, LV1RAS, LV3RAS, LM1RAS, LM3RAS, NVERT)
!     PURPOSE: PUT THE RAS CONSTRAINT TO THE DRT TABLE BY
!              CREATING A MASK
!
      IMPLICIT None
!
      Integer NVERT0, LV1RAS, LV3RAS, LM1RAS, LM3RAS, NVERT
      Integer IDRT0(NVERT0,5),IDOWN0(NVERT0,0:3),IVER(NVERT0)
!
      Integer, PARAMETER :: LTAB=1, NTAB=2
      Integer::  IOR(0:3,0:3)=reshape([ 0,1,2,3,1,1,3,3,2,3,2,3,3,3,3,3 ], [4,4])
      Integer:: IAND(0:3,0:3)=reshape([ 0,0,0,0,0,1,0,1,0,0,2,2,0,1,2,3 ], [4,4])
      Integer:: IV, LEV, N, IVV, IC, ID, MASK, IVD
!
!     LOOP OVER ALL VERTICES AND CHECK ON RAS CONDITIONS
!     CREATE MASK
!
      DO IV=1,NVERT0
        LEV=IDRT0(IV,LTAB)
        N=IDRT0(IV,NTAB)
        IVER(IV)=0
        IF((LEV.EQ.LV1RAS).AND.(N.GE.LM1RAS)) IVER(IV)=1
        IF((LEV.EQ.LV3RAS).AND.(N.GE.LM3RAS)) IVER(IV)=IVER(IV)+2
      END DO
!
!     NOW LOOP FORWARDS, MARKING THOSE VERTICES CONNECTED FROM ABOVE.
!     SINCE IVER WAS INITIALIZED TO ZERO, NO CHECKING IS NEEDED.
!
      DO IV=1,NVERT0-1
        IVV=IVER(IV)
        DO IC=0,3
          ID=IDOWN0(IV,IC)
          IF(ID.EQ.0) Cycle
          IVER(ID)=IOR(IVER(ID),IVV)
        END DO
      END DO
!
!     THEN LOOP BACKWARDS. SAME RULES, EXCEPT THAT CONNECTIVITY
!     SHOULD BE PROPAGATED ONLY ABOVE THE RESTRICTION LEVELS.
!
      DO IV=NVERT0-1,1,-1
        LEV=IDRT0(IV,LTAB)
        MASK=0
        IF(LEV.GT.LV1RAS) MASK=1
        IF(LEV.GT.LV3RAS) MASK=MASK+2
        IVV=IVER(IV)
        DO IC=0,3
          ID=IDOWN0(IV,IC)
          IF(ID.EQ.0) Cycle
          IVD=IVER(ID)
          IVV=IOR(IVV,IAND(MASK,IVD))
        END DO
        IVER(IV)=IVV
      END DO
!
!     WE ARE NOW INTERESTED ONLY IN VERTICES CONNECTED BOTH TO
!     ALLOWED VERTICES FOR RAS-SPACE 1 AND RAS-SPACE 3.
!     THOSE ARE NUMBERED IN ASCENDING ORDER, THE REST ARE ZEROED.
!
      NVERT=0
      DO IV=1,NVERT0
        IF(IVER(IV).EQ.3) THEN
          NVERT=NVERT+1
          IVER(IV)=NVERT
        ELSE
          IVER(IV)=0
        ENDIF
      END DO
      IF (NVERT.eq.0) Call SysAbendMsg('Restr','No configuration was found\n',    &
                                       'Check NACTEL, RAS1, RAS2, RAS3 values')
!
      END SUBROUTINE RESTR
