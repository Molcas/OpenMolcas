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
!#define _DEBUGPRINT_
      SUBROUTINE MKDAW(NVERT,IDOWN,IDAW)
!     PURPOSE: CONSTRUCT DIRECT ARC WEIGHTS TABLE
!
#ifdef _DEBUGPRINT_
      use Definitions, only: LF => u6
#endif
      IMPLICIT None
!
      Integer NVERT
      Integer IDOWN(NVERT,0:3),IDAW(NVERT,0:4)

      Integer IC, IV, ISUM, IDWN
!
!     BEGIN TO CONSTRUCT DOWN CHAIN TABLE
!
      DO IC=0,3
       IDAW(NVERT,IC)=0
      END DO
      IDAW(NVERT,4)=1
      DO IV=NVERT-1,1,-1
        ISUM=0
        DO IC=0,3
          IDAW(IV,IC)=0
          IDWN=IDOWN(IV,IC)
          IF(IDWN.EQ.0) Cycle
          IDAW(IV,IC)=ISUM
          ISUM=ISUM+IDAW(IDWN,4)
        END DO
        IDAW(IV,4)=ISUM
      END DO
!
#ifdef _DEBUGPRINT_
      Write(LF,*)
      Write(LF,*)' DIRECT ARC WEIGHTS:'
      DO IV=1,NVERT
        Write(LF,'(1X,I4,5X,5(1X,I6))') IV,(IDAW(IV,IC),IC=0,4)
      END DO
      Write(LF,*)
#endif
      END SUBROUTINE MKDAW
