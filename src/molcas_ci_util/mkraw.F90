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
      SUBROUTINE MKRAW(NVERT,IDOWN,IUP,IRAW)
!
!     PURPOSE: CONSTRUCT UPCHAIN INDEX TABLE AND REVERSE ARC WEIGHTS
!
#ifdef _DEBUGPRINT_
      use Definitions, only: LF => u6
#endif
      IMPLICIT None
!
!
      Integer NVERT
      Integer IDOWN(NVERT,0:3),IUP(NVERT,0:3),IRAW(NVERT,0:4)

      Integer IU, IC, IDWN, IV, ISUM
!
!     BEGIN BY CONSTRUCTING THE UPCHAIN TABLE IUP:
!
      IUP(:,:)=0
      DO IU=1,NVERT-1
        DO IC=0,3
          IDWN=IDOWN(IU,IC)
          IF(IDWN.EQ.0) Cycle
          IUP(IDWN,IC)=IU
        END DO
      END DO
!
#ifdef _DEBUGPRINT_
      Write(LF,*)
      Write(LF,*)' THE UPCHAIN TABLE IN MKRAW:'
      DO IV=1,NVERT
        Write(LF,'(1X,I4,5X,4(1X,I6))') IV,(IUP(IV,IC),IC=0,3)
      END DO
      Write(LF,*)
#endif
!
!     USE UPCHAIN TABLE TO CALCULATE THE REVERSE ARC WEIGHT TABLE:
!
      IRAW(1,0:3)=0
      IRAW(1,4)=1
      DO IV=2,NVERT
        ISUM=0
        DO IC=0,3
          IRAW(IV,IC)=0
          IU=IUP(IV,IC)
          IF(IU.EQ.0) Cycle
          IRAW(IV,IC)=ISUM
          ISUM=ISUM+IRAW(IU,4)
        END DO
        IRAW(IV,4)=ISUM
      END DO
!
#ifdef _DEBUGPRINT_
      Write(LF,*)
      Write(LF,*)' THE REVERSE ARC WEIGHT TABLE IN MKRAW:'
      DO IV=1,NVERT
        Write(LF,'(1X,I4,5X,5(1X,I6))') IV,(IRAW(IV,IC),IC=0,4)
      END DO
      Write(LF,*)
#endif

      END SUBROUTINE MKRAW
