!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1984,1989-1993, Jeppe Olsen                            *
!***********************************************************************
      SUBROUTINE CSDTMT(IDFTP,ICFTP,DTOC,PSSIGN,IPRNT)
!
! Construct list of prototype combinations in IDFTP
! Construct list of prototype CSF'S, in ICFTP
! Construct matrix expanding prototype CSF's in terms of
! prototype combinations in DTOC
!
!
      use stdalloc, only: mma_allocate, mma_deallocate
      use MCLR_Data, only: MULTSP,MS2P,NTYP,MINOP,NCPCNT,NDPCNT
      IMPLICIT None
      DIMENSION IDFTP(*),ICFTP(*),DTOC(*)
      REAL*8 PSSIGN
      Integer IPRNT

!     local variables
      Integer, Allocatable:: SCR7(:)
      Integer MULTS,MS2,IDTBS,ICSBS,ITP,IOPEN,IFLAG,IDFTP,ICFTP,        &
     &        ICDCBS,NNDET
      REAL*8 DTOC
!./SPINFO/
!

      MULTS = MULTSP
      MS2 = MS2P
!
!*.. Set up determinants and upper determinants
!
      IDTBS = 0 ! dummy initialize
      ICSBS = 0 ! dummy initialize
      DO 20 ITP = 1, NTYP
        IOPEN = MINOP+ITP-1
        IF( ITP .EQ. 1 ) THEN
          IDTBS = 1
          ICSBS = 1
        ELSE
          IDTBS = IDTBS + (IOPEN-1)*NDPCNT(ITP-1)
          ICSBS = ICSBS + (IOPEN-1)*NCPCNT(ITP-1)
        END IF
!
        IF( IOPEN .NE. 0 ) THEN
          CALL mma_allocate(SCR7,IOPEN+1,Label='SCR7')
!. Proto type determinants and upper determinants
          IF( MS2+1 .EQ. MULTS ) THEN
            IFLAG = 2
            CALL SPNCOM_MCLR(scr7,IOPEN,MS2,NNDET,IDFTP(IDTBS),         &
     &                  ICFTP(ICSBS),IFLAG,PSSIGN,IPRNT)
          ELSE
            IFLAG = 1
            CALL SPNCOM_MCLR(scr7,IOPEN,MS2,NNDET,IDFTP(IDTBS),         &
     &                  ICFTP(ICSBS),IFLAG,PSSIGN,IPRNT)
            IFLAG = 3
            CALL SPNCOM_MCLR(scr7,IOPEN,MULTS-1,NNDET,IDFTP(IDTBS),     &
     &                  ICFTP(ICSBS),IFLAG,PSSIGN,IPRNT)
          END IF
          CALL mma_deallocate(SCR7)
        END IF
   20 CONTINUE
!. Matrix expressing csf's in terms of combinations
      ICDCBS =0 ! dummy initialize
      DO 30 ITP = 1, NTYP
        IOPEN = MINOP+ITP-1
        IF( ITP .EQ. 1 ) THEN
          IDTBS = 1
          ICSBS = 1
          ICDCBS =1
        ELSE
          IDTBS = IDTBS + (IOPEN-1)*NDPCNT(ITP-1)
          ICSBS = ICSBS + (IOPEN-1)*NCPCNT(ITP-1)
          ICDCBS = ICDCBS + NDPCNT(ITP-1)*NCPCNT(ITP-1)
        END IF
        IF(NDPCNT(ITP)*NCPCNT(ITP).EQ.0) GOTO 30
        IF(IOPEN .EQ. 0 ) THEN
          DTOC(ICDCBS) = 1.0D0
        ELSE
          CALL CSFDET_MCLR(IOPEN,IDFTP(IDTBS),NDPCNT(ITP),              &
     &               ICFTP(ICSBS),NCPCNT(ITP),DTOC(ICDCBS),             &
     &               PSSIGN,IPRNT)
        END IF
   30 CONTINUE
!

      END SUBROUTINE CSDTMT
