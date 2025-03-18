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
      Integer FUNCTION IABNUS(IASTR,NAEL,IAORD,ITPFSA,ISMFSA,NOCTPA,ZA, &
     &                ISSOA,NSSOA,                                      &
     &                IBSTR,NBEL,IBORD,ITPFSB,ISMFSB,NOCTPB,ZB,         &
     &                ISSOB,NSSOB,                                      &
     &                IOOS,NORB,IGENSG,ISGNA,ISGNB,ISGNAB,              &
     &                PSSIGN,IPSFAC,IPRNT)
!
! A determinant is given by strings IASTR,IBSTR .
! Find number of this determinant
!
! If PSSIGN .ne. 0, the determinant with higher alpha number is picked
! and phase factor IPSFAC calculated. This corresponds to
! configuration order
      IMPLICIT NONE
      INTEGER NAEL
      INTEGER IASTR(NAEL),IAORD(*),ITPFSA(*),ISMFSA(*)
      INTEGER NOCTPA
      INTEGER ZA(*),ISSOA(NOCTPA,*),NSSOA(NOCTPA,*)
      INTEGER NBEL
      INTEGER IBSTR(NBEL),IBORD(*),ITPFSB(*),ISMFSB(*)
      INTEGER NOCTPB
      INTEGER ZB(*),ISSOB(NOCTPB,*),NSSOB(NOCTPB,*)
      INTEGER IOOS(NOCTPA,NOCTPB,*)
      INTEGER NORB,IGENSG
      INTEGER ISGNA(*),ISGNB(*)
      REAL*8 PSSIGN
      INTEGER IPSFAC,IPRNT

!     Local variables
      INTEGER NTEST,IANUM,IBNUM,ISGNAB,IASYM,IBSYM,IATP,IBTP,           &
     &        IAREL,IBREL,ISTRNM
!
! Jeppe Olsen
!
      NTEST =  000
      NTEST = MAX(NTEST,IPRNT)
      IF( NTEST .GT. 300) THEN
       WRITE(6,*) ' >>> IABNUS SPEAKING <<< '
       WRITE(6,*) ' NOCTPA,NOCTPB ', NOCTPA,NOCTPB
       WRITE(6,*) ' ALPHA AND BETA STRING '
       CALL IWRTMA(IASTR,1,NAEL,1,NAEL)
       CALL IWRTMA(IBSTR,1,NBEL,1,NBEL)
      END IF
!.Number of alpha- and beta-string
!             ISTRNM(IOCC,NORB,NEL,Z,NEWORD,IREORD)
      IANUM = ISTRNM(IASTR,NORB,NAEL,ZA,IAORD,1)
      IBNUM = ISTRNM(IBSTR,NORB,NBEL,ZB,IBORD,1)
      IF( NTEST .GE. 10 ) WRITE(6,*) ' IANUM AND IBNUM ',IANUM,IBNUM
!
      IF(IGENSG.NE.0) THEN
        ISGNAB = ISGNA(IANUM)*ISGNB(IBNUM)
      ELSE
        ISGNAB = 1
      END IF
!. Symmetries and types
      IASYM = ISMFSA(IANUM)
      IBSYM = ISMFSB(IBNUM)
!?    IF( NTEST .GE.10) WRITE(6,*) ' IASYM IBSYM ',IASYM,IBSYM
      IATP = ITPFSA(IANUM)
      IBTP = ITPFSB(IBNUM)
!?    IF(NTEST.GE.10) WRITE(6,*) ' IATP,IBTP ', IATP,IBTP
      IAREL = IANUM - ISSOA(IATP,IASYM)+1
      IBREL = IBNUM - ISSOB(IBTP,IBSYM)+1
!?    IF(NTEST .GE.10) WRITE(6,*) ' IAREL IBREL ', IAREL,IBREL
!
      IF(PSSIGN.EQ.0.0D0) THEN
!.      Normal determinant ordering
        IABNUS = IOOS(IATP,IBTP,IASYM)                                  &
     &         + (IBREL-1)*NSSOA(IATP,IASYM) + IAREL - 1
        IPSFAC = 1
      ELSE IF (PSSIGN .NE. 0.0D0 ) THEN
!.      Ensure mapping to proper determinant in combination
        IF(IANUM.GE.IBNUM) THEN
!.        No need for switching around so
          IF(IASYM.EQ.IBSYM .AND. IATP .EQ. IBTP ) THEN
!.          Lower triangular packed, column wise !
            IABNUS = IOOS(IATP,IBTP,IASYM)  -1                          &
     &             + (IBREL-1)*NSSOA(IATP,IASYM) + IAREL                &
     &             -  IBREL*(IBREL-1)/2
          ELSE
            IABNUS = IOOS(IATP,IBTP,IASYM)                              &
     &             + (IBREL-1)*NSSOA(IATP,IASYM) + IAREL - 1
          END IF
          IPSFAC = 1
        ELSE IF (IBNUM .GT. IANUM ) THEN
!. Switch alpha and beta string around
          IF(IASYM.EQ.IBSYM .AND. IATP .EQ. IBTP ) THEN
!. Lower triangular packed, column wise !
            IABNUS = IOOS(IBTP,IATP,IBSYM)  -1                          &
     &             + (IAREL-1)*NSSOB(IBTP,IBSYM) + IBREL                &
     &             -  IAREL*(IAREL-1)/2
          ELSE
            IABNUS = IOOS(IBTP,IATP,IBSYM)                              &
     &             + (IAREL-1)*NSSOB(IBTP,IBSYM) + IBREL                &
     &             -  1
          END IF
          IPSFAC = nInt(PSSIGN)
        END IF

      END IF
!
!OLD
!OLD    IABNUS = IOOS(IATP,IBTP,IASYM) + (IBREL-1)*NSSOA(IATP,IASYM)
!OLD &           + IAREL - 1
!?    IF(NTEST .GT. 10 ) then
!?      WRITE(6,*) ' IOOS NSSOA ',IOOS(IATP,IBTP,IASYM),
!?   &              NSSOA(IATP,IASYM)
!?    END IF
!
      IF ( NTEST .GE.200) THEN
         WRITE(6,*) ' ALPHA AND BETA STRING '
         CALL IWRTMA(IASTR,1,NAEL,1,NAEL)
         CALL IWRTMA(IBSTR,1,NBEL,1,NBEL)
         WRITE(6,*) ' Corresponding determinant number ', IABNUS
      END IF

!
      END FUNCTION IABNUS
