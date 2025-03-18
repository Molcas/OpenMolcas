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
      SUBROUTINE CISIZE(NORB1,NORB2,NORB3,                              &
     &                  NEL1MN,NEL3MX,NACTEL,                           &
     &                  MINOP,MAXOP,                                    &
     &                  MXPCNT,MXPCSM,                                  &
     &                  NCNATS,NCNASM,NDTASM,NCSASM,                    &
     &                  NDPCNT,NCPCNT,                                  &
     &                  IICL,IIOP,IIOC,IPRNT)
!
!   Number of configurations per per configuration type and symmetry
!
!  Jeppe Olsen
!         August 1990 : Improved handling of large RAS 3 space
!         Winter 1991 : Modified for LUCIA
      IMPLICIT NONE
      INTEGER NORB1,NORB2,NORB3,NEL1MN,NEL3MX,NACTEL,MINOP,MAXOP,       &
     &        MXPCNT,MXPCSM
!.Output
      INTEGER NCNATS(MXPCNT,*),NCNASM(*),NDTASM(*),NCSASM(*)
!.Input
      INTEGER NDPCNT(*),NCPCNT(*)
!. Scratch
      INTEGER IICL(*),IIOP(*),IIOC(NORB1+NORB2+NORB3)
      INTEGER IPRNT
!
!     Local variables
      Logical Test
      Integer NTEST,ILOOP,ILOOP2,NCNF,NORBT,IORB1F,IORB1L,IORB2F,IORB2L,&
     &        IORB3F,IORB3L,NORB,MINCL1,NOP,ITYPE,NCL,ICL,IFRSTC,IORB,  &
     &        IPLACE,IPRORB,NEWORB,IEL1C,IEL3C,ICL1,IIICHK,MXMPTY,IOP,  &
     &        IFRSTO,IEL1,IEL3,IR3CHK,IFSTR3,K,KEL,KORB,ISYM,I,NTYP,    &
     &        ICSM,ISYMCN_MCLR
!
      NTEST = 0000
      NTEST = MAX(NTEST,IPRNT)
      ILOOP = 0
      ILOOP2 = 0
      NCNF = 0
      NORBT=NORB1+NORB2+NORB3

!
      CALL iCOPY(MXPCSM*MXPCNT,[0],0,NCNATS,1)
      CALL iCOPY(MXPCSM,[0],0,NCSASM,1)
      CALL iCOPY(MXPCSM,[0],0,NDTASM,1)
      CALL iCOPY(MXPCSM,[0],0,NCNASM,1)

!
      IORB1F = 1
      IORB1L = IORB1F+NORB1-1
!
      IORB2F = IORB1L + 1
      IORB2L = IORB2F + NORB2 - 1
!
      IORB3F = IORB2L + 1
      IORB3L = IORB3F + NORB3 - 1
!
      NORB = NORB1 + NORB2 + NORB3
! Min number of doubly occupied orbitals in RAS 1
      MINCL1 = MAX(0,NEL1MN-NORB1)
      IF(NTEST.GE.1)  WRITE(6,*)                                        &
     &  ' Min number of doubly occupied orbitals in RAS 1',MINCL1
      DO 5000 NOP = MINOP,MAXOP,2
        ITYPE = NOP-MINOP+1
        NCL = (NACTEL-NOP)/2
        IF( NTEST .GE. 10 )                                             &
     &  WRITE(6,*) ' NOP NCL ITYPE',NOP,NCL,ITYPE
!. first combination of double occupied orbitals
        CALL iCOPY(NORB,[0],0,IIOC,1)
        DO 10 ICL = 1, NCL
          IICL(ICL) = ICL
          IIOC(ICL) = 2
   10   CONTINUE
        IFRSTC = 1
!.. Loop over double occupied orbital configurations
 2000   CONTINUE
!
!. next double occupied configuration
          IF ( IFRSTC .EQ. 1 .OR. NCL .EQ. 0 ) GOTO 801
!         IF ( IFRSTC .EQ. 0 .AND. NCL .NE. 0 ) THEN
!
          DO 50 IORB = 1, NORB
            IF(IIOC(IORB) .EQ. 1 ) IIOC(IORB) = 0
   50     CONTINUE
!
          IPLACE = 0
  800     IPLACE = IPLACE + 1

          IPRORB = IICL(IPLACE)
          IIOC(IPRORB) = 0
          NEWORB = IPRORB+1
          IF((IPLACE .LT. NCL .AND. NEWORB .LT. IICL(IPLACE+1))         &
     &      .OR.                                                        &
     &      IPLACE .EQ. NCL .AND.  NEWORB.LE. NORB ) THEN
            IICL(IPLACE) = NEWORB
            IIOC(NEWORB) = 2
          ELSE IF                                                       &
     &    (.NOT.(IPLACE.EQ.NCL.AND.NEWORB.GE.NORB)) THEN
            IF(IPLACE .EQ. 1 ) THEN
              IICL(1) = 1
              IIOC(1) = 2
            ELSE
              IICL(IPLACE) = IICL(IPLACE-1) + 1
              IIOC(IICL(IPLACE)) = 2
            END IF
            GOTO 800
          ELSE
!. No more inactive configurations
             GOTO 2001
          END IF
!OLD      END IF
  801     CONTINUE
          IFRSTC = 0
          IF( NTEST .GE.1500) THEN
            WRITE(6,*) ' Next inactive configuration '
            CALL IWRTMA(IICL,1,NCL,1,NCL)
          END IF
!..         CHECK RAS1 and RAS 3
            IEL1C = 0
            IEL3C = 0
            ICL1  = 0
            DO  20 ICL = 1,NCL
              IORB = IICL(ICL)
              IF(IORB1F.LE.IORB .AND. IORB.LE.IORB1L ) THEN
                 IEL1C = IEL1C + 2
                 ICL1 = ICL1 + 1
              ELSE IF (IORB3F .LE. IORB .AND. IORB .LE. IORB3L)THEN
                 IEL3C = IEL3C + 2
              END IF
   20       CONTINUE
            IIICHK = 1
            IF(ICL1 .LT.MINCL1.AND. IIICHK.EQ.1) THEN
! Next higher combination with a higher number of inactive orbitals
              DO 41 ICL = 1,ICL1+1
                IIOC(IICL(ICL)) = 0
                IICL(ICL) = ICL
                IIOC(ICL) = 2
   41         CONTINUE
              IPLACE=ICL1+1
              IF( IPLACE.GE.NCL) GOTO 2001
              GOTO 800
            END IF
            IF(IEL3C.GT.NEL3MX) GOTO 2000
!. Highest orbital not occupied
         MXMPTY = NORB
         IORB = NORB+1
!. begin while
   12      CONTINUE
           IORB = IORB - 1
           IF(IIOC(IORB) .EQ. 2 ) THEN
             MXMPTY = IORB-1
             IF( IORB .NE. 1 ) GOTO 12
           END IF
!. End while
!
!. first active configuration
          IORB = 0
          IOP = 0
          DO 30 IORB = 1, NORB
            IF(IIOC(IORB) .EQ. 0 ) THEN
              IOP = IOP + 1
              IF( IOP .GT. NOP ) GOTO 31
              IIOC(IORB) = 1
              IIOP(IOP) = IORB
            END IF
   30     CONTINUE
   31     CONTINUE
          IFRSTO = 1
!
!. Next open shell configuration
 1000     CONTINUE
            IF(IFRSTO.EQ.1 .OR. NOP .EQ. 0 ) goto 701
!OLD        IF(IFRSTO .EQ. 0 .AND. NOP .NE. 0 ) THEN
            IPLACE = 0
  700       CONTINUE
              IPLACE = IPLACE + 1
              IPRORB = IIOP(IPLACE)
              NEWORB = IPRORB + 1
              IIOC(IPRORB) = 0
  690         CONTINUE
                IF(NEWORB .LE. MXMPTY .AND.                             &
     &             IIOC(MIN(NORBT,NEWORB)) .NE. 0) THEN
                   NEWORB = NEWORB + 1
                   IF(NEWORB .LE. MXMPTY ) GOTO 690
                END IF
! 691         End of loop
            Test=IPLACE .LT. NOP
            If (Test) Test=NEWORB .LT. IIOP(IPLACE+1)
            IF(Test .OR.                                                &
     &        IPLACE .EQ. NOP .AND. NEWORB .LE. MXMPTY ) THEN
              IIOP(IPLACE) = NEWORB
              IIOC(NEWORB) = 1
           ELSE IF (  IPLACE .NE. NOP ) THEN
              IF(IPLACE.EQ.1) THEN
                NEWORB = 1 - 1
              ELSE
                NEWORB = IIOP(IPLACE-1)
              END IF
 671          CONTINUE
                NEWORB = NEWORB + 1
              IF(IIOC(NEWORB) .NE. 0 .AND. NEWORB.LE.MXMPTY ) GOTO 671
              IIOP(IPLACE) = NEWORB
              IIOC(NEWORB) = 1
              GOTO 700
           ELSE
!. No more active configurations , so
              IF( NCL .NE. 0 ) GOTO 2000
              IF( NCL .EQ. 0 ) GOTO 5001
           END IF
!          END IF
  701      CONTINUE
           IFRSTO = 0

          IF( NTEST .GE.1500) THEN
            WRITE(6,*) ' Next active configuration '
            CALL IWRTMA(IIOP,1,NOP,1,NOP)
          END IF
!        RAS  CONSTRAINTS
           IEL1 = IEL1C
           IEL3 = IEL3C
!..        CHECK RAS1 and RAS3
           DO  40 IOP = 1,NOP
             IORB = IIOP(IOP)
             IF(IORB1F.LE.IORB .AND. IORB.LE.IORB1L ) THEN
                IEL1 = IEL1 + 1
             ELSE IF (IORB3F .LE. IORB .AND. IORB .LE. IORB3L)THEN
                IEL3 = IEL3 + 1
             END IF
   40      CONTINUE
           IR3CHK  = 1
           IF(IEL3.GT.NEL3MX.AND.IR3CHK.EQ.1) THEN
!. Number of electrons in substring
             IFSTR3 = 0
             DO 5607 IOP = 1, NOP
               IF(IIOP(IOP).GE.IORB3F) THEN
                IFSTR3 = IOP
                GOTO 5608
               END IF
 5607        CONTINUE
 5608        CONTINUE
             IF(IFSTR3.NE.NOP) THEN

!. Lowest possible string with NOP electrons
               DO 5610 K = 1, IFSTR3
                 IIOC(IIOP(K)) = 0
 5610          CONTINUE
!
               KEL = 0
               KORB = 0
 5630          CONTINUE
                 KORB = KORB + 1
                 IF(IIOC(KORB).NE.2) THEN
                   KEL = KEL + 1
                   IIOC(KORB) = 1
                   IIOP(KEL) = KORB
                 END IF
               IF(KEL.NE.IFSTR3) GOTO  5630
               IPLACE = IFSTR3
               GOTO 700
             END IF
           END IF
           IF( IEL1 .LT. NEL1MN .OR. IEL3 .GT. NEL3MX ) GOTO  999
!. Spatial symmetry
         ISYM = ISYMCN_MCLR(IICL,IIOP,NCL,NOP)
!
           IF( NTEST .GE.2000)                                          &
     &     WRITE(6,*) ' ISYM : ', ISYM
           IF(NTEST.GE.1500)                                            &
     &     WRITE(6,1120) ( IIOC(I),I = 1,NORB )
 1120      FORMAT('0  configuration included ',15I3,                    &
     &               ('                         ',15I3))
           NCNF=NCNF+1

           NCNASM(ISYM) = NCNASM(ISYM)+1
           IF(NTEST.GE.1500 )                                           &
     &     WRITE(6,1311) NCNF,(IIOC(I),I=1,NORB)
 1311      FORMAT('  configuration ',I3,20I2,/,(1X,18X,20I2))

           NCNATS(ITYPE,ISYM)=NCNATS(ITYPE,ISYM)+1
           IF(NTEST.GE.2000) WRITE(6,3111) NCNF,ITYPE
 3111      FORMAT('0  CONFIGURATION..',I3,' IS TYPE..',I3)
!
!** LOOP OVER CONFIGURATIONS, end
!
  999  CONTINUE
       ILOOP = ILOOP + 1
       ILOOP2 = ILOOP2 + 1
       IF(ILOOP2 .EQ. 10000000 ) THEN
         WRITE(6,*) ' 10 million configurations generated '
         ILOOP2 = 0
       END IF
!
       IF( NOP .EQ. 0 .AND. NCL .EQ. 0 ) GOTO 5001
       IF( NOP .EQ. 0 ) GOTO 2000
       GOTO 1000
 2001 CONTINUE
 5000 CONTINUE
 5001 CONTINUE

      IF( NTEST .GE. 2 ) THEN
        WRITE(6,'(A,I8)')                                               &
     &  '  Total number of configurations generated ', NCNF
      END IF
! ================================
!. Total number of CSF's and SD's
! ================================
      NTYP = MAXOP - MINOP + 1
      DO 610 ISYM = 1, MXPCSM
        DO 600 ITYPE = 1, NTYP
          NDTASM(ISYM) = NDTASM(ISYM)                                   &
     &    + NDPCNT(ITYPE)*NCNATS(ITYPE,ISYM)
          NCSASM(ISYM) = NCSASM(ISYM)                                   &
     &    + NCPCNT(ITYPE)*NCNATS(ITYPE,ISYM)
  600   CONTINUE
  610 CONTINUE

!
      IF( NTEST .GE. 2 ) THEN
       WRITE(6,*)
       ICSM = 0
       WRITE(6,'(/A)') ' Information about actual configurations '
       WRITE(6,'( A)') ' ========================================'
       WRITE(6,'(/A)')                                                  &
     & '    Symmetry     Configurations     CSFs     Combinations  '
       WRITE(6,'(A)')                                                   &
     & '  =============  ============== ============ ============  '
       DO 570 ICSM = 1, MXPCSM
         IF(NCNASM(ICSM) .NE. 0 ) THEN
            WRITE(6,'(4X,I3,4X,6X,I8,6X,I8,6X,I9)')                     &
     &      ICSM ,NCNASM(ICSM),NCSASM(ICSM),NDTASM(ICSM)
         END IF
  570  CONTINUE
      END IF
!
      END SUBROUTINE CISIZE
