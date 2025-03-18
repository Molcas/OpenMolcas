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
      SUBROUTINE CONFG2 (NORB1,NORB2,NORB3,                             &
     &                   NEL1MN,NEL3MX,                                 &
     &                   MINOP,MAXOP,                                   &
     &                   IREFSM,NEL,ICONF,                              &
     &                   NCNFTP,IIOC,IIOP,IICL,IPRNT)
!
!  Generate array,ICONF,giving occupation of each configuration
!  for CI space of reference symmetry IREFSM.
!
!
!  Jeppe Olsen April 1989
!              August 1990 : Improved handling of large RAS 3 space
!
!  Turbo configuration generator
!  Iconf is ordered so all configuratiuns of the same type are
!  consecutively stored .
!  ICONF is written so closed orbitals are given first and then single
!  occupied orbitals
!
      IMPLICIT None
      Integer NORB1,NORB2,NORB3,NEL1MN,NEL3MX,MINOP,MAXOP,IREFSM,NEL
!.Output
      Integer, Intent(Out)::ICONF(*)
!.Input
      Integer, Intent(In):: NCNFTP(*)
!.Scratch
      Integer IIOC(*),IICL(*),IIOP(*)
      Integer IPRNT
!
! local variables
      Logical Test
      Integer NTEST,IORB1F,IORB1L,IORB2F,IORB2L,IORB3F,IORB3L,NORB,     &
     &        JCONF,ICFREE,MINCL1,NOP,ITYPE,NCL,ICL,IFRSTC,IORB,        &
     &        IPLACE,IPRORB,NEWORB,IEL1C,IEL3C,ICL1,IIICHK,MXMPTY,      &
     &        IOP,IFRSTO,IEL1,IEL3,IR3CHK,IFSTR3,K,KEL,KORB,ISYM,I,     &
     &        IBAS,IOPEN,IOC,LICONF,ISYMCN_MCLR
!
      NTEST = 0000

      NTEST = MAX(NTEST,IPRNT)
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
!
!.. Loop over types of configurations
!
      JCONF =  0
      ICFREE = 1
! Min number of doubly occupied orbitals in RAS 1
      MINCL1 = MAX(0,NEL1MN-NORB1)
      IF(NTEST.GE.1) WRITE(6,*)                                         &
     &  ' Min number of doubly occupied orbitals in RAS 1',MINCL1
      DO 5000 NOP = MINOP,MAXOP,2
        ITYPE = NOP-MINOP+1
        NCL = (NEL-NOP)/2
        IF( NTEST .GE. 10 )                                             &
     &  WRITE(6,*) ' NOP NCL ITYPE',NOP,NCL,ITYPE
!
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
          IF( IFRSTC .EQ. 1 .OR. NCL .EQ. 0 ) GOTO 801
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
!
            IICL(IPLACE) = NEWORB
            IIOC(NEWORB) = 2
          ELSE IF                                                       &
     &    (.NOT.(IPLACE.EQ.NCL.AND.NEWORB.GE.NORB)) THEN
!
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
  801     CONTINUE
          IFRSTC = 0
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
          IF( NTEST .GE.1500) THEN
            WRITE(6,*) ' Next inactive configuration '
            CALL IWRTMA(IICL,1,NCL,1,NCL)
          END IF
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
            IF(IFRSTO .EQ. 1 .OR. NOP .EQ. 0 ) GOTO 701
            IPLACE = 0
  700       CONTINUE
              IPLACE = IPLACE + 1
              IPRORB = IIOP(IPLACE)
              NEWORB = IPRORB + 1
              IIOC(IPRORB) = 0

! PAM 2013: Searching for next orbital with IIOC=0:
  690         CONTINUE
              Test = NEWORB.LE.MXMPTY
              If (Test) Test=IIOC(NEWORB) .NE. 0
              IF(Test) THEN
                NEWORB = NEWORB + 1
                GOTO 690
              END IF

            Test = IPLACE .LT. NOP
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
              IF(IIOC(NEWORB) .NE. 0 .AND. NEWORB.LT.MXMPTY ) GOTO 671
              IIOP(IPLACE) = NEWORB
              IIOC(NEWORB) = 1
              GOTO 700
           ELSE
!. No more active configurations , so
              IF( NCL .NE. 0 ) GOTO 2000
              IF( NCL .EQ. 0 ) GOTO 5001
           END IF
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
!. Faster routine for RAS 3, added august 1990
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
         IF(ISYM.EQ.IREFSM) THEN
           IF(NTEST.GE.100)                                             &
     &     WRITE(6,1120) ( IIOC(I),I = 1,NORB )
 1120      FORMAT('0  configuration included ',15I3)
           JCONF=JCONF+1
!
           DO 60 ICL = 1, NCL
             ICONF(ICFREE-1+ICL) = IICL(ICL)
   60      CONTINUE
           DO 61 IOP = 1, NOP
             ICONF(ICFREE-1+NCL+IOP) = IIOP(IOP)
   61      CONTINUE
           ICFREE = ICFREE + NOP + NCL
         END IF

!
!** LOOP OVER active configurations end
!
  999   CONTINUE
        IF( NOP .EQ. 0 .AND. NCL .EQ. 0 ) GOTO 5001
        IF( NOP .EQ. 0 ) GOTO 2000
        GOTO 1000
 2001 CONTINUE
 5000 CONTINUE
 5001 CONTINUE

!
      IF( NTEST .GE. 100)THEN
        WRITE(6,'(/A,I3)') '  Configurations of symmetry ', IREFSM
        WRITE(6,*)        ' ================================='
        IBAS = 0
        DO 1200 IOPEN = MINOP,MAXOP
          ITYPE = IOPEN - MINOP + 1
          ICL = (NEL-IOPEN)/2
          IOC = IOPEN + ICL
          LICONF = NCNFTP(ITYPE)
          WRITE(6,'(/A,2I3)')                                           &
     &    '  Type with number of closed and open orbitals ',ICL,IOPEN
          WRITE(6,'(A,I7)')                                             &
     &    '  Number of configurations of this type',LICONF
          DO 1180 JCONF = 1,LICONF
            WRITE(6,'(3X,20I3)') (ICONF(IBAS+IORB),IORB=1,IOC)
            IBAS = IBAS + IOC
 1180     CONTINUE
 1200   CONTINUE
      END IF

!
      END SUBROUTINE CONFG2
