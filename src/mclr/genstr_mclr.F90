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
! Copyright (C) 1990, Jeppe Olsen                                      *
!***********************************************************************
      SUBROUTINE GENSTR_MCLR(NEL,NELMN1,NELMX1,NELMN3,NELMX3,           &
     &                  ISTASO,NOCTYP,NSMST,Z,LSTASO,                   &
     &                  IREORD,STRING,IOC,IOTYP,IPRNT)
!
! Generate strings consisting of  NEL electrons fulfilling
!   1 : Between NELMN1 AND NELMX1 electrons in the first NORB1 orbitals
!   2 : Between NELMN3 AND NELMX3 electrons in the last  NORB3 orbitals
!
! In the present version the strings are directly ordered into
! symmetry and occupation type .
!
! Jeppe Olsen Winter of 1990
! ========
! Output :
! ========
! STRING(IEL,ISTRIN) : Occupation of strings.
! IREORD             : Reordering array going from lexical
!                      order to symmetry and occupation type order.
!
      use MCLR_Data, only: NACOB,NORB1,NORB2,NORB3
      IMPLICIT None
      INTEGER NEL,NELMN1,NELMX1,NELMN3,NELMX3,NOCTYP,NSMST
      INTEGER IOTYP,IPRNT
!. Input
      Integer ISTASO(NOCTYP,NSMST)
      INTEGER Z(NACOB,NEL)
!
!.Output
      INTEGER STRING(NEL,*),IREORD(*)
!.Scratch arrays
      Integer IOC(*),LSTASO(NOCTYP,NSMST)
!
!
!     Local variables
      Integer NTEST0,NTEST,NSTRIN,IORB1F,IORB1L,IORB2F,IORB2L,IORB3F,   &
     &        IORB3L,IEL1,IEL2,IEL3,IFRST1,IFRST2,IFRST3,NONEW1,NONEW2, &
     &        NONEW3,ISYM,ITYP,LEXCI,LACTU,NPR,ISTRIN,LSTRIN,KSTRIN,    &
     &        IEL,IOCTP2_MCLR,ISTRNM,ISYMST_MCLR,i
      NTEST0 = 0
      NTEST = MAX(NTEST0,IPRNT)
!
      Call iCopy(NOCTYP*NSMST,[0],0,LSTASO,1)
      NSTRIN = 0
      IORB1F = 1
      IORB1L = IORB1F+NORB1-1
      IORB2F = IORB1L + 1
      IORB2L = IORB2F+NORB2-1
      IORB3F = IORB2L + 1
      IORB3L = IORB3F+NORB3-1
! Loop over possible partitionings between RAS1,RAS2,RAS3
      DO 1001 IEL1 = NELMX1,NELMN1,-1
      DO 1003 IEL3 = NELMN3,NELMX3, 1
       IF(IEL1.GT. NORB1 ) GOTO 1001
       IF(IEL3.GT. NORB3 ) GOTO 1003
       IEL2 = NEL - IEL1-IEL3
       IF(IEL2 .LT. 0 .OR. IEL2 .GT. NORB2 ) GOTO 1003
       IFRST1 = 1
! Loop over RAS 1 occupancies
901    CONTINUE
         IF( IEL1 .NE. 0 ) THEN
           IF(IFRST1.EQ.1) THEN
            IOC(1:IEL1) = [(i,i=1,IEL1)]
            IFRST1 = 0
           ELSE
             CALL NXTORD(IOC,IEL1,IORB1F,IORB1L,NONEW1)
             IF(NONEW1 .EQ. 1 ) GOTO 1003
           END IF
         END IF
         IF( NTEST .GE. 500) THEN
           WRITE(6,*) ' RAS 1 string '
           CALL IWRTMA(IOC,1,IEL1,1,IEL1)
         END IF
         IFRST2 = 1
         IFRST3 = 1
! Loop over RAS 2 occupancies
902      CONTINUE
           IF( IEL2 .NE. 0 ) THEN
             IF(IFRST2.EQ.1) THEN
              IOC(IEL1+1:IEL1+IEL2) = [(i,i=IORB2F,IORB2F+IEL2-1)]
              IFRST2 = 0
             ELSE
               CALL NXTORD(IOC(IEL1+1),IEL2,IORB2F,IORB2L,NONEW2)
               IF(NONEW2 .EQ. 1 ) THEN
                 IF(IEL1 .NE. 0 ) GOTO 901
                 IF(IEL1 .EQ. 0 ) GOTO 1003
               END IF
             END IF
           END IF
           IF( NTEST .GE. 500) THEN
             WRITE(6,*) ' RAS 1 2 string '
             CALL IWRTMA(IOC,1,IEL1+IEL2,1,IEL1+IEL2)
           END IF
           IFRST3 = 1
! Loop over RAS 3 occupancies
903        CONTINUE
             IF( IEL3 .NE. 0 ) THEN
               IF(IFRST3.EQ.1) THEN
                IOC(IEL1+IEL2+1:IEL1+IEL2+IEL3) =                       &
     &            [(i,i=IORB3F,IORB3F+IEL3-1)]
                IFRST3 = 0
               ELSE
                 CALL NXTORD(IOC(IEL1+IEL2+1),                          &
     &           IEL3,IORB3F,IORB3L,NONEW3)
                 IF(NONEW3 .EQ. 1 ) THEN
                   IF(IEL2 .NE. 0 ) GOTO 902
                   IF(IEL1 .NE. 0 ) GOTO 901
                   GOTO 1003
                 END IF
               END IF
             END IF
             IF( NTEST .GE. 500 ) THEN
               WRITE(6,*) ' RAS 1 2 3 string '
               CALL IWRTMA(IOC,1,NEL,1,NEL)
             END IF
! Next string has been constructed , Enlist it !.
             NSTRIN = NSTRIN + 1
!. Symmetry
             ISYM = ISYMST_MCLR(IOC,NEL)
!. Occupation type
             ITYP = IOCTP2_MCLR(IOC,NEL,IOTYP)
!
             IF(ITYP.NE.0) THEN
               LSTASO(ITYP,ISYM) = LSTASO(ITYP,ISYM)+ 1
               LEXCI = ISTRNM(IOC,NACOB,NEL,Z,IREORD,0)
               LACTU = ISTASO(ITYP,ISYM)-1+LSTASO(ITYP,ISYM)
               IREORD(LEXCI) = LACTU
               If (NEL.ne.0) Call iCopy(NEL,IOC,1,STRING(1,LACTU),1)
             END IF
!
           IF( IEL3 .NE. 0 ) GOTO 903
           IF( IEL3 .EQ. 0 .AND. IEL2 .NE. 0 ) GOTO 902
           IF( IEL3 .EQ. 0 .AND. IEL2 .EQ. 0 .AND. IEL1 .NE. 0)         &
     &     GOTO 901
1003  CONTINUE
1001  CONTINUE
!
      IF(NTEST.GE.1 ) THEN
        WRITE(6,*) ' Number of strings generated   ', NSTRIN
      END IF
      IF(NTEST.GE.10)  THEN
        IF(NTEST.GE.100) THEN
          NPR = NSTRIN
        ELSE
          NPR = MIN(NSTRIN,50)
        END IF
        WRITE(6,*) ' Strings generated '
        WRITE(6,*) ' =================='
        ISTRIN = 0
        DO 100 ISYM = 1, NSMST
        DO 101 ITYP = 1,NOCTYP
          LSTRIN = MIN(LSTASO(ITYP,ISYM),NPR-ISTRIN)
          IF(LSTRIN.GT.0) THEN
            WRITE(6,*) ' Strings of type and symmetry ',ITYP,ISYM
            DO 90 KSTRIN = 1,LSTRIN
              ISTRIN = ISTRIN + 1
              WRITE(6,'(2X,I4,8X,(10I5))')                              &
     &        ISTRIN,(STRING(IEL,ISTRIN),IEL = 1,NEL)
90          CONTINUE
          END IF
101     CONTINUE
100     CONTINUE
!
        WRITE(6,*) ' Array giving actual place from lexical place'
        WRITE(6,*) ' ============================================'
        CALL IWRTMA(IREORD,1,NPR,1,NPR)
      END IF
!
      END SUBROUTINE GENSTR_MCLR
