************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1989,1993, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE DIHDJ2_MCLR(IASTR,IBSTR,NIDET,
     &                 JASTR,JBSTR,NJDET,
     &                 NAEL,NBEL,
     & jWORK,LWORK,NORB,HAMIL,ISYM,NINOB,ECORE,ICOMBI,PSIGN,
     & IASTRM,IBSTRM,JASTRM,JBSTRM,
     & IGENSG,IASGN,IBSGN,JASGN,JBSGN,LIA,LIB,NDIF0,NDIF1,NDIF2,
     & IPRT)
*
* A set of left hand side determinants defined by string numbers
* IASTR and IBSTR and a set of right hand side determinants
* defined by JASTR and JBSTR are given.
*
* Obtain Hamiltonian matrix  < IA IB ! H ! JA JB >
*
* If Icombi .NE. 0 Spin combinations are assumed  for alpha and
* beta strings with different orbital configurations
*   1/SQRT(2) * ( !I1A I2B! + PSIGN * !I2A I1B! )
*
* If ISYM .EQ. 0 FULL Hamiltonian is constructed
* If ISYM .NE. 0 LOWER half of hamiltonian is constructed
*
* JEPPE OLSEN JANUARY 1989
*
*. Modifed to work with string numbers instead of strings
*. March 93
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IASTR(*),IBSTR(*)
      DIMENSION JASTR(*),JBSTR(*)
      DIMENSION IASTRM(NAEL,*),IBSTRM(NBEL,*)
      DIMENSION JASTRM(NAEL,*),JBSTRM(NBEL,*)
      DIMENSION IASGN(*),IBSGN(*),JASGN(*),JBSGN(*)
*
      DIMENSION jWORK(*), HAMIL(*)
      DIMENSION LIA(NAEL),LIB(NBEL)
*
*
*. Scratch space : 4 vectors of length NORB
      KLFREE = 1
      KLIAE  = KLFREE
      KLFREE = KLIAE + NORB
      KLIBE  = KLFREE
      KLFREE = KLIBE + NORB
*
      KLJAE = KLFREE
      KLFREE = KLJAE + NORB
      KLJBE = KLFREE
      KLFREE = KLJBE + NORB
*
      IF( ISYM .EQ. 0 ) THEN
        LHAMIL = NIDET*NJDET
      ELSE
        LHAMIL = NIDET*(NIDET+1) / 2
      END IF
      CALL SETVEC(HAMIL,0.0D0,LHAMIL)
*
      NTERMS= 0
      NDIF0 = 0
      NDIF1 = 0
      NDIF2 = 0
*
*. Loop over J determinants
*
      JAEQJB = -1    ! dummy initialize
      CONST  = 0.0D0 ! dummy initialize
      IEL1   = -1    ! dummy initialize
      JEL1   = -1    ! dummy initialize
      SIGNA  = 0.0D0 ! dummy initialize
      SIGNB  = 0.0D0 ! dummy initialize
      IPERM  = 0     ! dummy initialize
      JPERM  = 0     ! dummy initialize
      SIGN   = 0.0D0 ! dummy initialize
      XVAL   = 0.0D0 ! dummy initialize
      DO 1000 JDET = 1,NJDET
* Expand JDET
        JASTAC =JASTR(JDET)
        JBSTAC =JBSTR(JDET)
*
        IF(IGENSG .GT. 0 ) THEN
         JXSGN = JASGN(JASTAC)*JBSGN(JBSTAC)
        ELSE
         JXSGN = 1
        END IF
*
        CALL ISETVC(jWORK(KLJAE),0,NORB)
        CALL ISETVC(jWORK(KLJBE),0,NORB)
        DO 40 IAEL = 1, NAEL
          jWORK(KLJAE-1+JASTRM(IAEL,JASTAC) ) = 1
   40   CONTINUE
*
        DO 50 IBEL = 1, NBEL
          jWORK(KLJBE-1+JBSTRM(IBEL,JBSTAC) ) = 1
   50   CONTINUE
*
        IF( ICOMBI .NE. 0 ) THEN
          IF(JASTAC .EQ. JBSTAC) THEN
             JAEQJB = 1
          ELSE
             JAEQJB = 0
          END IF
        END IF
*
*
*
        IF( ISYM .EQ. 0 ) THEN
          MINI = 1
        ELSE
          MINI = JDET
        END IF
*
* Loop over I determinants
*
        DO 900 IDET = MINI, NIDET
          IASTAC = IASTR(IDET)
          IBSTAC = IBSTR(IDET)
*
          IF(IGENSG .GT. 0 ) THEN
           IXSGN = IASGN(IASTAC)*IBSGN(IBSTAC)
          ELSE
           IXSGN = 1
          END IF
*
          IF(IASTAC.EQ.IBSTAC) THEN
            IAEQIB = 1
          ELSE
            IAEQIB = 0
          END IF

*
          IF(ICOMBI.EQ.1 .AND. IAEQIB+JAEQJB.EQ.0 ) THEN
              NLOOP = 2
          ELSE
              NLOOP = 1
          END IF
          DO 899 ILOOP = 1, NLOOP
           NTERMS = NTERMS + 1
* For second part of spin combinations strings should be swopped
           IF(ILOOP.EQ.1) THEN
             CALL ICOPVE(IASTRM(1,IASTAC),LIA,NAEL)
             CALL ICOPVE(IBSTRM(1,IBSTAC),LIB,NBEL)
           ELSE IF (ILOOP.EQ.2) THEN
             CALL ICOPVE(IASTRM(1,IASTAC),LIB,NAEL)
             CALL ICOPVE(IBSTRM(1,IBSTAC),LIA,NBEL)
           END IF
*
* ==============================
*. Number of orbital differences
* ==============================
*
           NACM = 0
           DO 61 IAEL = 1, NAEL
             NACM = NACM + jWORK(KLJAE-1+LIA(IAEL))
   61      CONTINUE
           NBCM = 0
           DO 62 IBEL = 1, NBEL
             NBCM = NBCM + jWORK(KLJBE-1+LIB(IBEL))
   62      CONTINUE
           NADIF = NAEL-NACM
           NBDIF = NBEL-NBCM
*
           IF(NADIF+NBDIF .GT. 2 ) GOTO 898
*. Factor for combinations
           IF( ICOMBI .EQ. 0 ) THEN
             CONST = 1.0D0
           ELSE
             IF((JAEQJB +IAEQIB) .EQ.2 ) THEN
               CONST = 1.0D0
             ELSE IF( (JAEQJB+IAEQIB) .EQ. 1 ) THEN
               CONST = 1.0D0/SQRT(2.0D0)*(1.0D0+PSIGN)
              ELSE IF( (JAEQJB+IAEQIB) .EQ. 0 ) THEN
               IF( ILOOP .EQ. 1)  THEN
                 CONST = 1.0D0
               ELSE
                 CONST = PSIGN
               END IF
             END IF
           END IF
*. External sign factor
           IF(IXSGN*JXSGN .EQ. -1 ) CONST = - CONST
*
* ==================================================
*.. Find differing orbitals and sign for permutation
* ==================================================
*
* Expand idet
           CALL ISETVC(jWORK(KLIAE),0,NORB)
           CALL ISETVC(jWORK(KLIBE),0,NORB)
*
           DO 42 IAEL = 1, NAEL
             jWORK(KLIAE-1+LIA(IAEL)) = 1
   42      CONTINUE
*
           DO 52 IBEL = 1, NBEL
             jWORK(KLIBE-1+LIB(IBEL) ) = 1
   52      CONTINUE
*
*. One pair of differing alpha electrons
*
           IF(NADIF .EQ. 1 ) THEN
             DO 120 IAEL = 1,NAEL
               IF(jWORK(KLJAE-1+LIA(IAEL)).EQ.0) THEN
                 IA = LIA(IAEL)
                 IEL1 = IAEL
                 GOTO 121
               END IF
  120        CONTINUE
  121        CONTINUE
*
             DO 130 JAEL = 1,NAEL
               IF(jWORK(KLIAE-1+JASTRM(JAEL,JASTAC)).EQ.0) THEN
                 JA = JASTRM(JAEL,JASTAC)
                 JEL1 = JAEL
                 GOTO 131
                END IF
  130        CONTINUE
  131        CONTINUE
             SIGNA = DBLE((-1)**(JEL1+IEL1))
           END IF
*
*. One pair of differing beta electrons
*
           IF(NBDIF .EQ. 1 ) THEN
             DO 220 IBEL = 1,NBEL
               IF(jWORK(KLJBE-1+LIB(IBEL) ).EQ.0) THEN
                 IB = LIB(IBEL)
                 IEL1 = IBEL
                 GOTO 221
                END IF
  220        CONTINUE
  221        CONTINUE
             DO 230 JBEL = 1,NBEL
               IF(jWORK(KLIBE-1+JBSTRM(JBEL,JBSTAC)).EQ.0) THEN
                 JB = JBSTRM(JBEL,JBSTAC)
                 JEL1 = JBEL
                 GOTO 231
                END IF
  230        CONTINUE
  231        CONTINUE
             SIGNB = DBLE((-1)**(JEL1+IEL1))
           END IF
*
*. Two pairs of differing alpha electrons
*
           IF(NADIF .EQ. 2 ) THEN
             IDIFF = 0
             DO 320 IAEL = 1,NAEL
               IF(jWORK(KLJAE-1+LIA(IAEL)          ).EQ.0) THEN
                 IF( IDIFF .EQ. 0 ) THEN
                   IDIFF = 1
                   I1 = LIA(IAEL)
                   IPERM = IAEL
                 ELSE
                   I2 = LIA(IAEL)
                   IPERM = IAEL + IPERM
                   GOTO 321
                 END IF
               END IF
  320        CONTINUE
  321        CONTINUE
*
             JDIFF = 0
             DO 330 JAEL = 1,NAEL
               IF(jWORK(KLIAE-1+JASTRM(JAEL,JASTAC)).EQ.0) THEN
                 IF( JDIFF .EQ. 0 ) THEN
                   JDIFF = 1
                   J1 = JASTRM(JAEL,JASTAC)
                   JPERM = JAEL
                 ELSE
                   J2 = JASTRM(JAEL,JASTAC)
                   JPERM = JAEL + JPERM
                   GOTO 331
                 END IF
               END IF
  330        CONTINUE
  331        CONTINUE
             SIGN = DBLE((-1)**(IPERM+JPERM))
           END IF
*
*. Two pairs of differing beta electrons
*
           IF(NBDIF .EQ. 2 ) THEN
             IDIFF = 0
             DO 420 IBEL = 1,NBEL
               IF(jWORK(KLJBE-1+LIB(IBEL)          ).EQ.0) THEN
                 IF( IDIFF .EQ. 0 ) THEN
                   IDIFF = 1
                   I1 = LIB(IBEL)
                   IPERM = IBEL
                 ELSE
                   I2 = LIB(IBEL)
                   IPERM = IBEL + IPERM
                   GOTO 421
                  END IF
               END IF
  420        CONTINUE
  421        CONTINUE
*
             JDIFF = 0
             DO 430 JBEL = 1,NBEL
               IF(jWORK(KLIBE-1+JBSTRM(JBEL,JBSTAC)).EQ.0) THEN
                 IF( JDIFF .EQ. 0 ) THEN
                   JDIFF = 1
                   J1 = JBSTRM(JBEL,JBSTAC)
                   JPERM = JBEL
                 ELSE
                   J2 = JBSTRM(JBEL,JBSTAC)
                   JPERM = JBEL + JPERM
                   GOTO 431
                 END IF
               END IF
  430        CONTINUE
  431        CONTINUE
             SIGN = DBLE((-1)**(IPERM+JPERM))
           END IF
*
* =======================
* Value of matrix element
* =======================
*
        IF( NADIF .EQ. 2 .OR. NBDIF .EQ. 2 ) THEN
* 2 differences in alpha or beta strings
          NDIF2 = NDIF2 + 1
* SIGN * (I1 J1 ! I2 J2 ) - ( I1 J2 ! I2 J1 )
          XVAL = SIGN*( GTIJKL_MCLR(I1,J1,I2,J2)-
     &                  GTIJKL_MCLR(I1,J2,I2,J1) )
        ELSE IF( NADIF .EQ. 1 .AND. NBDIF .EQ. 1 ) THEN
*. 1 difference in alpha strings and one difference in beta string
          NDIF2 = NDIF2 + 1
* SIGN * (IA JA ! IB JB )
          XVAL = SIGNA*SIGNB* GTIJKL_MCLR(IA,JA,IB,JB)
* 1 differences in alpha or beta strings
        ELSE IF( NADIF .EQ. 1 .AND. NBDIF .EQ. 0 .OR.
     &           NADIF .EQ. 0 .AND. NBDIF .EQ. 1 )THEN
          NDIF1 = NDIF1 + 1
* SIGN *
*(  H(I1 J1 ) +
*  (SUM OVER ORBITALS OF BOTH      SPIN TYPES  ( I1 J1 ! JORB JORB )
* -(SUM OVER ORBITALS OF DIFFERING SPIN TYPE   ( I1 JORB ! JORB J1 ) )
          IF( NADIF .EQ. 1 ) THEN
            I1 = IA
            J1 = JA
            SIGN = SIGNA
          ELSE
            I1 = IB
            J1 = JB
            SIGN = SIGNB
          END IF
*
          XVAL = GETH1I_MCLR(I1,J1)
          DO 520 JAEL = 1, NAEL
            JORB = JASTRM(JAEL,JASTAC)
            XVAL = XVAL + GTIJKL_MCLR(I1,J1,JORB,JORB)
  520     CONTINUE
          DO 521 JBEL = 1, NBEL
            JORB = JBSTRM(JBEL,JBSTAC)
            XVAL = XVAL + GTIJKL_MCLR(I1,J1,JORB,JORB)
  521     CONTINUE
          IF( NADIF .EQ. 1 ) THEN
            DO 522 JAEL = 1, NAEL
              JORB = JASTRM(JAEL,JASTAC)
              XVAL = XVAL - GTIJKL_MCLR(I1,JORB,JORB,J1)
  522       CONTINUE
          ELSE
            DO 523 JBEL = 1, NBEL
              JORB = JBSTRM(JBEL,JBSTAC)
              XVAL = XVAL - GTIJKL_MCLR(I1,JORB,JORB,J1)
  523       CONTINUE
          END IF
          XVAL = XVAL * SIGN
        ELSE IF( NADIF .EQ. 0 .AND. NBDIF .EQ. 0 ) THEN
*. Diagonal elements
          NDIF0 = NDIF0 + 1
          XVAL = ECORE
          DO 650 IAB = 1, 2
            IF(IAB .EQ. 1 ) THEN
              NIABEL = NAEL
            ELSE
              NIABEL = NBEL
            END IF
            DO 640 JAB = 1, 2
              IF(JAB .EQ. 1 ) THEN
                NJABEL = NAEL
              ELSE
                NJABEL = NBEL
              END IF
              DO 630 IEL = 1, NIABEL
                IF( IAB .EQ. 1 ) THEN
                  IORB = LIA(IEL)
                ELSE
                  IORB = LIB(IEL)
                END IF
                IF(IAB .EQ. JAB ) XVAL = XVAL + GETH1I_MCLR(IORB,IORB)
                DO 620 JEL = 1, NJABEL
                  IF( JAB .EQ. 1 ) THEN
                    JORB = LIA(JEL)
                  ELSE
                    JORB = LIB(JEL)
                  END IF
                  XVAL = XVAL + 0.5D0*GTIJKL_MCLR(IORB,IORB,JORB,JORB)
*. test

                  IF( IAB . EQ. JAB )
     &            XVAL = XVAL - 0.5D0*GTIJKL_MCLR(IORB,JORB,JORB,IORB)
*. test
          FAC = GTIJKL_MCLR(IORB,JORB,JORB,IORB)
  620           CONTINUE
  630         CONTINUE
  640       CONTINUE
  650     CONTINUE
        END IF

        IF( ISYM .EQ. 0 ) THEN
          HAMIL((JDET-1)*NIDET+IDET) =
     &    HAMIL((JDET-1)*NIDET+IDET) + CONST * XVAL
        ELSE
          HAMIL((IDET-1)*IDET/2 + JDET ) =
     &    HAMIL((IDET-1)*IDET/2 + JDET ) + CONST * XVAL
        END IF
  898 CONTINUE
  899 CONTINUE
  900 CONTINUE
 1000 CONTINUE

      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(LWORK)
        CALL Unused_integer(NINOB)
        CALL Unused_integer(IPRT)
      END IF
      END
