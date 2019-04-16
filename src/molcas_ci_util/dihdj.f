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
* Copyright (C) 1989,2003, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE DIHDJ_MOLCAS(IASTR,IBSTR,NIDET,JASTR,JBSTR,NJDET,
     &     NAEL,NBEL,
     &     IWORK,NORB,ONEBOD,HAMIL,ISYM,NINOB,ECORE,ICOMBI,PSIGN,
     &     IPRINT,TUVX,ExFac,IREOTS)
C
C A SET OF DETERMINANTS IA DEFINED BY ALPHA AND BETASTRINGS
C IASTR,IBSTR AND ANOTHER SET OF DETERMINATS DEFINED BY STRINGS
C JASTR AND JBSTR ARE GIVEN . OBTAIN CORRESPONDING HAMILTONIAN MATRIX
C
C IF ICOMBI .NE. 0 COMBINATIONS ARE USED FOR ALPHA AND BETA STRING
C THAT DIFFERS :
C   1/SQRT(2) * ( |I1A I2B| + PSIGN * |I2A I1B| )
C
C IF ISYM .EQ. 0 FULL HAMILTONIAN IS CONSTRUCTED
C IF ISYM .NE. 0 LOWER HALF OF HAMILTONIAN IS CONSTRUCTED
C
C JEPPE OLSEN JANUARY 1989
*             IREOTS added for new Molcas compatibility, August 2003
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION JASTR(NAEL,*),JBSTR(NBEL,*)
C
      DIMENSION IWORK(*), HAMIL(*), ONEBOD(NORB,NORB)
      DIMENSION TUVX(*), IREOTS(*)
C
      NTEST = 0
C  Initialization
      KLIAB = 0
      IAEQIB = 0
      JAEQJB = 0
      CONST  = 0.0D0
      IEL1   = 0
      JEL1   = 0
      IPERM  = 0
      JPERM  = 0
      SIGN   = 0.0D0
      SIGNA  = 0.0D0
      SIGNB  = 0.0D0
      XVAL   = 0.0D0
      IA     = 0
      IB     = 0
      JA     = 0
      JB     = 0
      I1     = 0
      I2     = 0
      J1     = 0
      J2     = 0
C SCATCH SPACE
C
C .. 1 : EXPANSION OF ALPHA AND BETA STRINGS OF TYPE I
C
      KLFREE = 1
      KLIAE  = KLFREE
      KLFREE = KLIAE + NORB
      KLIBE  = KLFREE
      KLFREE = KLIBE + NORB
C
      KLJAE = KLFREE
      KLFREE = KLJAE + NORB
      KLJBE = KLFREE
      KLFREE = KLJBE + NORB
      IF( ICOMBI .NE. 0 ) THEN
        KLIAB  = KLFREE
        KLFREE = KLFREE + NIDET
      END IF
C
      IF( ICOMBI .NE. 0 ) THEN
C SET UP ARRAY COMBARING ALPHA AND BETA STRINGS IN IDET LIST
        DO 13 IDET = 1, NIDET
          IAEQIB = 1
          DO 14 IEL = 1, NAEL
            IF(IASTR(IEL,IDET) .NE. IBSTR(IEL,IDET))IAEQIB = 0
14        CONTINUE
          IWORK(KLIAB-1+IDET) = IAEQIB
13      CONTINUE
      END IF
C
      IF( ISYM .EQ. 0 ) THEN
        LHAMIL = NIDET*NJDET
      ELSE
        LHAMIL = NIDET*(NIDET+1) / 2
      END IF
      CALL DCOPY_(LHAMIL,[0.0D0],0,HAMIL,1)
C
      NTERMS= 0
      NDIF0 = 0
      NDIF1 = 0
      NDIF2 = 0
C.. LOOP OVER J DETERMINANTS
C
      DO 1000 JDET = 1,NJDET
C
C EXPAND JDET
        CALL ICOPY(NORB,[0],0,IWORK(KLJAE),1)
        CALL ICOPY(NORB,[0],0,IWORK(KLJBE),1)
C
        IF( ICOMBI .NE. 0 ) THEN
          JAEQJB = 1
          DO 32 IEL = 1, NAEL
            IF(JASTR(IEL,JDET) .NE. JBSTR(IEL,JDET))JAEQJB = 0
32        CONTINUE
        END IF
C
        DO 40 IAEL = 1, NAEL
          IWORK(KLJAE-1+JASTR(IAEL,JDET) ) = 1
40      CONTINUE
C
        DO 50 IBEL = 1, NBEL
          IWORK(KLJBE-1+JBSTR(IBEL,JDET) ) = 1
50      CONTINUE
C
        IF( NTEST .GE. 10 ) THEN
          WRITE(6,*) ' LOOP 1000 JDET =  ',JDET
          WRITE(6,*) ' JASTR AND JBSTR '
          CALL IWRTMA(JASTR(1,JDET),1,NAEL,1,NAEL)
          CALL IWRTMA(JBSTR(1,JDET),1,NBEL,1,NBEL)
          WRITE(6,*) ' EXPANDED ALPHA AND BETA STRING '
          CALL IWRTMA(IWORK(KLJAE),1,NORB,1,NORB)
          CALL IWRTMA(IWORK(KLJBE),1,NORB,1,NORB)
        END IF
C
        IF( ISYM .EQ. 0 ) THEN
          MINI = 1
        ELSE
          MINI = JDET
        END IF
        DO 900 IDET = MINI, NIDET
        IF( ICOMBI .EQ. 0 ) THEN
            NLOOP = 1
        ELSE
          IAEQIB = IWORK(KLIAB-1+IDET)
          IF(IAEQIB+JAEQJB .EQ. 0 ) THEN
            NLOOP = 2
          ELSE
            NLOOP = 1
          END IF
        END IF
C
        DO 899 ILOOP = 1, NLOOP
         NTERMS = NTERMS + 1
C
C.. COMPARE DETERMINANTS
C
C SWAP IA AND IB FOR SECOND PART OF COMBINATIONS
        If ( ILOOP.eq.2 ) then
          Do iii = 1,NAEL
            jjj = IASTR(iii,IDET)
            IASTR(iii,IDET) = IBSTR(iii,IDET)
            IBSTR(iii,IDET) = jjj
          End Do
        End If
C
        NACM = 0
        DO 61 IAEL = 1, NAEL
          NACM = NACM + IWORK(KLJAE-1+IASTR(IAEL,IDET))
61      CONTINUE
C
        NBCM = 0
        DO 62 IBEL = 1, NBEL
          NBCM = NBCM + IWORK(KLJBE-1+IBSTR(IBEL,IDET))
62      CONTINUE
C
        NADIF = NAEL-NACM
        NBDIF = NBEL-NBCM
        IF( NTEST .GE. 10 ) THEN
          WRITE(6,*) '  LOOP 900 IDET ',IDET
          WRITE(6,*) ' COMPARISON , NADIF , NBDIF ', NADIF,NBDIF
        END IF
C
        IF(NADIF+NBDIF .GT. 2 ) GOTO 899
C
C
C  FACTOR FOR COMBINATIONS
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
C
C.. FIND DIFFERING ORBITALS AND SIGN FOR PERMUTATION
C
C EXPAND IDET
        CALL ICOPY(NORB,[0],0,IWORK(KLIAE),1)
        CALL ICOPY(NORB,[0],0,IWORK(KLIBE),1)
C
          DO 42 IAEL = 1, NAEL
            IWORK(KLIAE-1+IASTR(IAEL,IDET) ) = 1
42        CONTINUE
C
          DO 52 IBEL = 1, NBEL
            IWORK(KLIBE-1+IBSTR(IBEL,IDET) ) = 1
52        CONTINUE
C
        IF(NADIF .EQ. 1 ) THEN
          DO 120 IAEL = 1,NAEL
            IF(IWORK(KLJAE-1+IASTR(IAEL,IDET)).EQ.0) THEN
              IA = IASTR(IAEL,IDET)
              IEL1 = IAEL
              GOTO 121
             END IF
120       CONTINUE
121       CONTINUE
C
          DO 130 JAEL = 1,NAEL
            IF(IWORK(KLIAE-1+JASTR(JAEL,JDET)).EQ.0) THEN
              JA = JASTR(JAEL,JDET)
              JEL1 = JAEL
              GOTO 131
             END IF
130       CONTINUE
131       CONTINUE
          SIGNA = DBLE((-1)**(JEL1+IEL1))
        END IF
        IF(NBDIF .EQ. 1 ) THEN
          DO 220 IBEL = 1,NBEL
            IF(IWORK(KLJBE-1+IBSTR(IBEL,IDET)).EQ.0) THEN
              IB = IBSTR(IBEL,IDET)
              IEL1 = IBEL
              GOTO 221
             END IF
220       CONTINUE
221       CONTINUE
C
          DO 230 JBEL = 1,NBEL
            IF(IWORK(KLIBE-1+JBSTR(JBEL,JDET)).EQ.0) THEN
              JB = JBSTR(JBEL,JDET)
              JEL1 = JBEL
              GOTO 231
             END IF
230       CONTINUE
231       CONTINUE
          SIGNB = DBLE((-1)**(JEL1+IEL1))
        END IF
        IF(NADIF .EQ. 2 ) THEN
          IDIFF = 0
          DO 320 IAEL = 1,NAEL
            IF(IWORK(KLJAE-1+IASTR(IAEL,IDET)).EQ.0) THEN
              IF( IDIFF .EQ. 0 ) THEN
                IDIFF = 1
                I1 = IASTR(IAEL,IDET)
                IPERM = IAEL
              ELSE
                I2 = IASTR(IAEL,IDET)
                IPERM = IAEL + IPERM
                GOTO 321
              END IF
            END IF
320       CONTINUE
321       CONTINUE
C
          JDIFF = 0
          DO 330 JAEL = 1,NAEL
            IF(IWORK(KLIAE-1+JASTR(JAEL,JDET)).EQ.0) THEN
              IF( JDIFF .EQ. 0 ) THEN
                JDIFF = 1
                J1 = JASTR(JAEL,JDET)
                JPERM = JAEL
              ELSE
                J2 = JASTR(JAEL,JDET)
                JPERM = JAEL + JPERM
                GOTO 331
              END IF
            END IF
330       CONTINUE
331       CONTINUE
          SIGN = DBLE((-1)**(IPERM+JPERM))
        END IF
C
        IF(NBDIF .EQ. 2 ) THEN
          IDIFF = 0
          DO 420 IBEL = 1,NBEL
            IF(IWORK(KLJBE-1+IBSTR(IBEL,IDET)).EQ.0) THEN
              IF( IDIFF .EQ. 0 ) THEN
                IDIFF = 1
                I1 = IBSTR(IBEL,IDET)
                IPERM = IBEL
              ELSE
                I2 = IBSTR(IBEL,IDET)
                IPERM = IBEL + IPERM
                GOTO 421
               END IF
            END IF
420       CONTINUE
421       CONTINUE
C
          JDIFF = 0
          DO 430 JBEL = 1,NBEL
            IF(IWORK(KLIBE-1+JBSTR(JBEL,JDET)).EQ.0) THEN
              IF( JDIFF .EQ. 0 ) THEN
                JDIFF = 1
                J1 = JBSTR(JBEL,JDET)
                JPERM = JBEL
              ELSE
                J2 = JBSTR(JBEL,JDET)
                JPERM = JBEL + JPERM
                GOTO 431
              END IF
            END IF
430       CONTINUE
431       CONTINUE
          SIGN = DBLE((-1)**(IPERM+JPERM))
        END IF
C
C OBTAIN VALUE OF HAMILTONIAN ELEMENT
C
        IF( NADIF .EQ. 2 .OR. NBDIF .EQ. 2 ) THEN
          NDIF2 = NDIF2 + 1
C SIGN * (I1 J1 | I2 J2 ) - ( I1 J2 | I2 J1 )
          I1 = I1 + NINOB
          I2 = I2 + NINOB
          J1 = J1 + NINOB
          J2 = J2 + NINOB
*. Well, there is no reordering in integrals any more. so
          I1_REO = IREOTS(I1)
          J1_REO = IREOTS(J1)
          I2_REO = IREOTS(I2)
          J2_REO = IREOTS(J2)
          XVAL = SIGN*( GETH2A(I1_REO,J1_REO,I2_REO,J2_REO,TUVX)
     *                 -GETH2A(I1_REO,J2_REO,I2_REO,J1_REO,TUVX) )
        ELSE IF( NADIF .EQ. 1 .AND. NBDIF .EQ. 1 ) THEN
          NDIF2 = NDIF2 + 1
C SIGN * (IA JA | IB JB )
          IA = IA + NINOB
          IB = IB + NINOB
          JA = JA + NINOB
          JB = JB + NINOB
*
          IA_REO = IREOTS(IA)
          IB_REO = IREOTS(IB)
          JA_REO = IREOTS(JA)
          JB_REO = IREOTS(JB)
          XVAL = SIGNA*SIGNB* GETH2A(IA_REO,JA_REO,IB_REO,JB_REO,TUVX)
        ELSE IF( NADIF .EQ. 1 .AND. NBDIF .EQ. 0 .OR.
     *           NADIF .EQ. 0 .AND. NBDIF .EQ. 1 )THEN
          NDIF1 = NDIF1 + 1
C SIGN *
C(  H(I1 J1 ) +
C  (SUM OVER ORBITALS OF BOTH      SPIN TYPES  ( I1 J1 | JORB JORB )
C -(SUM OVER ORBITALS OF DIFFERING SPIN TYPE   ( I1 JORB | JORB J1 ) )
          IF( NADIF .EQ. 1 ) THEN
            I1 = IA + NINOB
            J1 = JA + NINOB
            SIGN = SIGNA
          ELSE
            I1 = IB + NINOB
            J1 = JB + NINOB
            SIGN = SIGNB
          END IF
C
          I1_REO = IREOTS(I1)
          J1_REO = IREOTS(J1)
          XVAL = ONEBOD(IREOTS(I1-NINOB),IREOTS(J1-NINOB))
          DO 520 JAEL = 1, NAEL
            JORB = JASTR(JAEL,JDET)+NINOB
            JORB_REO = IREOTS(JORB)
            XVAL = XVAL + GETH2A(I1_REO,J1_REO,JORB_REO,JORB_REO,TUVX)
520       CONTINUE
          DO 521 JBEL = 1, NBEL
            JORB = JBSTR(JBEL,JDET)+NINOB
            JORB_REO = IREOTS(JORB)
            XVAL = XVAL + GETH2A(I1_REO,J1_REO,JORB_REO,JORB_REO,TUVX)
521       CONTINUE
          IF( NADIF .EQ. 1 ) THEN
            DO 522 JAEL = 1, NAEL
              JORB = JASTR(JAEL,JDET)+NINOB
              JORB_REO = IREOTS(JORB)
              XVAL = XVAL-GETH2A(I1_REO,JORB_REO,JORB_REO,J1_REO,TUVX)
522         CONTINUE
          ELSE
            DO 523 JBEL = 1, NBEL
              JORB = JBSTR(JBEL,JDET)+NINOB
              JORB_REO = IREOTS(JORB)
              XVAL = XVAL-GETH2A(I1_REO,JORB_REO,JORB_REO,J1_REO,TUVX)
523         CONTINUE
          END IF
          XVAL = XVAL * SIGN
        ELSE IF( NADIF .EQ. 0 .AND. NBDIF .EQ. 0 ) THEN
          NDIF0 = NDIF0 + 1
C SUM(I,J OF JDET) H(I,J) + (I I | J J ) - (I J | J I )
C
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
                  IORB = IASTR(IEL,IDET)
                ELSE
                  IORB = IBSTR(IEL,IDET)
                END IF
                IORB_REO = IREOTS(IORB)
                IF(IAB .EQ. JAB ) XVAL = XVAL+ONEBOD(IORB_REO,IORB_REO)
                IORB = IORB + NINOB
                DO 620 JEL = 1, NJABEL
                  IF( JAB .EQ. 1 ) THEN
                    JORB = IASTR(JEL,IDET)+NINOB
                  ELSE
                    JORB = IBSTR(JEL,IDET)+NINOB
                  END IF
                  JORB_REO = IREOTS(JORB)
                  XVAL = XVAL
     &          + 0.5D0*GETH2A(IORB_REO,IORB_REO,JORB_REO,JORB_REO,TUVX)
                  IF ( IAB . EQ. JAB ) XVAL =
     &               XVAL - ExFac * 0.5D0 *
     &               GETH2A(IORB_REO,JORB_REO,JORB_REO,IORB_REO,TUVX)
620             CONTINUE
630           CONTINUE
640         CONTINUE
650       CONTINUE
        END IF
C
        IF( ISYM .EQ. 0 ) THEN
          HAMIL((JDET-1)*NIDET+IDET) =
     *    HAMIL((JDET-1)*NIDET+IDET) + CONST * XVAL
        ELSE
          HAMIL((IDET-1)*IDET/2 + JDET ) =
     *    HAMIL((IDET-1)*IDET/2 + JDET ) + CONST * XVAL
        END IF
C RESTORE ORDER |||
        If ( ILOOP.eq.2 ) then
          Do iii = 1,NAEL
            jjj = IASTR(iii,IDET)
            IASTR(iii,IDET) = IBSTR(iii,IDET)
            IBSTR(iii,IDET) = jjj
          End Do
        End If
899   CONTINUE
900   CONTINUE
1000  CONTINUE
C
      IF( IPRINT .GE. 2 ) THEN
        WRITE(6,*) '  HAMILTONIAN MATRIX '
        IF( ISYM .EQ. 0 ) THEN
          CALL WRTMAT(HAMIL,NIDET,NJDET,NIDET,NJDET)
        ELSE
          CALL PRSYM(HAMIL,NIDET)
        END IF
      END IF
C
      RETURN
      END
