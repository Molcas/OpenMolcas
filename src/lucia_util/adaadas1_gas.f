************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE ADAADAS1_GAS(     NK,     I1,   XI1S,    LI1,   IORB,
     &                          NIORB,    IAC,   JORB,  NJORB,    JAC,
     &                           KSTR,   NKEL,  NKSTR,   IREO,     IZ,
     &                          NOCOB,   KMAX,   KMIN,   IEND, SCLFAC,
     &                          NSTRI)
*
*
* Obtain I1(KSTR) = +/- a+/a  IORB a+/a JORB !KSTR>
* Only orbital pairs IOB .gt. JOB are included
*
* KSTR is restricted to strings with relative numbers in the
* range KMAX to KMIN
* =====
* Input
* =====
* IORB : First I orbital to be added
* NIORB : Number of orbitals to be added : IORB to IORB-1+NIORB
*        are used. They must all be in the same TS group
* JORB : First J orbital to be added
* LORB : Number of orbitals to be added : JORB to JORB-1+NJORB
*        are used. They must all be in the same TS group
* KMAX : Largest allowed relative number for K strings
* KMIN : Smallest allowed relative number for K strings
*
* ======
* Output
* ======
*
* NK      : Number of K strings
* I1(KSTR,JORB) : ne. 0 =>  a+IORB a+JORB !KSTR> = +/-!ISTR>
* XI1S(KSTR,JORB) : above +/-
*          : eq. 0    a + JORB !KSTR> = 0
* Offset is KMIN
*
* L.R. Jan 20, 1998
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER KSTR(NKEL,NKSTR)
      INTEGER IREO(*), IZ(NOCOB,*)
*.Output
      INTEGER I1(LI1,*)
      DIMENSION XI1S(LI1,*)
*. Local scratch, atmost 1000 orbitals in a given TS block)
      DIMENSION ISCR(1000)
*PAM2009 Array of values, replacing expressions such as ''DBLE((-1)**INT8)'':
      DIMENSION SGNARR(0:63)
      DO I=0,62,2
       SGNARR(I)  = 1.0D0
       SGNARR(I+1)=-1.0D0
      END DO
*
* Some dummy initializations
      SIGN = 0.0D0

      NTEST = 00
      IF(NTEST.NE.0) THEN
       WRITE(6,*) ' ==================== '
       WRITE(6,*) ' ADADS1_GAS speaking '
       WRITE(6,*) ' ==================== '
       WRITE(6,*) ' IORB,NIORB,IAC ', IORB,NIORB,IAC
       WRITE(6,*) ' JORB,NJORB,JAC ', JORB,NJORB,JAC
*
C      WRITE(6,*) ' Kstrings in action (el,string) '
C      WRITE(6,*) ' ==============================='
C      CALL IWRTMA(KSTR,NKEL,NKSTR,NKEL,NKSTR)
       WRITE(6,*) ' 24 elements of reorder array'
       CALL IWRTMA(IREO,1,24,1,24)
*
      END IF
*
      IORBMIN = IORB
      IORBMAX = IORB + NIORB - 1
*
      JORBMIN = JORB
      JORBMAX = JORB + NJORB - 1
*
      NIJ = NIORB*NJORB
*
      KEND = MIN(NKSTR,KMAX)
      IF(KEND.LT.NKSTR) THEN
        IEND = 0
      ELSE
        IEND = 1
      END IF
      NK = KEND-KMIN+1
*
      IF(IAC.EQ.2.AND.JAC.EQ.2) THEN
*
* ==========================
* Creation- creation mapping
* ==========================
*
        DO KKSTR = KMIN,KEND
         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' Occupation of string ', KKSTR
           CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
         END IF
*. Loop over electrons after which JORB can be added
         DO JEL = 0, NKEL
*
           IF(JEL.EQ.0 ) THEN
             JORB1 = JORBMIN - 1
           ELSE
             JORB1 = MAX(JORBMIN-1,KSTR(JEL,KKSTR))
           END IF
           IF(JEL.EQ.NKEL) THEN
             JORB2 = JORBMAX + 1
           ELSE
             JORB2 = MIN(JORBMAX+1,KSTR(JEL+1,KKSTR))
           END IF
           IF(NTEST.GE.1000)
     &      WRITE(6,*) ' JEL JORB1 JORB2 ',JEL,JORB1,JORB2
*
           IF(JEL.GT.0.AND.JORB1.GE.JORBMIN.AND.
     &                     JORB1.LE.JORBMAX) THEN
*. vanishing for any IORB
             IJOFF = (JORB1-JORBMIN)*NIORB
             DO IIORB = 1, NIORB
               IJ = IJOFF + IIORB
               I1(KKSTR-KMIN+1,IJ) = 0
               XI1S(KKSTR-KMIN+1,IJ) = 0.0D0
             END DO
           END IF
*
           IF(JORB1.LT.JORBMAX.AND.JORB2.GT.JORBMIN) THEN
*. Orbitals JORB1+1 - JORB2-1 can be added after electron JEL
*PAM2009             SIGNJ = DBLE((-1)**JEL) * SCLFAC
             SIGNJ = SGNARR(JEL) * SCLFAC
*. reverse lexical number of the first JEL ELECTRONS
             ILEX0 = 1
             DO JJEL = 1, JEL
               ILEX0 = ILEX0 + IZ(KSTR(JJEL,KKSTR),JJEL)
             END DO
             DO JJORB = JORB1+1, JORB2-1
* And electron JEL + 1
               ILEX1 = ILEX0 + IZ(JJORB,JEL+1)
*. Add electron IORB
               DO IEL = JEL, NKEL
                 IF(IEL.EQ.JEL) THEN
                   IORB1 = MAX(JJORB,IORBMIN-1)
                 ELSE
                   IORB1 = MAX(IORBMIN-1,KSTR(IEL,KKSTR))
                 END IF
                   IF(IEL.EQ.NKEL) THEN
                   IORB2 = IORBMAX+1
                 ELSE
                   IORB2 = MIN(IORBMAX+1,KSTR(IEL+1,KKSTR))
                 END IF
                 IF(NTEST.GE.5000)
     &            WRITE(6,*) ' IEL IORB1 IORB2 ',IEL,IORB1,IORB2
                 IF(IEL.GT.JEL.AND.IORB1.GE.IORBMIN.AND.
     &                             IORB1.LE.IORBMAX) THEN
                   IJ = (JJORB-JORBMIN)*NIORB+IORB1-IORBMIN+1
                   I1(KKSTR-KMIN+1,IJ) = 0
                   XI1S(KKSTR-KMIN+1,IJ) = 0.0D0
                 END IF
                 IF(IORB1.LT.IORBMAX.AND.IORB2.GT.IORBMIN) THEN
*. Orbitals IORB1+1 - IORB2 -1 can be added after ELECTRON IEL in KSTR
*. Reverse lexical number of the first IEL+1 electrons
                   ILEX2 = ILEX1
                   DO IIEL = JEL+1,IEL
                     ILEX2 = ILEX2 + IZ(KSTR(IIEL,KKSTR),IIEL+1)
                   END DO
*. add terms for the last electrons
                   DO IIEL = IEL+1,NKEL
                     ILEX2 = ILEX2 + IZ(KSTR(IIEL,KKSTR),IIEL+2)
                   END DO
                   IJOFF = (JJORB-JORBMIN)*NIORB
*PAM2009           SIGNIJ =  SIGNJ*(-1.0D0) ** (IEL+1)
                   SIGNIJ =  SIGNJ*SGNARR(IEL+1)
                   DO IIORB = IORB1+1, IORB2-1
                     IJ = IJOFF + IIORB - IORBMIN + 1
                     ILEX = ILEX2 + IZ(IIORB,IEL+2)
                     IACT = IREO(ILEX)
                     IF(NTEST.GE.1000) THEN
                       WRITE(6,*) 'IIORB JJORB', IIORB,JJORB
                       WRITE(6,*) ' ILEX IACT ', ILEX,IACT
                     END IF
                     I1(KKSTR-KMIN+1,IJ) = IACT
                     XI1S(KKSTR-KMIN+1,IJ) = SIGNIJ
                   END DO
                 END IF
               END DO
             END DO
           END IF
         END DO
        END DO
      ELSE IF(IAC.EQ.1.AND.JAC.EQ.1) THEN
*
* ===========================================
* annihilation - annihilation mapping (a i a j)
* ===========================================
*
        DO KKSTR = KMIN,KEND
*. Active range for electrons
         IIELMIN = 0
         IIELMAX = 0
         JJELMIN = 0
         JJELMAX = 0
*
         DO KEL = 1, NKEL
          KKORB = KSTR(KEL,KKSTR)
          IF(IIELMIN.EQ.0.AND.KKORB.GE.IORBMIN)IIELMIN = KEL
          IF(JJELMIN.EQ.0.AND.KKORB.GE.JORBMIN)JJELMIN = KEL
          IF(KKORB.LE.IORBMAX) IIELMAX = KEL
          IF(KKORB.LE.JORBMAX) JJELMAX = KEL
         END DO
         IF(IIELMIN.EQ.0) IIELMIN = NKEL  + 1
         IF(JJELMIN.EQ.0) JJELMIN = NKEL  + 1

         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' Occupation of string ', KKSTR
           CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
         END IF
*. Loop over first electron to be removed
C        DO JEL = 1, NKEL
         DO JEL = JJELMIN,JJELMAX
           JJORB = KSTR(JEL,KKSTR)
*. Loop over second electron to be removed
C          DO IEL = JEL+1, NKEL
           DO IEL = MAX(JEL+1,IIELMIN),IIELMAX
             IIORB = KSTR(IEL,KKSTR)
             IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' IEL JEL IORB JORB ',
     &        IEL,JEL,IORB,JORB
             END IF
             IF(IIORB.GE.IORBMIN.AND.IIORB.LE.IORBMAX.AND.
     &          JJORB.GE.JORBMIN.AND.JJORB.LE.JORBMAX     )THEN
*PAM09               SIGN = DBLE((-1)**(IEL+JEL-1)) * SCLFAC
               SIGN = SGNARR(IEL+JEL+1) * SCLFAC
*. reverse lexical number of the double annihilated string
               ILEX = 1
               DO JJEL = 1, JEL-1
                 ILEX = ILEX + IZ(KSTR(JJEL,KKSTR),JJEL)
               END DO
               DO JJEL = JEL+1,IEL-1
                 ILEX = ILEX + IZ(KSTR(JJEL,KKSTR),JJEL-1)
               END DO
               DO JJEL = IEL+1, NKEL
                 ILEX = ILEX + IZ(KSTR(JJEL,KKSTR),JJEL-2)
               END DO
               IACT = IREO(ILEX)
               IF(IACT.LE.0.OR.IACT.GT.NSTRI) THEN
                 WRITE(6,*) ' IACT out of bounds, IACT =  ', IACT
*                 STOP       ' IACT out of bounds '
        CALL SYSABENDMSG('lucia_util/adaadas1_gas',
     &                    'Internal error',' ')
               END IF
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' ILEX and IACT ', ILEX, IACT
               END IF
*
               IJ = (JJORB-JORB)*NIORB + IIORB-IORB+1
               I1(KKSTR-KMIN+1,IJ) = IACT
               XI1S(KKSTR-KMIN+1,IJ) = SIGN
             END IF
*            ^ End if orbitals are in correct range
           END DO
*          ^ End of loop over IEL
         END DO
*        ^ End of loop over JEL
        END DO
*       ^ End of loop over Kstrings
*
      ELSE IF(IAC.EQ.2.AND.JAC.EQ.1) THEN
*
* ===================================
* Creation-annihilation map a+ i a j
* ===================================
*
C       DO KKSTR = 1, NKSTR
        DO KKSTR = KMIN,KEND
*. Indicate where a given orbital i should be added in KKSTR
         IZERO = 0
         CALL ISETVC(ISCR(IORBMIN),IZERO,NIORB)
         IIEL = 1
         DO IIORB = IORBMIN,IORBMAX
 2810      CONTINUE
           IF(IIEL.LE.NKEL) THEN
             IF(IIORB.LT.KSTR(IIEL,KKSTR)) THEN
               ISCR(IIORB)=IIEL
             ELSE IF (IIORB.EQ.KSTR(IIEL,KKSTR)) THEN
               ISCR(IIORB) = - IIEL
               IIEL = IIEL+1
             ELSE IF (IIORB.GT.KSTR(IIEL,KKSTR)) THEN
               IIEL = IIEL + 1
               GOTO 2810
             END IF
           ELSE IF(IIEL.EQ.NKEL+1) THEN
              ISCR(IIORB) = IIEL
           END IF
         END DO
         IF(NTEST.GE.10000) THEN
           WRITE(6,*) ' ISCR from IORBMIN array for KKSTR = ', KKSTR
           WRITE(6,*) ' IORBMIN, NIORB', IORBMIN,NIORB
           CALL IWRTMA(ISCR(IORBMIN),1,NIORB,1,NIORB)
         END IF
         DO JEL = 1, NKEL
           JJORB = KSTR(JEL,KKSTR)
           IF(JJORB.GE.JORBMIN.AND.JJORB.LE.JORBMAX)THEN
             DO IIORB = IORBMIN,IORBMAX
               IEL = ISCR(IIORB)
C?             write(6,*) ' JEL IEL JJORB IIORB',JEL,IEL,JJORB,IIORB
               IACT = 0
               IF(IEL.GT.0.AND.IIORB.GT.JJORB) THEN
*. New string is  a+1 ... a+ jel-1 a+jel+1 ..a+iel-1 a+iiorb a+iel+1 ...
*. Lexical number of new string
                 ILEX = 1
                 DO KEL = 1, JEL-1
                  ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                 END DO
                 DO KEL = JEL+1, IEL-1
                  ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL-1)
                 END DO
                 ILEX = ILEX + IZ(IIORB,IEL-1)
                 DO KEL = IEL, NKEL
                  ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                 END DO
                 IACT = IREO(ILEX)
                 IF(IACT.LE.0.OR.IACT.GT.NSTRI) THEN
                   WRITE(6,*) ' 1: IACT out of bounds, IACT =  ', IACT
                   WRITE(6,*) ' NSTRI = ', NSTRI
                   WRITE(6,*) 'IIORB,JJORB ',IIORB,JJORB
                   WRITE(6,*) ' Kstring : '
                   CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
                   WRITE(6,*) ' ILEX = ', ILEX
                   WRITE(6,*) 'IZ matrix'
                   CALL IWRTMA(IZ,NOCOB,NKEL,NOCOB,NKEL)
*                   STOP ' IACT out of bounds'
        CALL SYSABENDMSG('lucia_util/adaadas1_gas',
     &                    'Internal error',' ')
                 END IF
*PAM2009                 SIGN = DBLE((-1)**(IEL+JEL-1)) * SCLFAC
                 SIGN = SGNARR(IEL+JEL+1) * SCLFAC
               ELSE IF(IEL.GT.0 .AND. IIORB.LT.JJORB) THEN
*. New string is  a+1 ... a+ iel-1 a+ iiorb a+iel+1 ..a+jel-1 a+jel+1 ...
                 ILEX = 1
                 DO KEL = 1, IEL-1
                   ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                 END DO
                 ILEX = ILEX + IZ(IIORB,IEL)
                 DO KEL = IEL,JEL-1
                   ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL+1)
                 END DO
                 DO KEL = JEL + 1, NKEL
                   ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                 END DO
C?               write(6,*) ' ILEX =' , ILEX
                 IACT = IREO(ILEX)
                 IF(IACT.LE.0.OR.IACT.GT.NSTRI) THEN
                   WRITE(6,*) '2 IACT out of bounds, IACT =  ', IACT
*                   STOP       ' IACT out of bounds '
        CALL SYSABENDMSG('lucia_util/adaadas1_gas',
     &                    'Internal error',' ')
                 END IF
C?               write(6,*) ' IACT = ', IACT
*PAM2009                 SIGN = DBLE((-1)**(IEL+JEL  )) * SCLFAC
                 SIGN = SGNARR(IEL+JEL) * SCLFAC
               ELSE IF (IEL.LT.0. AND. IIORB .EQ. JJORB) THEN
*. Diagonal excitation
                 SIGN = SCLFAC
                 IACT = KKSTR
               END IF
               IF(IACT.NE.0) THEN
                 IJ = (JJORB-JORB)*NIORB + IIORB-IORB+1
                 I1(KKSTR-KMIN+1,IJ) = IACT
                 XI1S(KKSTR-KMIN+1,IJ) = SIGN
               END IF
             END DO
*            ^ End of loop over IIORB
           END IF
*          ^ End of  active cases
         END DO
*        ^ End of loop over electrons to be annihilated
        END DO
*       ^ End of loop over Kstrings
      ELSE IF(IAC.EQ.1.AND.JAC.EQ.2) THEN
*
* ======================================
* Annihilation-creation  map a i a+ j
* ======================================
*
*. Diagonal excitations ?
        IF(IORBMIN.EQ.JORBMIN) THEN
         IDIAG = 1
        ELSE
         IDIAG = 0
        END IF
C       DO KKSTR = 1, NKSTR
        DO KKSTR = KMIN,KEND
*. Indicate where a given orbital j should be added in KKSTR
         IZERO = 0
         CALL ISETVC(ISCR(JORBMIN),IZERO,NJORB)
         JJEL = 1
         DO JJORB = JORBMIN,JORBMAX
 0803      CONTINUE
           IF(JJEL.LE.NKEL) THEN
             IF(JJORB.LT.KSTR(JJEL,KKSTR)) THEN
               ISCR(JJORB)=JJEL
             ELSE IF (JJORB.EQ.KSTR(JJEL,KKSTR)) THEN
               ISCR(JJORB) = - JJEL
               JJEL = JJEL+1
             ELSE IF (JJORB.GT.KSTR(JJEL,KKSTR)) THEN
               JJEL = JJEL + 1
               GOTO 0803
             END IF
           ELSE IF(JJEL.EQ.NKEL+1)THEN
              ISCR(JJORB) = JJEL
           END IF
         END DO
         IF(NTEST.GE.10000) THEN
           WRITE(6,*) ' ISCR from JORBMIN array for KKSTR = ', KKSTR
           CALL IWRTMA(ISCR(JORBMIN),1,NJORB,1,NJORB)
         END IF
         DO IEL = 1, NKEL+IDIAG
*. IEL = NKEL + 1 will be used for excitations a+j aj
           IF(IEL.LE.NKEL) IIORB = KSTR(IEL,KKSTR)
           IF((IIORB.GE.IORBMIN.AND.IIORB.LE.IORBMAX).OR.
     &         IEL.EQ.NKEL+1 )THEN
             DO JJORB = JORBMIN,JORBMAX
               JEL = ISCR(JJORB)
C?             WRITE(6,*) ' IEL IIORB JEL JJORB ',
C?   &                      IEL,IIORB,JEL,JJORB
               IACT = 0
               IF(IEL.LE.NKEL) THEN
                 IF(JEL.GT.0.AND.JJORB.GT.IIORB) THEN
*. New string is  a+1 ... a+ iel-1 a+iel+1 ..a+jel-1 a+jjorb a+jel+1 ...
*. Lexical number of new string
                   ILEX = 1
                   DO KEL = 1, IEL-1
                    ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                   END DO
                   DO KEL = IEL+1, JEL-1
                    ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL-1)
                   END DO
                   ILEX = ILEX + IZ(JJORB,JEL-1)
                   DO KEL = JEL, NKEL
                    ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                   END DO
                   IACT = IREO(ILEX)
                   IF(IACT.LE.0.OR.IACT.GT.NSTRI) THEN
                     WRITE(6,*) '3 IACT out of bounds, IACT =  ', IACT
                     WRITE(6,*) ' ILEX,KKSTR ', ILEX, KKSTR
                     WRITE(6,*) ' occupation of KSTR '
                     CALL IWRTMA(KSTR,1,NKEL,1,NKEL)
                     WRITE(6,*) ' IEL JEL ', IEL,JEL
                     WRITE(6,*) ' IIORB,JJORB',IIORB,JJORB
*                     STOP       ' IACT out of bounds '
        CALL SYSABENDMSG('lucia_util/adaadas1_gas',
     &                    'Internal error',' ')
                   END IF
*PAM2009                   SIGN = DBLE((-1)**(IEL+JEL)) * SCLFAC
                   SIGN = SGNARR(IEL+JEL) * SCLFAC
                 ELSE IF(JEL.GT.0 .AND. JJORB.LT.IIORB) THEN
*. New string is  a+1 ... a+ jel-1 a+ jjorb a+jel+1 ..a+iel-1 a+iel+1 ...
                   ILEX = 1
                   DO KEL = 1, JEL-1
                     ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                   END DO
                   ILEX = ILEX + IZ(JJORB,JEL)
                   DO KEL = JEL,IEL-1
                     ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL+1)
                   END DO
                   DO KEL = IEL + 1, NKEL
                     ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                   END DO
                   IACT = IREO(ILEX)
*PAM2009                   SIGN = DBLE((-1)**(IEL+JEL-1)) * SCLFAC
                   SIGN = SGNARR(IEL+JEL+1) * SCLFAC
                   IF(IACT.LE.0.OR.IACT.GT.NSTRI) THEN
                     WRITE(6,*) '4 IACT out of bounds, IACT =  ', IACT
                     WRITE(6,*) ' NSTRI = ', NSTRI
                     WRITE(6,*) 'IIORB,JJORB ',IIORB,JJORB
                     WRITE(6,*) ' Kstring : '
                     CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
                     WRITE(6,*) ' ILEX = ', ILEX
                     WRITE(6,*) 'IZ matrix'
                     CALL IWRTMA(IZ,NOCOB,NKEL,NOCOB,NKEL)
*                     STOP       ' IACT out of bounds '
        CALL SYSABENDMSG('lucia_util/adaadas1_gas',
     &                    'Internal error',' ')
                   END IF
                 END IF
                 IF(IACT.NE.0) THEN
                   IJ = (JJORB-JORB)*NIORB + IIORB-IORB+1
                   I1(KKSTR-KMIN+1,IJ) = IACT
                   XI1S(KKSTR-KMIN+1,IJ) = SIGN
                 END IF
               ELSE IF(IEL.EQ.NKEL+1.AND.JEL.GT.0) THEN
*. Diagonal excitations aja+j
                 JJ = (JJORB-JORB)*NJORB + JJORB-JORB+1
                 I1(KKSTR-KMIN+1,JJ) = KKSTR
                 XI1S(KKSTR-KMIN+1,JJ) = SCLFAC
               END IF
             END DO
*            ^ End of loop over JJORB
           END IF
*          ^ End of  active cases
         END DO
*        ^ End of loop over electrons to be annihilated
        END DO
*       ^ End of loop over Kstrings
      END IF
*.    ^ End of types of creation mappings
*
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Output from ADADST1_GAS '
        WRITE(6,*) ' ===================== '
        WRITE(6,*) ' Number of K strings accessed ', NK
        IF(NK.NE.0) THEN
          IJ = 0
          DO  JJORB = JORB,JORB+NJORB-1
            JJORBR = JJORB-JORB+1
            DO  IIORB = IORB, IORB + NIORB - 1
              IJ = IJ + 1
C?            WRITE(6,*) ' IJ = ', IJ
C?            IF(IIORB.GT.JJORB) THEN
                IIORBR = IIORB - IORB + 1
                WRITE(6,*)
     &          ' Info for orbitals (iorb,jorb) ', IIORB,JJORB
                WRITE(6,*) ' Excited strings and sign '
                CALL IWRTMA(I1(1,IJ),1,NK,1,NK)
                CALL WRTMAT(XI1S(1,IJ),1,NK,1,NK)
C?            END IF
            END DO
          END DO
        END IF
      END IF
*
      RETURN
      END
