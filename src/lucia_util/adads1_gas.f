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
      SUBROUTINE ADADS1_GAS(     NK,     I1,   XI1S,    LI1,   IORB,
     &                        NIORB,   JORB,  NJORB,   KSTR,   NKEL,
     &                        NKSTR,   IREO,     IZ,  NOCOB,   KMAX,
     &                         KMIN,   IEND, SCLFAC)
*
* Obtain I1(KSTR) = +/- A+ IORB A+ JORB !KSTR>
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
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER KSTR(NKEL,NKSTR)
      INTEGER IREO(*), IZ(NOCOB,*)
*.Output
      INTEGER I1(LI1,*)
      DIMENSION XI1S(LI1,*)
*
      NTEST = 000
      IF(NTEST.NE.0) THEN
       WRITE(6,*) ' ==================== '
       WRITE(6,*) ' ADADS1_GAS speaking '
       WRITE(6,*) ' ==================== '
       WRITE(6,*) ' IORB,NIORB ', IORB,NIORB
       WRITE(6,*) ' JORB,NJORB ', JORB,NJORB
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
      DO KKSTR = KMIN,KEND
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Occupation of string ', KKSTR
         CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
       END IF
*. Loop over electrons after which JORB can be added
*PAM2009 Added variable ODDJEL to replace ''DBLE((-1)**JEL)''
       ODDJEL=-1.0D0
       DO JEL = 0, NKEL
         ODDJEL=-ODDJEL
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
     &    WRITE(6,*) ' JEL JORB1 JORB2 ',JEL,JORB1,JORB2
*
         IF(JEL.GT.0.AND.JORB1.GE.JORBMIN.AND.
     &                   JORB1.LE.JORBMAX) THEN
*. vanishing for any IORB
           IJOFF = (JORB1-JORBMIN)*NIORB
           DO IIORB = 1, NIORB
             IJ = IJOFF + IIORB
             if(ij.gt.nij) then
               write(6,*) ' ij .gt. nij '
               write(6,*) ' JORB1 IIORB' , JORB1,IIORB
               write(6,*) ' ijoff ', ijoff
*               stop
        CALL SYSABENDMSG('lucia_util/adads1_gas','Internal error',' ')
             end if
             I1(KKSTR-KMIN+1,IJ) = 0
             XI1S(KKSTR-KMIN+1,IJ) = 0.0D0
           END DO
         END IF
*
         IF(JORB1.LT.JORBMAX.AND.JORB2.GT.JORBMIN) THEN
*. Orbitals JORB1+1 - JORB2-1 can be added after electron JEL
*PAM2009           SIGNJ = DBLE((-1)**JEL) * SCLFAC
           SIGNJ=ODDJEL*SCLFAC
*. reverse lexical number of the first JEL ELECTRONS
           ILEX0 = 1
           DO JJEL = 1, JEL
             ILEX0 = ILEX0 + IZ(KSTR(JJEL,KKSTR),JJEL)
           END DO
           DO JJORB = JORB1+1, JORB2-1
* And electron JEL + 1
             ILEX1 = ILEX0 + IZ(JJORB,JEL+1)
*. Add electron IORB
*PAM2009 Added ODDIEL to replace ''(-1.0)**IEL''
             ODDIEL=-ODDJEL
             DO IEL = JEL, NKEL
               ODDIEL=-ODDIEL
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
               IF(NTEST.GE.1000)
     &          WRITE(6,*) ' IEL IORB1 IORB2 ',IEL,IORB1,IORB2
               IF(IEL.GT.JEL.AND.IORB1.GE.IORBMIN.AND.
     &                           IORB1.LE.IORBMAX) THEN
                 IJ = (JJORB-JORBMIN)*NIORB+IORB1-IORBMIN+1
             if(ij.gt.nij) then
               write(6,*) ' ij .gt. nij '
               write(6,*) ' JJORB IORB1' , JJORB,IORB1
               write(6,*) ' ijoff ', ijoff
*               stop
               CALL SYSABENDMSG('lucia_util/adads1_gas',
     &                          'Internal error',' ')
             end if
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
*PAM2009                 SIGNIJ =  SIGNJ*(-1.0D0) ** (IEL+1)
                 SIGNIJ =  SIGNJ*(-ODDIEL)
                 DO IIORB = IORB1+1, IORB2-1
                   IJ = IJOFF + IIORB - IORBMIN + 1
                   IF(IJ.LE.0.OR.IJ.GT.NIJ) THEN
                     WRITE(6,*) ' PROBLEMO ADADS1 : IJ : ', IJ
                     WRITE(6,*) ' IJOFF IORBMIN ', IJOFF,IORBMIN
                     WRITE(6,*) ' IIORB JJORB ', IIORB,JJORB
*                     stop
                     CALL SYSABENDMSG('lucia_util/adads1_gas',
     &                          'Internal error',' ')
                     NTEST = 1000
                   END IF
                   ILEX = ILEX2 + IZ(IIORB,IEL+2)
                   IACT = IREO(ILEX)
                   IF(NTEST.GE.1000)
     &             WRITE(6,*) ' IIORB JJORB ', IIORB,JJORB
                   IF(NTEST.GE.1000)
     &             WRITE(6,*) ' IJ ILEX,IACT',IJ, ILEX,IACT
                   IF(NTEST.GE.1000)
     &             WRITE(6,*) ' ILEX0 ILEX1 ILEX2 ILEX ',
     &                          ILEX0,ILEX1,ILEX2,ILEX
                   I1(KKSTR-KMIN+1,IJ) = IACT
                   XI1S(KKSTR-KMIN+1,IJ) = SIGNIJ
                   IF(IJ.LT.0) THEN
                     WRITE(6,*) ' NEGATIVE IJ in ADADS1 '
*                     STOP ' NEGATIVE IJ in ADADS1 '
                     CALL SYSABENDMSG('lucia_util/adads1_gas',
     &                          'Internal error',' ')
                   END IF
                 END DO
               END IF
             END DO
           END DO
         END IF
       END DO
      END DO
*
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Output from ADADST1_GAS '
        WRITE(6,*) ' ===================== '
        WRITE(6,*) ' Number of K strings accessed ', NK
        IF(NK.NE.0) THEN
          IJ = 0
          DO  JJORB = JORB,JORB+NJORB-1
            DO  IIORB = IORB, IORB + NIORB - 1
              IJ = IJ + 1
C?            WRITE(6,*) ' IJ = ', IJ
C?            IF(IIORB.GT.JJORB) THEN
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
