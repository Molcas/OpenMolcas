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
      SUBROUTINE ADADS1(NK,I1,XI1S,IOBSM,IOBTP,IOBOFF,NIOB,
     &                  JOBSM,JOBTP,JOBOFF,NJOB,IJORD,NKDIM,
     &                  ICLS,ISM,I2MAPO,I2MAPS,
     &                 I2MPF,L2MP,I2MPO,I2MPL,
     &                 I1MAPO,I1MAPS,
     &                 I1MPF,L1MP,I1MPO,I1MPL,
     &                 IEL1,IEL3,I2EL1,I2EL3,
     &                 ISSO,NSSO,I2SSO,N2SSO,NOCTP,N2OCTP,
     &                 NORB1,NORB2,NORB3,NORB,KMAX,KMIN,IEND)
*
* Obtain I1(KSTR) = +/- A+ IORB !KSTR>
*
* KSTR is restricted to strings with relative numbers in the
* range KMAX to KMIN
* =====
* Input
* =====
* ICLS,ISM : Class and symmetry of string with added electron
* I2MAPO,I2MAPS : map N-2 strings to N-1 strings
* I1MAPO,I1MAPS : map N-1 strings to  strings
* IEL1(3) : Number of electrons in RAS1(3) for I strings
* I2EL1(3) : Number of electrons in RAS1(3) for K strings
* ISSO : TS symmetry offset for I strings
* NSSO : Number of TS strings for I strings
* I2SSO : TS symmetry offset for K strings
* N2SSO : Number of TS strings for K strings
* NOCTP : Number of occupation types for I strings
* N2OCTP : Number of occupation types for K strings
* NORB1(2,3) : Number of RAS1(2,3) orbitals
* IORBSM : Orbital symmety array
* NORB : Number of active  orbitals
* KMAX : Largest allowed relative number for K strings
* KMIN : Smallest allowed relative number for K strings
*
* ======
* Output
* ======
*
* NK      : Number of K strings
* I1(KSTR) : ne. 0 => a + IORB a+ JORB !KSTR> = +/-!ISTR>
* XI1S(KSTR) : above +/-
*          : eq. 0    a + IORB !KSTR> = 0
* Offset is KMIN
* IEND : = 0 => end of N-2 strings has not been encountered
* IEND : = 1 => end of N-2 strings has     been encountered
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER IEL1(*),IEL3(*),I2EL1(*),I2EL3(*)
      INTEGER ISSO(NOCTP,*),NSSO(NOCTP,*)
      INTEGER I2SSO(N2OCTP,*),N2SSO(N2OCTP,*)
      INTEGER I2MAPO(*),I2MAPS(*)
      INTEGER I1MAPO(*),I1MAPS(*)
      INTEGER I2MPL(*),I2MPO(*)
*eaw
      Integer I1mpo(*),i1mpl(*)
*eaw
*.Output
      INTEGER I1(NKDIM,*)
      DIMENSION XI1S(NKDIM,*)
*
      NIJ = 0 ! dummy initialize
      iprstr=0
      NTEST =  IPRSTR !100
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' ================ '
        WRITE(6,*) ' Info from ADADS1 '
        WRITE(6,*) ' ================ '
        WRITE(6,*)
     &' Iobsm Iobtp Ioboff, Niob',IOBSM,IOBTP,IOBOFF,NIOB
        WRITE(6,*)
     &' Jobsm Jobtp Joboff, NJob',JOBSM,JOBTP,JOBOFF,NJOB
        WRITE(6,*) ' icls ism',ICLS,ISM
        WRITE(6,*) ' I1MPF, I2MPF ', I1MPF,I2MPF
      END IF
*
      NK = KMAX - KMIN + 1
*. Type of kstrings
      IF(IOBTP.EQ.1) THEN
        KEL1 = IEL1(ICLS) - 1
        KEL3 = IEL3(ICLS)
      ELSE IF(IOBTP.EQ.2) THEN
        KEL1 = IEL1(ICLS)
        KEL3 = IEL3(ICLS)
      ELSE
        KEL1 = IEL1(ICLS)
        KEL3 = IEL3(ICLS) - 1
      END IF
      IF(JOBTP.EQ.1) THEN
        KEL1 = KEL1 - 1
        KEL3 = KEL3
      ELSE IF(JOBTP .EQ. 2 ) THEN
        KEL1 = KEL1
        KEL3 = KEL3
      ELSE
        KEL1 = KEL1
        KEL3 = KEL3 - 1
      END IF
C?    write(6,*) ' icls ', icls
C?    write(6,*) ' iel1 iel3 ',iel1(icls),iel3(icls)
C?    write(6,*) ' kel1 kel3 ',kel1,kel3
      KTYPE = 0
      DO 10 KKTYPE = 1, N2OCTP
       IF(I2EL1(KKTYPE).EQ.KEL1.AND.
     &    I2EL3(KKTYPE).EQ.KEL3) KTYPE = KKTYPE
   10 CONTINUE
C?    write(6,*) ' ktype ', ktype
      IF(KTYPE.EQ.0) THEN
        NK = 0
        IEND = 1
        GOTO 101
      END IF
*. Symmetry of K strings
C          SYMCOM(ITASK,IOBJ,I1,I2,I12)
      CALL SYMCOM_MCLR(2,4,IOBSM,JKSM,ISM)
      IF(JKSM.EQ.0) THEN
        CALL SYMCOM_MCLR(2,4,JOBSM      ,IKSM,ISM)
        IF(IKSM.EQ.0) THEN
          NK = 0
          IEND = 1
          GOTO 101
        ELSE
          CALL SYMCOM_MCLR(2,4,IOBSM,KSM,IKSM)
        END IF
      ELSE
        CALL SYMCOM_MCLR(2,4,JOBSM,KSM,JKSM)
      END IF
C?    write(6,*) ' JOBSM,KSM,JKSM ',JOBSM,KSM,JKSM
      IF(KSM.EQ.0) THEN
        NK = 0
        IEND = 1
        GOTO 101
      END IF
      KOFF = I2SSO(KTYPE,KSM)
      KEND = MIN(KMAX,N2SSO(KTYPE,KSM))
      IF( KEND .EQ. N2SSO(KTYPE,KSM)) THEN
        IEND = 1
      ELSE
        IEND = 0
      END IF
      NK = KEND - KMIN + 1
      IOFF = ISSO(ICLS,ISM)
*. Loop over iorb,jorb
      IF(IJORD.EQ.0) THEN
        NIJ = NIOB* NJOB
      ELSE
        NIJ = NIOB*(NIOB+1)/2
      END IF
      If (NKDim.gt.0) Then
         DO IJ = 1, NIJ
            Call ICopy(NKDIM,[0],0,I1(1,IJ),1)
            call dcopy_(NKDIM,[0.0d0],0,XI1S(1,IJ),1)
         END DO
      End If
*
      DO J = 1, NJOB
        IF(IJORD.EQ.1) THEN
          IMIN = J
        ELSE
          IMIN = 1
        END IF
c       DO I = IMIN,NIOB
c       IJ = IJ + 1
C     DO IJ = 1, NIJ
C       CALL NXTIJ(I,J,NIOB,NJOB,IJORD,NONEW)
c       IORB = IOBOFF-1+I
        JORB = JOBOFF-1+J
        IF(IJORD.EQ.1) THEN
          IJOFF = (J-1)*NIOB-(J-1)*(J-2)/2+1-imin+1
        ELSE
          IJOFF = (J-1)*NIOB+1
        END IF
        DO 100 KSTR = KOFF+KMIN-1 , KOFF+KEND- 1
* N-2 => N-1
          JKSTR = 0
          IF(I2MPF.EQ.1) THEN
            IF(I2MAPO((KSTR-1)*L2MP+JORB).EQ.JORB) THEN
               JKSTR = I2MAPS((KSTR-1)*L2MP+JORB)
            END IF
          ELSE IF (I2MPF.EQ.0) THEN
           IFST = MAX(1,JORB+I2MPL(KSTR)-NORB)
           DO JJORB = IFST, MIN(JORb,I2MPL(KSTR))
             IF(I2MAPO(I2MPO(KSTR)-1+JJORB).EQ.JORB) THEN
               JKSTR = I2MAPS(I2MPO(KSTR)-1+JJORB)
             END IF
           END DO
         END IF
         IF(JKSTR.EQ.0) GOTO 100
            IF(JKSTR.GT.0 ) THEN
              SIGN = 1.0D0
            ELSE
              JKSTR = - JKSTR
              SIGN = - 1.0D0
            END IF
* N-1 => N
            do i = imin,niob
            ij = ijoff + i - 1
            iorb = ioboff-1+i
            ISTR = 0
            IF(I1MPF.EQ.1) THEN
              IF(I1MAPO((JKSTR-1)*L1MP+IORB).EQ.IORB) THEN
               ISTR = I1MAPS((JKSTR-1)*L1MP+IORB)
              END IF
            ELSE IF (I1MPF.EQ.0) THEN
             IFST = MAX(1,IORB+NORB-I1MPL(JKSTR))
             DO IIORB = IFST, MIN(IORB,I1MPL(JKSTR))
               IF(I1MAPO(I1MPO(JKSTR)-1+IIORB).EQ.IORB) THEN
                 ISTR = I1MAPS(I1MPO(JKSTR)-1+IIORB)
               END IF
             END DO
           END IF
           IF(ISTR.EQ.0) GOTO 99
*. Synthesis
           IF(ISTR.GT.0) THEN
              I1(KSTR-KOFF-KMIN+2,IJ) = ISTR - IOFF + 1
              XI1S(KSTR-KOFF-KMIN+2,IJ) = SIGN
           ELSE
             I1(KSTR-KOFF-KMIN+2,IJ) = -ISTR - IOFF + 1
             XI1S(KSTR-KOFF-KMIN+2,IJ) = -SIGN
           END IF
   99      CONTINUE
           END DO
*
  100   CONTINUE
*
C     END DO
      END DO
  101 CONTINUE
*
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Output from ADADS1 '
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' Number of K strings accessed ', NK
        IF(NK.NE.0) THEN
          DO IJ = 1, NIJ
            WRITE(6,*) ' Excited strings and sign for ij = ',IJ
            CALL IWRTMA(I1(1,IJ),1,NK,1,NK)
            CALL WRTMAT(XI1S(1,IJ),1,NK,1,NK)
          END DO
        END IF
      END IF
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer_array(NSSO)
        CALL Unused_integer(NORB1)
        CALL Unused_integer(NORB2)
        CALL Unused_integer(NORB3)
      END IF
      END
