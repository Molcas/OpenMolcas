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
      SUBROUTINE ADS1(NK,I1,XI1S,LI1,IORB,LORB,
     &                ICLS,ISM,IMAPO,IMAPS,IMPL,IMPO,IMPF,LMAP,
     &                IEL1,IEL3,
     &                I1EL1,I1EL3,ISSO,NSSO,I1SSO,N1SSO,NOCTP,N1OCTP,
     &                NORB1,NORB2,NORB3,ORBSM,NORB,KMAX,KMIN,IEND)
*
* Obtain I1(KSTR) = +/- A+ IORB !KSTR>
*
* KSTR is restricted to strings with relative numbers in the
* range KMAX to KMIN
* =====
* Input
* =====
* IORB : Firat orbital to be added
* LORB : Number of orbitals to be added : IORB to IORB-1+LORB
*        are used. They must all be in the same TS group
* ICLS,ISM : Class and symmetry of string with added electron
* IMAPO,IMAPS : map from Kstrings to Istrings
* IEL1(3) : Number of electrons in RAS1(3) for I strings
* I1EL1(3) : Number of electrons in RAS1(3) for K strings
* ISSO : TS symmetry offset for I strings
* NSSO : Number of TS strings for I strings
* I1SSO : TS symmetry offset for K strings
* N1SSO : Number of TS strings for K strings
* NOCTP : Number of occupation types for I strings
* N1OCTP : Number of occupation types for K strings
* NORB1(2,3) : Number of RAS1(2,3) orbitals
* IORBSM : Orbital symmety array
* NORB : Number of active  orbitals
* KMAX : Largest allowed relative number for K strings
*        If Kmax is set to -1 all strings are searched
* KMIN : Smallest allowed relative number for K strings
*
* ======
* Output
* ======
*
* NK      : Number of K strings
* I1(KSTR,JORB) : ne. 0 => a + JORB !KSTR> = +/-!ISTR>
* XI1S(KSTR,JORB) : above +/-
*          : eq. 0    a + JORB !KSTR> = 0
* Offset is KMIN
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER IEL1(*),IEL3(*),I1EL1(*),I1EL3(*)
      INTEGER ISSO(NOCTP,*),NSSO(NOCTP,*)
      INTEGER I1SSO(N1OCTP,*),N1SSO(N1OCTP,*)
      INTEGER ORBSM(*)
C     INTEGER IMAPO(NORB,*),IMAPS(NORB,*)
      INTEGER IMAPO(*),IMAPS(*)
      INTEGER IMPL(*),IMPO(*)
*.Output
      INTEGER I1(*)
      DIMENSION XI1S(*)
*
      LDIM = 0 ! dummy initilize
      NTEST = 000
      IF(NTEST.NE.0) THEN
       WRITE(6,*) ' ============== '
       WRITE(6,*) ' ADSTS speaking '
       WRITE(6,*) ' ============== '
       WRITE(6,*) ' IORB,ISM,ICLS',IORB,ISM,ICLS
       WRITE(6,*) ' IMPF, LMAP ', IMPF, LMAP
       WRITE(6,*) ' N1SSO : '
       CALL IWRTMA(N1SSO,N1OCTP,8,N1OCTP,8)
      END IF
      NK = KMAX - KMIN + 1
*. Type of kstrings
      IF(IORB.LE.NORB1) THEN
        KEL1 = IEL1(ICLS) - 1
        KEL3 = IEL3(ICLS)
      ELSE IF(IORB.LE.NORB1+NORB2) THEN
        KEL1 = IEL1(ICLS)
        KEL3 = IEL3(ICLS)
      ELSE
        KEL1 = IEL1(ICLS)
        KEL3 = IEL3(ICLS) - 1
      END IF
      KTYPE = 0
C      write(6,*) ' N1OCTP ', N1OCTP
      DO 10 KKTYPE = 1, N1OCTP
       IF(I1EL1(KKTYPE).EQ.KEL1.AND.
     &    I1EL3(KKTYPE).EQ.KEL3) KTYPE = KKTYPE
   10 CONTINUE
C      write(6,*) ' kel1 kel3 ktype ',KEL1,KEL3,KTYPE
      IF(KTYPE.EQ.0) THEN
        NK = 0
        IEND = 1
        GOTO 101
      END IF
*. Symmetry of K strings
C          SYMCOM(ITASK,IOBJ,I1,I2,I12)
      CALL SYMCOM_MCLR(2,4,ORBSM(IORB),KSM,ISM)
      IF(KSM.EQ.0) THEN
        NK = 0
        IEND = 1
        GOTO 101
      END IF
      KOFF = I1SSO(KTYPE,KSM)
C?    WRITE(6,*) ' KTYPE KSM ', KTYPE,KSM
      IF(KMAX.EQ.-1) THEN
        KEND = N1SSO(KTYPE,KSM)
      ELSE
        KEND = MIN(N1SSO(KTYPE,KSM),KMAX)
      END IF
      IF(KEND.LT.N1SSO(KTYPE,KSM)) THEN
        IEND = 0
      ELSE
        IEND = 1
      END IF
      NK = KEND-KMIN+1
      IF(KMAX.EQ.-1) THEN
       LDIM = NK
      ELSE
       LDIM = LI1
      END IF
C?    IF(KMAX.EQ.-1) WRITE(6,*) ' KMAX = -1, LDIM=',LDIM
      IOFF = ISSO(ICLS,ISM)
      KSUB = KOFF+KMIN-2
      DO 110 IIORB = IORB,IORB+LORB-1
      IORBR = IIORB-IORB+1
      DO 100 KSTR = KOFF+KMIN-1 , KOFF+KEND-1
C        write(6,*) ' KSTR = ',KSTR
        KREL = KSTR - KSUB
*
        ISTR = 0
        IF(IMPF.EQ.1) THEN
          IF(IMAPO((KSTR-1)*LMAP+IIORB).EQ.IIORB)
     &    ISTR = IMAPS((KSTR-1)*LMAP+IIORB)
        ELSE
C          write(6,*) ' IMPL = ',IMPL(KSTR)
C          write(6,*) ' IMPO = ',IMPO(KSTR)
          DO IIIORB = 1, IMPL(KSTR)
           IF(IMAPO(IMPO(KSTR)-1+IIIORB).EQ.IIORB)
     &     ISTR = IMAPS(IMPO(KSTR)-1+IIIORB)
          END DO
        END IF
        IF(ISTR.GT.0 ) THEN
          I1(KREL+(IORBR-1)*LDIM) = ISTR - IOFF + 1
          XI1S(KREL+(IORBR-1)*LDIM) = 1.0D0
        ELSE IF (ISTR.LT.0 ) THEN
          I1(KREL+(IORBR-1)*LDIM) = -ISTR - IOFF + 1
          XI1S(KREL+(IORBR-1)*LDIM) = -1.0D0
        ELSE  IF( ISTR.EQ.0) THEN
          I1(KREL+(IORBR-1)*LDIM) = 0
          XI1S(KREL+(IORBR-1)*LDIM) = 0.0D0
        END IF
  100 CONTINUE
  110 CONTINUE
  101 CONTINUE
*
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Output from ASTR '
        WRITE(6,*) ' ================ '
        WRITE(6,*) ' Number of K strings accessed ', NK
        IF(NK.NE.0) THEN
          DO 200 IIORB = IORB,IORB+LORB-1
            IIORBR = IIORB-IORB+1
            WRITE(6,*) ' Info for orbital ', IIORB
            WRITE(6,*) ' Excited strings and sign '
            CALL IWRTMA(I1  (1+(IIORBR-1)*LDIM),1,NK,1,NK)
            CALL WRTMAT(XI1S(1+(IIORBR-1)*LDIM),1,NK,1,NK)
  200     CONTINUE
        END IF
      END IF
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer_array(NSSO)
        CALL Unused_integer(NORB3)
        CALL Unused_integer(NORB)
      END IF
      END
