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
      SUBROUTINE PRCMAT(NSS,XMATR,XMATI)
      use Definitions, only: u6
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NSS
      REAL*8 XMATR(NSS,NSS),XMATI(NSS,NSS)

      INTEGER JSTA, JEND, ISS, JSS
C Write out matrix elements over states as a complex matrix
C in square format
      DO JSTA=1,NSS,2
       JEND=MIN(NSS,JSTA+1)
       WRITE(u6,*)
       WRITE(u6,'(1X,A8,12X,I3,35X,I3)')' STATE  ',(JSS,JSS=JSTA,JEND)
       DO ISS=1,NSS
       WRITE(u6,'(1X,I4,2x,2(A1,F10.6,A1,F10.6,A1,3x))')
     &           ISS,('(',XMATR(ISS,JSS),',',XMATI(ISS,JSS),
     &           ')',JSS=JSTA,JEND)
       END DO
      END DO
      END SUBROUTINE PRCMAT

      SUBROUTINE PRCMAT2(INPUT,NSS,XMATR,XMATI)
      use Cntrl, only: ISOCMP, SOPRNM
      IMPLICIT NONE
      INTEGER INPUT, NSS
      REAL*8 XMATR(NSS,NSS),XMATI(NSS,NSS)

      CHARACTER(LEN=8) PROPERTY
      CHARACTER(LEN=1) DIRECTION
      CHARACTER(LEN=200) FILENAME
      INTEGER LU, JSTA, ISS
      Integer, External:: IsFreeUnit
C Write out matrix elements over states as a complex matrix
C in parsable format
      if(INPUT.gt.0) THEN
        PROPERTY = SOPRNM(INPUT)
      ELSE
        PROPERTY = 'EIGVEC'
      ENDIF

      WRITE(DIRECTION,'(I1)') ISOCMP(INPUT)
      IF (PROPERTY(1:5).EQ."MLTPL") THEN
          IF (PROPERTY(8:8).EQ.'0') THEN
            FILENAME = 'monopole-'//DIRECTION//'.txt'
          ELSE IF (PROPERTY(8:8).EQ.'1') THEN
            FILENAME = 'dipole-'//DIRECTION//'.txt'
          ELSE IF (PROPERTY(8:8).EQ.'2') THEN
            FILENAME = 'quadrupole-'//DIRECTION//'.txt'
          ELSE
            RETURN
          END IF
      ELSE IF (PROPERTY(1:5).EQ."MLTPV") THEN
          IF (PROPERTY(8:8).EQ.'2') THEN
            FILENAME = 'velocity_quadrupole-'//DIRECTION//'.txt'
          ELSE
            RETURN
          END IF
      ELSE IF (PROPERTY(1:4).EQ."VELO") THEN
          FILENAME = 'velocity_dipole-'//DIRECTION//'.txt'
      ELSE IF (PROPERTY(1:4).EQ."ANGM") THEN
          FILENAME = 'angmom-'//DIRECTION//'.txt'
      ELSE IF (PROPERTY(1:6).EQ."EIGVEC") THEN
          FILENAME = "eigvectors.txt"
      ELSE
          RETURN
      END IF
      Lu = 88
      Lu = IsFreeUnit(Lu)
      OPEN(UNIT=LU,FILE=FILENAME,STATUS='REPLACE')
      WRITE(LU,*) "#NROW NCOL REAL IMAG"
      DO JSTA=1,NSS
        DO ISS=1,NSS
        WRITE(LU,'(I4,I4,A1,ES25.16,A1,ES25.16)') ISS,JSTA,' ',
     &   XMATR(ISS,JSTA),' ',XMATI(ISS,JSTA)
        END DO
      END DO
      CLOSE(LU)
      END SUBROUTINE PRCMAT2

      SUBROUTINE PRCMAT3(NSS,SMATR,SMATI,DIR)
      IMPLICIT None
      INTEGER NSS
      REAL*8 SMATR(NSS,NSS), SMATI(NSS,NSS)
      INTEGER DIR
      CHARACTER(LEN=1) DIRECTION
      CHARACTER(LEN=200) FILENAME
      Integer LU, JSTA, ISS
      Integer, External:: IsFreeUnit
C Write out spin matrix elements in parsable format
      WRITE(DIRECTION,'(I1)') DIR
      FILENAME = 'spin-'//DIRECTION//'.txt'
      Lu = 88
      Lu = IsFreeUnit(Lu)
      OPEN(UNIT=Lu,FILE=FILENAME,STATUS='REPLACE')
      WRITE(Lu,*) "#NROW NCOL REAL IMAG"
      DO JSTA=1,NSS
        DO ISS=1,NSS
        WRITE(Lu,'(I4,I4,A1,ES25.16,A1,ES25.16)') ISS,JSTA,' ',
     &   SMATR(ISS,JSTA),' ',SMATI(ISS,JSTA)
        END DO
      END DO
      CLOSE(Lu)
      END SUBROUTINE PRCMAT3

      SUBROUTINE MULMAT(NSS,XMATR,XMATI,ee,Z)
      IMPLICIT None
      INTEGER NSS
      REAL*8 XMATR(NSS,NSS),XMATI(NSS,NSS), ee
      COMPLEX*16 Z(NSS,NSS)

      Integer ISS,JSS
      ee=0.d0
      Z(:,:)=(0.0d0,0.0d0)

      DO ISS=1,NSS
      DO JSS=1,NSS
      ee=ee+XMATR(ISS,JSS)*XMATR(ISS,JSS)+
     & XMATI(ISS,JSS)*XMATI(ISS,JSS)
      Z(ISS,JSS)=Z(ISS,JSS)+
     &CMPLX(XMATR(ISS,JSS),XMATI(ISS,JSS),kind=8)
      enddo
      enddo
      END SUBROUTINE MULMAT

