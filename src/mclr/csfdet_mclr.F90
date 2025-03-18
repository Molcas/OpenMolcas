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
      SUBROUTINE CSFDET_MCLR(NOPEN,IDET,NDET,ICSF,NCSF,CDC,PSSIGN,      &
     &                  IPRCSF)

! Expand csf's in terms of combinations with
! the use of the Graebenstetter method ( I.J.Q.C.10,P142(1976) )
!
! Input :
!         NOPEN : NUMBER OF OPEN ORBITALS
!         IDET  : OCCUPATION OF combinations
!         NDET  : NUMBER OF combinations
!         ICSF  : INTERMEDIATE SPIN COUPLINGS OF
!                 CSF'S IN BRANCHING DIAGRAM
! Output :
!         CDC :  NDET X NCSF MATRIX
!                GIVING EXPANSION FROM COMB'S TO CSF,S
!                CSF BASIS = Comb basis *CDC
!
! If combinations are use ( signaled by PSSIGN .ne. 0 )
! the factors are multiplies with sqrt(2), corresponding to
! a combination being 1/sqrt(2) times the sum or difference of two
! determinants
!
! The terms are not mutiplied with any sqrt(2), so the transformation is to
! the determinant normalization
!
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT NONE
      INTEGER NOPEN,NDET,NCSF
      INTEGER IDET(NOPEN,NDET),ICSF(NOPEN,NCSF)
      REAL*8 CDC(NDET,NCSF)
      Real*8 PSSIGN
      Integer IPRCSF

!     Local variables
      Real*8, Allocatable:: LMDET(:), lSCSF(:)
      INTEGER NTEST,JDET,JDADD,IOPEN,JCSF
      REAL*8 CMBFAC,COEF,SIGN

      NTEST = 0000
      NTEST = MAX(IPRCSF,NTEST)
      IF(PSSIGN.EQ.0.0D0) THEN
       CMBFAC = 1.0D0
      ELSE
       CMBFAC = SQRT(2.0D0)
      END IF
      CALL mma_allocate(LMDET,NDET*NOPEN,Label='LMDET')
      CALL mma_allocate(LSCSF,NDET*NOPEN,Label='LSCSF')
!
!.. OBTAIN INTERMEDIATE VALUES OF MS FOR ALL DETERMINANTS
      DO 10 JDET = 1, NDET
        CALL MSSTRN_MCLR(IDET(1,JDET),LMDET(1+(JDET-1)*NOPEN),NOPEN)
   10 CONTINUE

!
      DO 1000 JCSF = 1, NCSF
       IF( NTEST .GE. 105 ) WRITE(6,*) ' ....Output for CSF ',JCSF
!
! OBTAIN INTERMEDIATE COUPLINGS FOR CSF
      CALL MSSTRN_MCLR(ICSF(1,JCSF),LSCSF,NOPEN)
!
      DO 900 JDET = 1, NDET
! EXPANSION COEFFICIENT OF DETERMINANT JDET FOR CSF JCSF
      COEF = 1.0D0
      SIGN = 1.0D0
      JDADD = (JDET-1)*NOPEN
      DO 700 IOPEN = 1, NOPEN
!
! + + CASE
        IF(ICSF(IOPEN,JCSF).EQ.1.AND.IDET(IOPEN,JDET).EQ.1) THEN
          COEF = COEF * (LSCSF(IOPEN)+LMDET(JDADD+IOPEN) )              &
     &         / (2.0D0*LSCSF(IOPEN) )
        ELSE IF(ICSF(IOPEN,JCSF).EQ.1.AND.IDET(IOPEN,JDET).EQ.0) THEN
! + - CASE
          COEF = COEF *(LSCSF(IOPEN)-LMDET(JDADD+IOPEN) )               &
     &         / (2.0D0*LSCSF(IOPEN) )
        ELSE IF(ICSF(IOPEN,JCSF).EQ.0.AND.IDET(IOPEN,JDET).EQ.1) THEN
! - + CASE
          COEF = COEF * (LSCSF(IOPEN)-LMDET(JDADD+IOPEN) +1.0D0)        &
     &         / (2.0D0*LSCSF(IOPEN)+2.0D0 )
          SIGN  = - SIGN
        ELSE IF(ICSF(IOPEN,JCSF).EQ.0.AND.IDET(IOPEN,JDET).EQ.0) THEN
! - - CASE
          COEF = COEF * (LSCSF(IOPEN)+LMDET(JDADD+IOPEN) +1.0D0)        &
     &         / (2.0D0*LSCSF(IOPEN)+2.0D0 )
        END IF
  700 CONTINUE
       CDC(JDET,JCSF) = SIGN * CMBFAC * SQRT(COEF)
  900 CONTINUE
 1000 CONTINUE
!
      CALL mma_deallocate(LSCSF)
      CALL mma_deallocate(LMDET)
!
      IF( NTEST .GE. 5 ) THEN
        WRITE(6,*)
        WRITE(6,'(A,2I2)')                                              &
     &  '  The CDC array for  NOPEN ',NOPEN
        WRITE(6,*)
        CALL WRTMAT(CDC,NDET,NCSF,NDET,NCSF)
      END IF

!
      END SUBROUTINE CSFDET_MCLR
