!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine CSFDET_LUCIA(NOPEN,IDET,NDET,ICSF,NCSF,CDC,SCR,nSCR,PSSIGN,IPRCSF)
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
! Scratch :
!          SCR, SHOULD AT LEAST BE ???
!
! If combinations are use (signaled by PSSIGN /= 0)
! the factors are multiplies with sqrt(2), corresponding to
! a combination being 1/sqrt(2) times the sum or difference of two
! determinants

implicit real*8(A-H,O-Z)
integer nOpen, nDet, nCSF, nSCR
integer IDET(NOPEN,NDET), ICSF(NOPEN,NCSF)
real*8 CDC(NDET,NCSF)
real*8 SCR(nSCR)

NTEST = 0
NTEST = max(IPRCSF,NTEST)
if (PSSIGN == 0.0d0) then
  CMBFAC = 1.0d0
else
  CMBFAC = sqrt(2.0d0)
end if

KLFREE = 1
KLMDET = KLFREE
KLFREE = KLMDET+NDET*NOPEN
KLSCSF = KLFREE
KLFREE = KLSCSF+NOPEN

! OBTAIN INTERMEDIATE VALUES OF MS FOR ALL DETERMINANTS
do JDET=1,NDET
  call MSSTRN_LUCIA(IDET(1,JDET),SCR(KLMDET+(JDET-1)*NOPEN),NOPEN,IPRCSF)
end do

do JCSF=1,NCSF
  if (NTEST >= 105) write(6,*) ' ....Output for CSF ',JCSF

  ! OBTAIN INTERMEDIATE COUPLINGS FOR CSF
  call MSSTRN_LUCIA(ICSF(1,JCSF),SCR(KLSCSF),NOPEN,IPRCSF)

  do JDET=1,NDET
    ! EXPANSION COEFFICIENT OF DETERMINANT JDET FOR CSF JCSF
    COEF = 1.0d0
    SIGN = 1.0d0
    JDADD = (JDET-1)*NOPEN
    do IOPEN=1,NOPEN

      if ((ICSF(IOPEN,JCSF) == 1) .and. (IDET(IOPEN,JDET) == 1)) then
        ! + + CASE
        COEF = COEF*(SCR(KLSCSF-1+IOPEN)+SCR(KLMDET-1+JDADD+IOPEN))/(2.0d0*SCR(KLSCSF-1+IOPEN))
      else if ((ICSF(IOPEN,JCSF) == 1) .and. (IDET(IOPEN,JDET) == 0)) then
        ! + - CASE
        COEF = COEF*(SCR(KLSCSF-1+IOPEN)-SCR(KLMDET-1+JDADD+IOPEN))/(2.0d0*SCR(KLSCSF-1+IOPEN))
      else if ((ICSF(IOPEN,JCSF) == 0) .and. (IDET(IOPEN,JDET) == 1)) then
        ! - + CASE
        COEF = COEF*(SCR(KLSCSF-1+IOPEN)-SCR(KLMDET-1+JDADD+IOPEN)+1.0d0)/(2.0d0*SCR(KLSCSF-1+IOPEN)+2.0d0)
        SIGN = -SIGN
      else if ((ICSF(IOPEN,JCSF) == 0) .and. (IDET(IOPEN,JDET) == 0)) then
        ! - - CASE
        COEF = COEF*(SCR(KLSCSF-1+IOPEN)+SCR(KLMDET-1+JDADD+IOPEN)+1.0d0)/(2.0d0*SCR(KLSCSF-1+IOPEN)+2.0d0)
      end if
    end do
    CDC(JDET,JCSF) = SIGN*CMBFAC*sqrt(COEF)
  end do
end do

if (NTEST >= 5) then
  write(6,*)
  write(6,'(A,2I2)') '  The CDC array for  NOPEN ',NOPEN
  write(6,*) ' NDET, NCSF = ',NDET,NCSF
  write(6,*)
  call WRTMAT(CDC,NDET,NCSF,NDET,NCSF)
end if

end subroutine CSFDET_LUCIA
