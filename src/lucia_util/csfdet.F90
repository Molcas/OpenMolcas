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

!#define _DEBUGPRINT_
subroutine CSFDET(NOPEN,IDET,NDET,ICSF,NCSF,CDC,SCR,nSCR,PSSIGN)
! Expand csf's in terms of combinations with
! the use of the Grabenstetter method (I.J.Q.C. 10, P142 (1976))
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

use lucia_data, only: IPRCIX
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nOpen, nDet, IDET(NOPEN,NDET), nCSF, ICSF(NOPEN,NCSF), nSCR
real(kind=wp), intent(out) :: CDC(NDET,NCSF), SCR(nSCR)
real(kind=wp), intent(in) :: PSSIGN
integer(kind=iwp) :: IOPEN, JCSF, JDADD, JDET, KLFREE, KLMDET, KLSCSF
real(kind=wp) :: CMBFAC, COEF, SGN

if (PSSIGN == Zero) then
  CMBFAC = One
else
  CMBFAC = sqrt(Two)
end if

KLFREE = 1
KLMDET = KLFREE
KLFREE = KLMDET+NDET*NOPEN
KLSCSF = KLFREE
KLFREE = KLSCSF+NOPEN

! OBTAIN INTERMEDIATE VALUES OF MS FOR ALL DETERMINANTS
do JDET=1,NDET
  call MSSTRN(IDET(:,JDET),SCR(KLMDET+(JDET-1)*NOPEN),NOPEN)
end do

do JCSF=1,NCSF
# ifdef _DEBUGPRINT
  write(u6,*) ' ....Output for CSF ',JCSF
# endif

  ! OBTAIN INTERMEDIATE COUPLINGS FOR CSF
  call MSSTRN(ICSF(:,JCSF),SCR(KLSCSF),NOPEN)

  do JDET=1,NDET
    ! EXPANSION COEFFICIENT OF DETERMINANT JDET FOR CSF JCSF
    COEF = One
    SGN = One
    JDADD = (JDET-1)*NOPEN
    do IOPEN=1,NOPEN

      if ((ICSF(IOPEN,JCSF) == 1) .and. (IDET(IOPEN,JDET) == 1)) then
        ! + + CASE
        COEF = COEF*(SCR(KLSCSF-1+IOPEN)+SCR(KLMDET-1+JDADD+IOPEN))/(Two*SCR(KLSCSF-1+IOPEN))
      else if ((ICSF(IOPEN,JCSF) == 1) .and. (IDET(IOPEN,JDET) == 0)) then
        ! + - CASE
        COEF = COEF*(SCR(KLSCSF-1+IOPEN)-SCR(KLMDET-1+JDADD+IOPEN))/(Two*SCR(KLSCSF-1+IOPEN))
      else if ((ICSF(IOPEN,JCSF) == 0) .and. (IDET(IOPEN,JDET) == 1)) then
        ! - + CASE
        COEF = COEF*(SCR(KLSCSF-1+IOPEN)-SCR(KLMDET-1+JDADD+IOPEN)+One)/(Two*SCR(KLSCSF-1+IOPEN)+Two)
        SGN = -SGN
      else if ((ICSF(IOPEN,JCSF) == 0) .and. (IDET(IOPEN,JDET) == 0)) then
        ! - - CASE
        COEF = COEF*(SCR(KLSCSF-1+IOPEN)+SCR(KLMDET-1+JDADD+IOPEN)+One)/(Two*SCR(KLSCSF-1+IOPEN)+Two)
      end if
    end do
    CDC(JDET,JCSF) = SGN*CMBFAC*sqrt(COEF)
  end do
end do

if (IPRCIX >= 5) then
  write(u6,*)
  write(u6,'(A,2I2)') '  The CDC array for  NOPEN ',NOPEN
  write(u6,*) ' NDET, NCSF = ',NDET,NCSF
  write(u6,*)
  call WRTMAT(CDC,NDET,NCSF,NDET,NCSF)
end if

end subroutine CSFDET
