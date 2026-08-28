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

!#define _DEBUGPRINT_
subroutine CSFDET_MCLR(NOPEN,IDET,NDET,ICSF,NCSF,CDC,PSSIGN)
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
!
! If combinations are use (signaled by PSSIGN /= 0)
! the factors are multiplies with sqrt(2), corresponding to
! a combination being 1/sqrt(2) times the sum or difference of two
! determinants
!
! The terms are not mutiplied with any sqrt(2), so the transformation is to
! the determinant normalization

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NOPEN, NDET, IDET(NOPEN,NDET), NCSF, ICSF(NOPEN,NCSF)
real(kind=wp), intent(out) :: CDC(NDET,NCSF)
real(kind=wp), intent(in) :: PSSIGN
integer(kind=iwp) :: IOPEN, JCSF, JDET
real(kind=wp) :: CMBFAC, COEF, SGN
real(kind=wp), allocatable :: LMDET(:,:), lSCSF(:)

if (PSSIGN == Zero) then
  CMBFAC = One
else
  CMBFAC = sqrt(Two)
end if
call mma_allocate(LMDET,NOPEN,NDET,Label='LMDET')
call mma_allocate(LSCSF,NOPEN,Label='LSCSF')

! OBTAIN INTERMEDIATE VALUES OF MS FOR ALL DETERMINANTS
do JDET=1,NDET
  call MSSTRN_MCLR(IDET(:,JDET),LMDET(:,JDET),NOPEN)
end do

do JCSF=1,NCSF
# ifdef _DEBUGPRINT_
  write(u6,*) ' ....Output for CSF ',JCSF
# endif

  ! OBTAIN INTERMEDIATE COUPLINGS FOR CSF
  call MSSTRN_MCLR(ICSF(1,JCSF),LSCSF,NOPEN)

  do JDET=1,NDET
    ! EXPANSION COEFFICIENT OF DETERMINANT JDET FOR CSF JCSF
    COEF = One
    SGN = One
    do IOPEN=1,NOPEN

      if (ICSF(IOPEN,JCSF) == 1) then
        if (IDET(IOPEN,JDET) == 1) then
          ! + + CASE
          COEF = COEF*(LSCSF(IOPEN)+LMDET(IOPEN,JDET))/(Two*LSCSF(IOPEN))
        else if (IDET(IOPEN,JDET) == 0) then
          ! + - CASE
          COEF = COEF*(LSCSF(IOPEN)-LMDET(IOPEN,JDET))/(Two*LSCSF(IOPEN))
        end if
      else if (ICSF(IOPEN,JCSF) == 0) then
        if (IDET(IOPEN,JDET) == 1) then
          ! - + CASE
          COEF = COEF*(LSCSF(IOPEN)-LMDET(IOPEN,JDET)+One)/(Two*LSCSF(IOPEN)+Two)
          SGN = -SGN
        else if (IDET(IOPEN,JDET) == 0) then
          ! - - CASE
          COEF = COEF*(LSCSF(IOPEN)+LMDET(IOPEN,JDET)+One)/(Two*LSCSF(IOPEN)+Two)
        end if
      end if
    end do
    CDC(JDET,JCSF) = SGN*CMBFAC*sqrt(COEF)
  end do
end do

call mma_deallocate(LSCSF)
call mma_deallocate(LMDET)

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,'(A,2I2)') '  The CDC array for  NOPEN ',NOPEN
write(u6,*)
call WRTMAT(CDC,NDET,NCSF,NDET,NCSF)
#endif

end subroutine CSFDET_MCLR
