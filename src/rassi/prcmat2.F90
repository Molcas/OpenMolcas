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

subroutine PRCMAT2(INPUT,NSS,XMATR,XMATI)
! Write out matrix elements over states as a complex matrix
! in parsable format

use Definitions, only: iwp, wp
use Cntrl, only: ISOCMP, SOPRNM

implicit none
integer(kind=iwp), intent(in) :: INPUT, NSS
real(kind=wp), intent(in) :: XMATR(NSS,NSS), XMATI(NSS,NSS)
character(len=8) PROPERTY
character(len=1) DIRECTION
character(len=200) FILENAME
integer(kind=iwp) LU, JSTA, ISS
integer(kind=iwp), external :: IsFreeUnit

if (INPUT > 0) then
  PROPERTY = SOPRNM(INPUT)
else
  PROPERTY = 'EIGVEC'
end if

write(DIRECTION,'(I1)') ISOCMP(INPUT)
if (PROPERTY(1:5) == 'MLTPL') then
  if (PROPERTY(8:8) == '0') then
    FILENAME = 'monopole-'//DIRECTION//'.txt'
  else if (PROPERTY(8:8) == '1') then
    FILENAME = 'dipole-'//DIRECTION//'.txt'
  else if (PROPERTY(8:8) == '2') then
    FILENAME = 'quadrupole-'//DIRECTION//'.txt'
  else
    return
  end if
else if (PROPERTY(1:5) == 'MLTPV') then
  if (PROPERTY(8:8) == '2') then
    FILENAME = 'velocity_quadrupole-'//DIRECTION//'.txt'
  else
    return
  end if
else if (PROPERTY(1:4) == 'VELO') then
  FILENAME = 'velocity_dipole-'//DIRECTION//'.txt'
else if (PROPERTY(1:4) == 'ANGM') then
  FILENAME = 'angmom-'//DIRECTION//'.txt'
else if (PROPERTY(1:6) == 'EIGVEC') then
  FILENAME = 'eigvectors.txt'
else
  return
end if
Lu = IsFreeUnit(88)
open(unit=LU,file=FILENAME,status='REPLACE')
write(LU,*) '#NROW NCOL REAL IMAG'
do JSTA=1,NSS
  do ISS=1,NSS
    write(LU,'(I6,1X,I6,A1,ES25.16,A1,ES25.16)') ISS,JSTA,' ',XMATR(ISS,JSTA),' ',XMATI(ISS,JSTA)
  end do
end do
close(LU)

end subroutine PRCMAT2
