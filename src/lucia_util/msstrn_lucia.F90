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

subroutine MSSTRN_LUCIA(INSTRN,UTSTRN,NOPEN,IPRCSF)
! A STRING IS GIVEN IN FORM A SEQUENCE OF ZEROES AND ONES
!
! REINTERPRET THIS AS :
!
! 1 : THE INPUT STRING IS A DETERMINANT AND THE
!     1'S INDICATE ALPHA ELECTRONS AND THE
!     0'S INDICATE BETA ELECTRONS.
!     UTSTRN IS THE MS-VALUES ATE EACH VERTEX
!
! 2 : THE INPUT STRING IS A CSF GIVEN IN A
!     BRANCHING DIAGRAM, WHERE
!     1'S INDICATE UPWARDS SPIN COUPLING
!     WHILE THE 0'S INDICATES DOWNWARDS SPIN COUPLING,
!     REEXPRESS THIS AS S VALUES OF ALL COUPLINGS
!
! THE TWO PROCEDURES ARE IDENTICAL.

implicit real*8(A-H,O-Z)
dimension INSTRN(NOPEN), UTSTRN(NOPEN)

UTSTRN(1) = dble(INSTRN(1))-0.5d0
do IOPEN=2,NOPEN
  UTSTRN(IOPEN) = UTSTRN(IOPEN-1)+dble(INSTRN(IOPEN))-0.5d0
end do

NTEST = 0
NTEST = max(NTEST,IPRCSF)
if (NTEST >= 10) then
  write(6,*) ' ... Output from MSSTRN'
  write(6,*) ' INSTRN AND UTSTRN'
  call IWRTMA(INSTRN,1,NOPEN,1,NOPEN)
  call WRTMAT(UTSTRN,1,NOPEN,1,NOPEN)
end if

end subroutine MSSTRN_LUCIA
