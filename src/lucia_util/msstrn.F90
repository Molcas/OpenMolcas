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
subroutine MSSTRN(INSTRN,UTSTRN,NOPEN)
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

use Constants, only: Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NOPEN, INSTRN(NOPEN)
real(kind=wp), intent(out) :: UTSTRN(NOPEN)
integer(kind=iwp) :: IOPEN

UTSTRN(1) = real(INSTRN(1),kind=wp)-Half
do IOPEN=2,NOPEN
  UTSTRN(IOPEN) = UTSTRN(IOPEN-1)+real(INSTRN(IOPEN),kind=wp)-Half
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' ... Output from MSSTRN'
write(u6,*) ' INSTRN AND UTSTRN'
call IWRTMA(INSTRN,1,NOPEN,1,NOPEN)
call WRTMAT(UTSTRN,1,NOPEN,1,NOPEN)
#endif

end subroutine MSSTRN
