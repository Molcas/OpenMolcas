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

subroutine ECP_shells(iAtmNr,List)

use Definitions, only: iwp, u6

implicit none
#include "itmax.fh"
integer(kind=iwp), intent(in) :: iAtmNr
integer(kind=iwp), intent(out) :: List(0:iTabMx)

!                                                                      *
!***********************************************************************
!                                                                      *
List(:) = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Dummy
if (iAtmNr == 0) then
  list(0) = 0
  list(1) = 0
  list(2) = 0
  list(3) = 0
! H-He
else if ((iAtmNr > 2) .and. (iAtmNr <= 4)) then
  list(0) = 1
  list(1) = 0
  list(2) = 0
  list(3) = 0
else if (iAtmNr <= 10) then
  list(0) = 1
  list(1) = 1
  list(2) = 0
  list(3) = 0
else if (iAtmNr <= 12) then
  list(0) = 1
  list(1) = 0
  list(2) = 0
  list(3) = 0
else if (iAtmNr <= 18) then
  list(0) = 1
  list(1) = 1
  list(2) = 0
  list(3) = 0
else if (iAtmNr <= 20) then
  list(0) = 1
  list(1) = 0
  list(2) = 0
  list(3) = 0
else if (iAtmNr <= 30) then
  list(0) = 1
  list(1) = 0
  list(2) = 1
  list(3) = 0
else if (iAtmNr <= 36) then
  list(0) = 1
  list(1) = 1
  list(2) = 1
  list(3) = 0
else if (iAtmNr <= 38) then
  list(0) = 1
  list(1) = 0
  list(2) = 0
  list(3) = 0
else if (iAtmNr <= 48) then
  list(0) = 1
  list(1) = 0
  list(2) = 1
  list(3) = 0
else if (iAtmNr <= 54) then
  list(0) = 1
  list(1) = 1
  list(2) = 1
  list(3) = 0
else if (iAtmNr <= 56) then
  list(0) = 1
  list(1) = 0
  list(2) = 0
  list(3) = 0
else if (iAtmNr <= 70) then
  list(0) = 1
  list(1) = 0
  list(2) = 0
  list(3) = 1
else if (iAtmNr <= 80) then
  list(0) = 1
  list(1) = 0
  list(2) = 1
  list(3) = 1
else if (iAtmNr <= 86) then
  list(0) = 1
  list(1) = 1
  list(2) = 1
  list(3) = 1
else if (iAtmNr <= 88) then
  list(0) = 1
  list(1) = 0
  list(2) = 0
  list(3) = 0
else if (iAtmNr <= 102) then
  list(0) = 1
  list(1) = 0
  list(2) = 0
  list(3) = 1
else if (iAtmNr <= 112) then
  list(0) = 1
  list(1) = 0
  list(2) = 1
  list(3) = 1
else if (iAtmNr <= 118) then
  list(0) = 1
  list(1) = 1
  list(2) = 1
  list(3) = 1
else if (iAtmNr <= 120) then
  list(0) = 1
  list(1) = 0
  list(2) = 0
  list(3) = 0
else
  write(u6,*) 'ECP_shells can not handle atom numbers beyond 112.'
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine ECP_Shells
