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

function arSin(Arg)

use Constants, only: One
use Definitions, only: wp

implicit none
real(kind=wp) :: arSin
real(kind=wp), intent(in) :: Arg
real(kind=wp) :: A
character(len=72) :: Warning
real(kind=wp), parameter :: Delta = 1.0e-12_wp

A = Arg
if (abs(A) > One) then
  write(Warning,3) A
  if (abs(A) >= One+Delta) then
    call WarningMessage(2,Warning)
    call Abend()
  end if
  !call WarningMessage(1,Warning)
  A = sign(One,A)
end if

ArSin = asin(A)

return

3 format(1X,'Warning argument of aSin= ',1F21.18)

end function arSin
