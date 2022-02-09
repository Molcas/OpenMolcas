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

function iRnge(Val,Bin,nBin)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iRnge
integer(kind=iwp), intent(in) :: nBin
real(kind=wp), intent(in) :: Val, Bin(nBin)
integer(kind=iwp) :: iBin

iRnge = nBin
do iBin=1,nBin-1
  if (Val > Bin(iBin)) then
    iRnge = iBin
    exit
  end if
end do

end function iRnge
