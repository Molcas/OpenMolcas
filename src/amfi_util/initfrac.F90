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

subroutine initfrac(nprimit1,nprimit2,nprimit3,nprimit4,quot1,quot2,expo1,expo2,expo3,expo4)
!bs initialize some arrays with factors  needed for cfunct(x)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nprimit1, nprimit2, nprimit3, nprimit4
real(kind=wp), intent(out) :: quot1(nprimit1,nprimit2,nprimit3,nprimit4), quot2(nprimit1,nprimit2,nprimit3,nprimit4)
real(kind=wp), intent(in) :: expo1(*), expo2(*), expo3(*), expo4(*)
integer(kind=iwp) :: irun2, irun3, irun4
real(kind=wp) :: sum24

do irun4=1,nprimit4
  do irun3=1,nprimit3
    do irun2=1,nprimit2
      sum24 = expo2(irun2)+expo4(irun4)
      quot1(1:nprimit1,irun2,irun3,irun4) = One/(One+(expo1(1:nprimit1)+expo3(irun3))/sum24)
    end do
  end do
end do
do irun4=1,nprimit4
  do irun3=1,nprimit3
    do irun2=1,nprimit2
      sum24 = expo2(irun2)+expo4(irun4)
      quot2(1:nprimit1,irun2,irun3,irun4) = One/(One+sum24/(expo1(1:nprimit1)+expo3(irun3)))
    end do
  end do
end do

return

end subroutine initfrac
