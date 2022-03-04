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

subroutine getpow(pmax,quot,quotpow,nprim1,nprim2,nprim3,nprim4)
!bs generates some powers of for the prefactors of cfunct(X)
!bs look out for details there and in initfrac

use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: pmax, nprim1, nprim2, nprim3, nprim4
real(kind=wp) :: quot(nprim1,nprim2,nprim3,nprim4), quotpow(nprim1,nprim2,nprim3,nprim4)

quotpow(:,:,:,:) = quot(:,:,:,:)**(real(pmax,kind=wp)-Half)

return

end subroutine getpow
