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

subroutine put_temp_data_on_intgrl(LUINTMZ_,NSYMZ_,NORBZ_)

use Intgrl, only: IAD2M, LUINTMZ, NORBZ, NSYMZ
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: LUINTMZ_, NSYMZ_, NORBZ_(8)
integer(kind=iwp) :: iAddress

iAddress = 0
IAD2M(:,:) = 0
! read the address list from the existing file
call iDaFile(LUINTMZ_,2,IAD2M,size(IAD2M),iAddress)
NSYMZ = NSYMZ_
LUINTMZ = LUINTMZ_
NORBZ(1:NSYMZ_) = NORBZ_(1:NSYMZ_)

end subroutine put_temp_data_on_intgrl
