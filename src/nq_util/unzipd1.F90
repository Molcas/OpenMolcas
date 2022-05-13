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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************

! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Dec. 08, 2021, created this file.               *
! ****************************************************************
subroutine UnzipD1(D1Unzip,D1MO,nD1MO)

use nq_Info, only: NASHT
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: D1Unzip(NASHT,NASHT)
integer(kind=iwp), intent(in) :: nD1MO
real(kind=wp), intent(in) :: D1MO(nD1MO)
integer(kind=iwp) :: iLoc, iv, ix

D1Unzip(:,:) = Zero
do iv=1,NASHT
  do ix=1,iv-1
    iLoc = (iv-1)*iv/2+ix
    D1Unzip(ix,iv) = Half*D1MO(iLoc)
    D1Unzip(iv,ix) = D1Unzip(ix,iv)
  end do
  ix = iv
  iLoc = (iv+1)*iv/2
  D1Unzip(ix,iv) = Half*D1MO(iLoc)
end do

return

end subroutine UnzipD1
