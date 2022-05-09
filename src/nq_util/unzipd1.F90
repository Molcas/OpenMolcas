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
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nD1MO
real(kind=wp) :: D1Unzip(NASHT**2), D1MO(nD1MO)
! Input: nD1MO D1MO
! Output: D1Unzip
integer(kind=iwp) :: iLoc1, iLoc2, iLoc3, iv, ix

call FZero(D1Unzip,NASHT**2)
do iv=1,NASHT
  do ix=1,iv-1
    iLoc1 = (iv-1)*NASHT+ix
    iLoc2 = (ix-1)*NASHT+iv
    iLoc3 = (iv-1)*iv/2+ix
    D1Unzip(iLoc1) = Half*D1MO(iLoc3)
    D1Unzip(iLoc2) = D1Unzip(iLoc1)
  end do
  ix = iv
  iLoc1 = (iv-1)*NASHT+ix
  iLoc3 = (iv+1)*iv/2
  D1Unzip(iLoc1) = Half*D1MO(iLoc3)
end do

return

end subroutine UnzipD1
