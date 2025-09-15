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
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************

subroutine G2qtoG2r(G2r,G2q,nG2,nG2r)

use Index_Functions, only: iTri
use input_mclr, only: ntAsh
use Constants, only: One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nG2, nG2r
real(kind=wp), intent(out) :: G2r(nG2r)
real(kind=wp), intent(in) :: G2q(nG2)
integer(kind=iwp) :: iB, iDij, iDkl, iijkl, iRij, iRijkl, iRkl, jB, kB, lB
real(kind=wp) :: Fact

do iB=1,ntash
  do jB=1,ntash
    iDij = iTri(ib,jB)
    iRij = jb+(ib-1)*ntash
    do kB=1,ntash
      do lB=1,ntash
        iDkl = iTri(kB,lB)
        iRkl = lb+(kb-1)*ntash
        fact = One
        if ((iDij >= iDkl) .and. (kB == lB)) fact = Two
        if ((iDij < iDkl) .and. (iB == jB)) fact = Two
        iijkl = iTri(iDij,iDkl)
        iRijkl = iTri(iRij,iRkl)
        G2r(iRijkl) = Fact*G2q(iijkl)
      end do
    end do
  end do
end do

end subroutine G2qtoG2r
