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
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
!*****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 07, 2022, created this file.               *
!*****************************************************************

subroutine AntiOneDFoil(TwoD,OneD,m,n)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: m, n
real(kind=wp) :: TwoD(m,n), OneD(m*n)
integer(kind=iwp) :: I, iLoc, J

iLoc = 1
do J=1,N
  do I=1,M
    TwoD(I,J) = OneD(iLoc)
    iLoc = iLoc+1
  end do
end do
return

end subroutine AntiOneDFoil
