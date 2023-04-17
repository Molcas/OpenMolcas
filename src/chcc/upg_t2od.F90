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

subroutine UpG_T2od(T2,dima,adda,dimb,addb)
! upgrade T2 (off-diagonal)

use chcc_global, only: no, T2c
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dima, adda, dimb, addb
real(kind=wp) :: T2(dima,dimb,no,no)
integer(kind=iwp) :: a, b, i, j

do j=1,no
  do i=1,no
    do b=1,dimb
      do a=1,dima
        T2c(a+adda,b+addb,i,j) = T2(a,b,i,j)
        T2c(b+addb,a+adda,j,i) = T2(a,b,i,j)
      end do
    end do
  end do
end do

return

end subroutine UpG_T2od
