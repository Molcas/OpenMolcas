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

subroutine UpG_T1(T1)
! upgrade T1

use chcc_global, only: no, nv, T1c
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: T1(nv,no)
integer(kind=iwp) :: a, i

do i=1,no
  do a=1,nv
    T1c(a,i) = T1(a,i)
  end do
end do

return

end subroutine UpG_T1
