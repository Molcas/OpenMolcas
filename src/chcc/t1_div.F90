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

subroutine T1_div(T1n,OE,no,nv)
! this routine does:
! T1n(a,i) = T1n(a,i)/(e(i)-e(a))
!
! division of T1n amplitides by denominator

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: no, nv
real(kind=wp) :: T1n(nv,no), OE(no+nv)
integer(kind=iwp) :: a, i
real(kind=wp) :: ei

do i=1,no
  ei = OE(i)
  do a=1,nv
    t1n(a,i) = t1n(a,i)/(ei-OE(no+a))
  end do
end do

return

end subroutine T1_div
