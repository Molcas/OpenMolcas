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

implicit none
integer no, nv
real*8 T1n(1:nv,1:no)
real*8 OE(1:no+nv)
! help variables
integer i, a
real*8 ei

do i=1,no
  ei = OE(i)
  do a=1,nv
    t1n(a,i) = t1n(a,i)/(ei-OE(no+a))
  end do
end do

return

end subroutine T1_div
