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

subroutine GetTauHlp2(Tau,T1,dima,adda,no,nv)
! Make Tau for aGrp == bGrp

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dima, adda, no, nv
real(kind=wp) :: Tau(nTri_Elem(dima),no,no), T1(nv,no)
integer(kind=iwp) :: a, ab, b, i, j
real(kind=wp) :: c

do j=1,no
  do i=1,no
    ab = 0
    do a=1,dima
      c = t1(adda+a,i)
      do b=1,a
        ab = ab+1
        Tau(ab,i,j) = Tau(ab,i,j)+c*t1(adda+b,j)
      end do
    end do
  end do
end do

return

end subroutine GetTauHlp2
