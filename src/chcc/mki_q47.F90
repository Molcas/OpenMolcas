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

subroutine MkI_Q47(Va,V,dimb,dima,no)
! this routine does:
! Create Va(B',o_A,o_B,A')  needed in step Q48
! from following available array (permuted as given):
! V1(B',o_A,o_B,A') = (A',o_A|B',o_B)
! N.B.
!
! N.B. Kvajt odflaknute, aj koment k rutine odflaknuty

use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dima, dimb, no
real(kind=wp) :: Va(dimb,no,no,dima), V(dimb,no,no,dima)
integer(kind=iwp) :: a, b, i, j

do a=1,dima
  do j=1,no
    do i=1,no
      do b=1,dimb
        Va(b,i,j,a) = Two*V(b,j,i,a)-V(b,i,j,a)
      end do
    end do
  end do
end do

return

end subroutine MkI_Q47
