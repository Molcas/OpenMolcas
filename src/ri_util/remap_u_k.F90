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

subroutine ReMap_U_k(U_k,nU_k,U_k_New,nU_k_New,iSO_ab)

use Index_Functions, only: nTri_Elem
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nU_k, nU_k_New, iSO_ab(2,nU_k)
real(kind=wp), intent(in) :: U_k(nU_k)
real(kind=wp), intent(out) :: U_k_New(nU_k_New)
integer(kind=iwp) :: i, ij, j, k

do k=1,nU_k
  i = iSO_ab(1,k)
  j = iSO_ab(2,k)
  ij = nTri_Elem(i-1)+j
  if (i == j) then
    U_k_New(ij) = U_k(k)
  else
    U_k_New(ij) = Half*U_k(k)
  end if
end do

return

end subroutine ReMap_U_k
