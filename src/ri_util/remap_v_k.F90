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

subroutine ReMap_V_k(iSym,V_k,nV_k,V_k_New,nV_k_New,iSO_ab,ij2K,m_ij2K)

use Index_Functions, only: nTri_Elem
use Constants, only: Half
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: iSym, nV_k, nV_k_New, iSO_ab(2,nV_k), m_ij2K
real(kind=wp), intent(in) :: V_k(nV_k)
real(kind=wp), intent(_OUT_) :: V_k_New(nV_k_New)
integer(kind=iwp), intent(out) :: ij2K(m_ij2K)
integer(kind=iwp) :: i, ij, j, k

if (iSym == 1) then
  do k=1,nV_k
    i = iSO_ab(1,k)
    j = iSO_ab(2,k)
    ij = nTri_Elem(i-1)+j
    if (i == j) then
      V_k_New(ij) = V_k(k)
    else
      V_k_New(ij) = Half*V_k(k)
    end if
    ij2K(ij) = k
  end do

  !write(u6,*) 'Triang <Vk|Vk> : ',ddot_(nV_k_New,V_k_New,1,V_k_New,1)

else

  do k=1,nV_k
    i = iSO_ab(1,k)
    j = iSO_ab(2,k)
    ij = nTri_Elem(i-1)+j
    ij2K(ij) = k
  end do
end if

return

end subroutine ReMap_V_k
