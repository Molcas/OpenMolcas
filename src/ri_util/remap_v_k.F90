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

implicit real*8(A-H,O-Z)
real*8 V_k(nV_k), V_k_New(nV_k_New)
integer iSym, iSO_ab(2,nV_k), ij2K(m_ij2K)

if (iSym == 1) then
  do k=1,nV_k
    i = iSO_ab(1,k)
    j = iSO_ab(2,k)
    ij = i*(i-1)/2+j
    if (i == j) then
      V_k_New(ij) = V_k(k)
    else
      V_k_New(ij) = 0.5d0*V_k(k)
    end if
    ij2K(ij) = k
  end do

  !write(6,*) 'Triang <Vk|Vk> : ',ddot_(nV_k_New,V_k_New,1,V_k_New,1)

else

  do k=1,nV_k
    i = iSO_ab(1,k)
    j = iSO_ab(2,k)
    ij = i*(i-1)/2+j
    ij2K(ij) = k
  end do
end if

return

end subroutine ReMap_V_k
