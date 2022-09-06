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

implicit real*8(A-H,O-Z)
real*8 U_k(nU_k), U_k_New(nU_k_New)
integer iSO_ab(2,nU_k)

do k=1,nU_k
  i = iSO_ab(1,k)
  j = iSO_ab(2,k)
  ij = i*(i-1)/2+j
  if (i == j) then
    U_k_New(ij) = U_k(k)
  else
    U_k_New(ij) = 0.5d0*U_k(k)
  end if
end do

return

end subroutine ReMap_U_k
