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

subroutine Gen_Prim_Density_Matrix(nBas,nPrim,D_p,nOcOb,oNum,oCof)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nBas, nPrim, nOcOb
real(kind=wp), intent(out) :: D_p(nPrim*(nPrim+1)/2)
real(kind=wp), intent(in) :: oNum(nBas), oCof(nBas,nPrim)
integer(kind=iwp) :: i, k, l

D_p(:) = Zero
do k=1,nPrim
  do l=1,k
    do i=1,nOcOb
      D_p(k*(k-1)/2+l) = D_p(k*(k-1)/2+l)+oNum(i)*oCof(i,k)*oCof(i,l)
    end do
  end do
end do

return

end subroutine Gen_Prim_Density_Matrix
