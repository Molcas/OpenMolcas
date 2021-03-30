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

subroutine Get_Prim_Density_Matrix(D,nBas,D_p,nPrim,TM)

use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nBas, nPrim
real(kind=wp), intent(inout) :: D(nBas*(nBas+1)/2)
real(kind=wp), intent(out) :: D_p(nPrim*(nPrim+1)/2)
real(kind=wp), intent(in) :: TM(nPrim,nBas)
integer(kind=iwp) :: i, j, k, l
real(kind=wp) :: Djl, TMij, TMkl, TmpDensity

D(:) = Half*D(:)
do i=1,nPrim
  do k=1,i
    TmpDensity = Zero
    do j=1,nBas
      do l=1,nBas
        TMij = TM(i,j)
        TMkl = TM(k,l)
        if (j < l) then
          Djl = D(l*(l-1)/2+j)
        else
          Djl = D(j*(j-1)/2+l)
        end if
        TmpDensity = TmpDensity+TMij*TMkl*Djl
      end do
    end do
    D_p(i*(i-1)/2+k) = TmpDensity
  end do
end do

return

end subroutine Get_Prim_Density_Matrix
