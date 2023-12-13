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

subroutine set_nnBSF(nSym,nBas,nnBSF,n2BSF)

use Index_Functions, only: nTri_Elem
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(8)
integer(kind=iwp), intent(out) :: nnBSF(8,8), n2BSF(8,8)
integer(kind=iwp) :: i, j

do j=1,nSym
  do i=j,nSym

    if (i == j) then
      nnBSF(i,j) = nTri_Elem(nBas(i))
    else
      nnBSF(i,j) = nBas(i)*nBas(j)
    end if

    nnBSF(j,i) = nnBSF(i,j)

    n2BSF(i,j) = nBas(i)*nBas(j)

    n2BSF(j,i) = n2BSF(i,j)

  end do
end do

return

end subroutine set_nnBSF
