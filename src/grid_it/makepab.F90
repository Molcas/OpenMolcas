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

subroutine MakePab(cmo,occ,cout,nCMO,nMOs,nIrrep,nBas)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nCMO, nMOs, nIrrep, nBas(0:7)
real(kind=wp), intent(in) :: cmo(nCMO), occ(nMOs)
real(kind=wp), intent(out) :: cout(nMOs)
integer(kind=iwp) :: i, id, id2, iIrr, j, nd, nd2

id = 0
id2 = 0
cout(:) = Zero
do iIrr=0,nIrrep-1
  nd = nBas(iIrr)
  nd2 = nd*nd
  do i=1,nd
    !if (occ(i+id) /= Zero) then
    do j=1,nd
      cout(i+id) = cout(i+id)+occ(i+id)*(cmo(j+i*(nd-1)+id2)**2)
    end do
  end do
  id = id+nd
  id2 = id2+nd2
end do

return

end subroutine MakePab
