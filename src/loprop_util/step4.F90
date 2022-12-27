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

subroutine Step4(SMatrix,nDim,TMatrix,iType)
! Step 4. LW S3 ->S4

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDim, iType(nDim)
real(kind=wp), intent(inout) :: SMatrix(nDim,nDim)
real(kind=wp), intent(out) :: TMatrix(nDim,nDim)
integer(kind=iwp) :: i, j

!lg  write(u6,*) 'Step 4', nDim
!lg  call RecPrt('T before LW 4',' ',TMatrix,nDim,nDim)
!lg  call RecPrt('S in step4 ',' ',SMatrix,nDim,nDim)
!lg  write(u6,*)
do i=1,nDim
  do j=1,nDim
    if ((i /= j) .and. (iType(i) /= iType(j))) then
      SMatrix(j,i) = Zero
    end if
  end do
end do
!lg call RecPrt('S before LW 4',' ',SMatrix,nDim,nDim)

call unitmat(TMatrix,nDim)
call Lowdin_LP(SMatrix,TMatrix,nDim)

return

end subroutine Step4
