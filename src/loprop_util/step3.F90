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

subroutine Step3(iCenter,SMatrix,nDim,TMatrix,iType)
! Step 3. GS S2 ->S3

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDim, iCenter(nDim), iType(nDim)
real(kind=wp), intent(inout) :: SMatrix(nDim,nDim)
real(kind=wp), intent(out) :: TMatrix(nDim,nDim)

!lg write(u6,*) 'Step 3', nDim
!lg call RecPrt('T before GS 3',' ',TMatrix,nDim,nDim)
!lg write(u6,*)
call unitmat(TMatrix,nDim)
call GramSchmidt(SMatrix,TMatrix,nDim,iType,iCenter,1)

return

end subroutine Step3
