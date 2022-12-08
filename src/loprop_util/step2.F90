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

subroutine Step2(iMatrix,SMatrix,nDim,TMatrix,iType,SMatrix_Save,Temp)
! Step 2. LO S0 ->S1
! occupied virtual on the same center
! and orthonormalization of the original AO basis
! The call to lowdin_lp gives me a transformation A1, and a
! transformed S matrix

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDim, iMatrix(nDim), iType(nDim)
real(kind=wp), intent(inout) :: SMatrix(nDim,nDim), SMatrix_Save(nDim,nDim)
real(kind=wp), intent(out) :: TMatrix(nDim,nDim), Temp(nDim,nDim)
integer(kind=iwp) :: i, j

!lg write(u6,*) 'Step 2', nDim
!lg call RecPrt('Save before LW 2',' ',SMatrix_Save,nDim,nDim)
!lg write(u6,*)
do i=1,nDim
  do j=1,nDim
    !write(u6,*) iMatrix(i),iMatrix(j),SMatrix(j,i)
    if ((iMatrix(i) /= iMatrix(j)) .and. (iType(i) /= iType(j))) then
      SMatrix(j,i) = Zero
    end if
  end do
end do
!lg call RecPrt('SMatrix before LW 2',' ',SMatrix,nDim,nDim)
!lg call RecPrt('SMatrix_Save before LW 2',' ',SMatrix_Save,nDim,nDim)

call unitmat(TMatrix,nDim)
call Lowdin_LP(SMatrix,TMatrix,nDim)
! Pick up S2
!lg call RecPrt('SMatrix after LW 2',' ',SMatrix,nDim,nDim)
!lg call RecPrt('TMatrix after LW 2',' ',TMatrix,nDim,nDim)
SMatrix(:,:) = SMatrix_Save(:,:)

! Now apply T2 to S2:  S3=T2(T)*S2*T2

call DGEMM_('N','N',nDim,nDim,nDim,One,SMatrix,nDim,TMatrix,nDim,Zero,Temp,nDim)
call DGEMM_('T','N',nDim,nDim,nDim,One,TMatrix,nDim,Temp,nDim,Zero,SMatrix,nDim)
!call RecPrt('S3',' ',SMatrix,nBas(1),nBas(1))
SMatrix_Save(:,:) = SMatrix(:,:)

return

end subroutine Step2
