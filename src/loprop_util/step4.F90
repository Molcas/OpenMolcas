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

implicit real*8(A-H,O-Z)
real*8 SMatrix(nDim*nDim), TMatrix(nDim*nDim)
integer iType(nDim)
#include "real.fh"

!lg  write(6,*) 'Step 4', nDim
!lg  call RecPrt('T before LW 4',' ',TMatrix,nDim,nDim)
!lg  call RecPrt('S in step4 ',' ',SMatrix,nDim,nDim)
!lg  write(6,*)
k = 0
do i=1,nDim
  do j=1,nDim
    k = k+1
    if ((i /= j) .and. (iType(i) /= iType(j))) then
      SMatrix(k) = Zero
    end if
  end do
end do
!lg call RecPrt('S before LW 4',' ',SMatrix,nDim,nDim)

call dcopy_(nDim**2,[Zero],0,TMatrix,1)
call dcopy_(nDim,[One],0,TMatrix,nDim+1)
call Lowdin_LP(SMatrix,TMatrix,nDim)

return

end subroutine Step4
