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

subroutine Step1(iCenter,Matrix,nDim,TMatrix,iType,Matrix0,Temp)
! Step 1. LO S0 ->S1
! occupied virtual on the same center
! and orthonormalization of the original AO basis
! The call to lowdin gives me a transformation A1, and a transformed
! S matrix

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDim, iCenter(nDim), iType(nDim)
real(kind=wp), intent(inout) :: Matrix(nDim,nDim), TMatrix(nDim,nDim), Matrix0(nDim,nDim)
real(kind=wp), intent(out) :: Temp(nDim,nDim)
integer(kind=iwp) :: i, j

do i=1,nDim
  do j=1,nDim
    !write(u6,*) iCenter(i),iCenter(j), Matrix(j,i)
    if (iCenter(i) /= iCenter(j)) then
      Matrix(j,i) = Zero
    end if
  end do
end do
call GramSchmidt(Matrix,TMatrix,nDim,iType,iCenter,0)
Matrix(:,:) = Matrix0(:,:)

! Now apply T1 to original S: S2=T1(T)*S*T1

call DGEMM_('N','N',nDim,nDim,nDim,One,Matrix,nDim,TMatrix,nDim,Zero,Temp,nDim)
call DGEMM_('T','N',nDim,nDim,nDim,One,TMatrix,nDim,Temp,nDim,Zero,Matrix,nDim)

return

end subroutine Step1
