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
! The call to lowdin gives me a transforation A1, and a transformed
! S matrix

implicit real*8(A-H,O-Z)
real*8 Matrix(nDim*nDim), TMatrix(nDim*nDim), Matrix0(nDim*nDim),Temp(nDim*nDim)
integer iCenter(nDim), iType(nDim)
#include "real.fh"

k = 0
do i=1,nDim
  do j=1,nDim
    k = k+1
    !write(6,*) iCenter(i),iCenter(j), Matrix(k)
    if (iCenter(i) /= iCenter(j)) then
      Matrix(k) = Zero
    end if
  end do
end do
call GramSchmidt(Matrix,TMatrix,nDim,iType,iCenter,0)
call dcopy_(nDim**2,Matrix0,1,Matrix,1)

! Now apply T1 to original S: S2=T1(T)*S*T1

call DGEMM_('N','N',nDim,nDim,nDim,1.0d0,Matrix,nDim,TMatrix,nDim,0.0d0,Temp,nDim)
call DGEMM_('T','N',nDim,nDim,nDim,1.0d0,TMatrix,nDim,Temp,nDim,0.0d0,Matrix,nDim)

return

end subroutine Step1
