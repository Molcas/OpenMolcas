!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2015, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine TMatrix(TMx,mInt)
!***********************************************************************
! subroutine to get the T matrix that defines the constrained and      *
! unconstrained subspaces.                                             *
!***********************************************************************

use Slapaf_Info, only: Coor, iRow_c, iter, nLambda, nStab
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mInt
real(kind=wp), intent(inout) :: TMx(mInt,mInt)
integer(kind=iwp) :: Lambda1, Lambda2
logical(kind=iwp) :: Invert
real(kind=wp), allocatable, dimension(:,:) :: C1, C2, CT

! Get the global constraint vectors

call mma_Allocate(C1,mInt,nLambda)
call get_drdq(C1,mInt,nLambda,Lambda1,Iter)

! Get the NG constraint vectors

call Merge_Constraints('UDC.NG','','UDC',nLambda,iRow_c)
call Fix_UDC(iRow_c,nLambda,size(Coor,2),nStab,.true.)
call mma_Allocate(C2,mInt,nLambda)
call get_drdq(C2,mInt,nLambda,Lambda2,Iter)

! Combine both sets of constraints and get the T matrix

nLambda = Lambda1+Lambda2
call mma_Allocate(CT,mInt,mInt)
CT(:,1:Lambda1) = C1(:,:)
CT(:,Lambda1+1:nLambda) = C2(:,:)
CT(:,nLambda+1:) = Zero
if (nLambda > 0) then
  call GS(CT,nLambda,TMx,mInt,.false.,.true.)
else
  call unitmat(TMx,mInt)
end if

! If NG constraints are to be inverted, combine the complement
! with the global constraints instead

call Qpg_iScalar('Invert constraints',Invert)
if (Invert) call Get_lScalar('Invert constraints',Invert)
if (Invert) then
  Lambda2 = mInt-nLambda
  CT(:,1:Lambda1) = C1(:,:)
  CT(:,Lambda1+1:Lambda1+Lambda2) = TMx(:,nLambda+1:nLambda+Lambda2)
  nLambda = Lambda1+Lambda2
  CT(:,nLambda+1:) = Zero
  if (nLambda > 0) then
    call GS(CT,nLambda,TMx,mInt,.false.,.true.)
  else
    call unitmat(TMx,mInt)
  end if
end if

call mma_Deallocate(C1)
call mma_Deallocate(C2)
call mma_Deallocate(CT)

return

end subroutine TMatrix
