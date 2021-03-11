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

use Slapaf_Info, only: nStab, Coor
use Slapaf_Parameters, only: iRow_c, nLambda, iter
implicit none
!***********************************************************************
! subroutine to get the T matrix that defines the constrained and      *
! unconstrained subspaces.                                             *
!***********************************************************************
#include "real.fh"
#include "stdalloc.fh"
integer, intent(In) :: mInt
real*8, intent(InOut) :: TMx(mInt,mInt)
integer Lambda1, Lambda2
logical Invert
real*8, allocatable, dimension(:,:) :: C1, C2, CT

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
call dCopy_(mInt*Lambda1,C1,1,CT(1,1),1)
call dCopy_(mInt*Lambda2,C2,1,CT(1,Lambda1+1),1)
call FZero(CT(1,nLambda+1),mInt*(mInt-nLambda))
if (nLambda > 0) then
  call GS(CT,nLambda,TMx,mInt,.false.,.true.)
else
  call FZero(TMx,mInt**2)
  call dCopy_(mInt,[One],0,TMx,mInt+1)
end if

! If NG constraints are to be inverted, combine the complement
! with the global constraints instead

call Qpg_iScalar('Invert constraints',Invert)
if (Invert) call Get_lScalar('Invert constraints',Invert)
if (Invert) then
  Lambda2 = mInt-nLambda
  call dCopy_(mInt*Lambda1,C1,1,CT(1,1),1)
  call dCopy_(mInt*Lambda2,TMx(1,nLambda+1),1,CT(1,Lambda1+1),1)
  nLambda = Lambda1+Lambda2
  call FZero(CT(1,nLambda+1),mInt*(mInt-nLambda))
  if (nLambda > 0) then
    call GS(CT,nLambda,TMx,mInt,.false.,.true.)
  else
    call FZero(TMx,mInt**2)
    call dCopy_(mInt,[One],0,TMx,mInt+1)
  end if
end if

call mma_Deallocate(C1)
call mma_Deallocate(C2)
call mma_Deallocate(CT)

return

end subroutine TMatrix
