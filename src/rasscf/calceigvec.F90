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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************

subroutine CalcEigVec(Matrix,NDIM,EigVec)

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: NDim
real(kind=wp) :: Matrix(NDIM,NDIM), EigVec(NDIM,NDIM)
integer(kind=iwp) :: ICol, INFO, IRIC, IRow, NElem, NScr
real(kind=wp) :: WGRONK(2)
logical(kind=iwp) :: UseJacob
real(kind=wp), allocatable :: Mat(:), Scr(:), Val(:,:)

UseJacob = .true.
EigVec(:,:) = Zero

if (UseJacob) then
  NElem = nTri_Elem(NDim)
  call mma_allocate(Mat,nElem,Label='Mat')
  call mma_allocate(Val,nDim,nDim,Label='Val')
  IRIC = 0
  do IRow=1,NDIM
    do ICol=1,IRow
      IRIC = IRIC+1
      Mat(IRIC) = Matrix(IRow,ICol)
    end do
  end do
  call unitmat(Val,NDim)
  !write(u6,*) 'eigenvector matrix before diag'
  !call RECPRT(' ',' ',Val,NDIM,NDIM)
  !write(u6,*) 'matrix to be diagonalized'
  !call TriPrt(' ',' ',Mat,NDIM)
  call JACOB(Mat,Val,NDim,NDim)
  !write(u6,*)'eigenvector matrix'
  !call RECPRT(' ',' ',Val,NDIM,NDIM)
  !do IRow=1,NDIM
  !  write(u6,*) (EigVec(IRow,ICol),ICol=1,NDim)
  !end do
  do ICol=1,NDIM
    EigVec(:,ICol) = Val(iCol,:)
  end do
  call mma_deallocate(Val)
  call mma_deallocate(Mat)

else
  NElem = NDim**2
  call mma_allocate(Mat,nElem,Label='Mat')
  call mma_allocate(Val,nDim,nDim,Label='Val')
  do ICol=1,NDIM
    Mat((ICol-1)*NDIM+1:ICol*NDIM) = Matrix(:,ICol)
  end do
  Val(:,:) = Zero
  call Dsyev_('V','U',NDIM,Mat,NDIM,Val,WGRONK,-1,INFO)
  NScr = int(WGRONK(1))
  call mma_allocate(Scr,nScr,Label='Scr')
  call Dsyev_('V','U',NDIM,Mat,NDIM,Val,Scr,NScr,INFO)
  do IRow=1,NDIM
    EigVec(IRow,:) = Mat((IRow-1)*NDIM+1:IRow*NDIM)
  end do
  call mma_deallocate(Scr)
  call mma_deallocate(Val)
  call mma_deallocate(Mat)
end if

end subroutine CalcEigVec
