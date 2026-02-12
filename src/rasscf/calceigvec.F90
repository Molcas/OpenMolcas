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

use stdalloc, only: mma_allocate, mma_deallocate

!*****Input
integer NDim
!*****Input & Output
real*8, dimension(NDIM,NDIM) :: Matrix, EigVec
!*****Calculating rotation matrix
real*8, allocatable :: Mat(:), Val(:,:), Scr(:)
integer NScr, INFO
real*8, dimension(2) :: WGRONK
!*****Auxiliary quantities
integer NElem ! NElem=NDim**2
integer IRow, ICol, IRIC, NI
logical UseJacob

UseJacob = .true.
EigVec(:,:) = 0.0d0

if (UseJacob) then
  NElem = NDim*(NDim+1)/2
  call mma_allocate(Mat,nElem,Label='Mat')
  call mma_allocate(Val,nDim,nDim,Label='Val')
  IRIC = 0
  do IRow=1,NDIM
    do ICol=1,IRow
      IRIC = IRIC+1
      Mat(IRIC) = Matrix(IRow,ICol)
    end do
  end do
  Val(:,:) = 0.0d0
  do NI=1,NDim
    Val(NI,NI) = 1.0d0
  end do
  !write(6,*) 'eigenvector matrix before diag'
  !call RECPRT(' ',' ',Val,NDIM,NDIM)
  !write(6,*) 'matrix to be diagonalized'
  !call TriPrt(' ',' ',Mat,NDIM)
  call JACOB(Mat,Val,NDim,NDim)
  !write(6,*)'eigenvector matrix'
  !call RECPRT(' ',' ',Val,NDIM,NDIM)
  !do IRow=1,NDIM
  !  write(6,*) (EigVec(IRow,ICol),ICol=1,NDim)
  !end do
  do ICol=1,NDIM
    do IRow=1,NDIM
      EigVec(IRow,ICol) = Val(iCol,iRow)
    end do
  end do
  call mma_deallocate(Val)
  call mma_deallocate(Mat)

else
  NElem = NDim**2
  call mma_allocate(Mat,nElem,Label='Mat')
  call mma_allocate(Val,nDim,nDim,Label='Val')
  do ICol=1,NDIM
    do IRow=1,NDIM
      Mat((ICol-1)*NDIM+IRow) = Matrix(IRow,ICol)
    end do
  end do
  Val(:,:) = 0.0d0
  call Dsyev_('V','U',NDIM,Mat,NDIM,Val,WGRONK,-1,INFO)
  NScr = int(WGRONK(1))
  call mma_allocate(Scr,nScr,Label='Scr')
  call Dsyev_('V','U',NDIM,Mat,NDIM,Val,Scr,NScr,INFO)
  do ICol=1,NDIM
    do IRow=1,NDIM
      EigVec(IRow,ICol) = Mat((IRow-1)*NDIM+ICol)
    end do
  end do
  call mma_deallocate(Scr)
  call mma_deallocate(Val)
  call mma_deallocate(Mat)
end if

end subroutine CalcEigVec
