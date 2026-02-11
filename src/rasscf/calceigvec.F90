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
      Subroutine CalcEigVec(Matrix,NDIM,EigVec)
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None

!*****Input
      INTEGER NDim
!*****Input & Output
      Real*8,DIMENSION(NDIM,NDIM)::Matrix,EigVec
!*****Calculating rotation matrix
      Real*8, Allocatable:: Mat(:),Val(:,:), Scr(:)
      INTEGER NScr,INFO
      Real*8,DIMENSION(2)::WGRONK
!*****Auxiliary quantities
      INTEGER NElem ! NElem=NDim**2
      INTEGER IRow,ICol,IRIC,NI
      LOGICAL UseJacob

      UseJacob=.true.
      EigVec(:,:)=0.0d0

      IF(UseJacob) THEN
       NElem=NDim*(NDim+1)/2
       Call mma_allocate(Mat,nElem,Label='Mat')
       Call mma_allocate(Val,nDim,nDim,Label='Val')
       IRIC=0
       DO IRow=1,NDIM
        Do ICol=1,IRow
         IRIC=IRIC+1
         Mat(IRIC)=Matrix(IRow,ICol)
        End Do
       END DO
       Val(:,:)=0.0D0
       DO NI=1,NDim
        Val(NI,NI)=1.0D0
       END DO
!       write(6,*)'eigenvector matrix before diag'
!       CALL RECPRT(' ',' ',Val,NDIM,NDIM)
!       write(6,*)'matrix to be diagonalized'
!       CALL TriPrt(' ',' ',Mat,NDIM)
       CALL JACOB(Mat,Val,NDim,NDim)
!       write(6,*)'eigenvector matrix'
!       CALL RECPRT(' ',' ',Val,NDIM,NDIM)
!       DO IRow=1,NDIM
!         write(6,*) (EigVec(IRow,ICol), ICol=1,NDim)
!       END DO
       DO ICol=1,NDIM
        Do IRow=1,NDIM
         EigVec(IRow,ICol)=Val(iCol,iRow)
        End Do
       END DO
       Call mma_deallocate(Val)
       Call mma_deallocate(Mat)

      ELSE
       NElem=NDim**2
       Call mma_allocate(Mat,nElem,Label='Mat')
       Call mma_allocate(Val,nDim,nDim,Label='Val')
       DO ICol=1,NDIM
        Do IRow=1,NDIM
         Mat((ICol-1)*NDIM+IRow)=Matrix(IRow,ICol)
        End Do
       END DO
       Val(:,:)=0.0D0
       Call Dsyev_('V','U',NDIM,Mat,NDIM,Val,WGRONK,-1,INFO)
       NScr=Int(WGRONK(1))
       Call mma_allocate(Scr,nScr,Label='Scr')
       Call Dsyev_('V','U',NDIM,Mat,NDIM,Val,Scr,NScr,INFO)
       DO ICol=1,NDIM
        Do IRow=1,NDIM
         EigVec(IRow,ICol)=Mat((IRow-1)*NDIM+ICol)
        End Do
       END DO
       Call mma_deallocate(Scr)
       Call mma_deallocate(Val)
       Call mma_deallocate(Mat)
       END IF

      End Subroutine CalcEigVec
