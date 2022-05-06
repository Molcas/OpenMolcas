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
! Copyright (C) 1991,2021,2022, Roland Lindh                           *
!***********************************************************************

subroutine Do_NIntX()

use nq_Grid, only: Grid_AO, TabAO
use nq_Grid, only: AOInt => Dens_AO
use nq_Grid, only: iBfn_Index
use nq_Info

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "stdalloc.fh"
real*8, allocatable :: A1(:,:), A2(:,:), A_tri(:)
real*8, allocatable :: A3(:,:,:), A4(:,:,:)

!                                                                      *
!***********************************************************************
!                                                                      *
mGrid = size(TabAO,2)
nBfn = size(iBfn_Index,2)
nD = size(Grid_AO,4)

!#define _ANALYSIS_
#ifdef _ANALYSIS_
mAO = size(TabAO,1)
nFn = size(Grid_AO,1)
write(6,*)
write(6,*) ' Analysing Grid_AO'
Thr = 1.0D-14
do iD=1,nD
  do iFn=1,nFn
    lBfn = 0
    Total = 0.0d0
    do iBfn=1,nBfn
      lGrid = 0
      do iGrid=1,mGrid
        if (abs(Grid_AO(iFn,iGrid,iBfn,iD)) < Thr) lGrid = lGrid+1
      end do
      if (lGrid == mGrid) lBfn = lBfn+1
      Total = Total+dble(lGrid)/dble(mGrid)
    end do
    Total = Total/dble(nBfn)
    write(6,*) 'Sparsity analysis, iD, iFn',iD,iFn
    write(6,*) ' Total sparsity in %:',1.0d2*Total
    write(6,*) ' Complete Bfn sparsity in %:',1.0d2*dble(lBfn)/dble(nBfn)
    write(6,*)
  end do
end do
write(6,*)
write(6,*) ' Analysing TabAO'
do iAO=1,mAO
  lBfn = 0
  Total = 0.0d0
  do iBfn=1,nBfn
    lGrid = 0
    do iGrid=1,mGrid
      if (abs(TabAO(iAO,iGrid,iBfn)) < Thr) lGrid = lGrid+1
    end do
    if (lGrid == mGrid) lBfn = lBfn+1
    Total = Total+dble(lGrid)/dble(mGrid)
  end do
  Total = Total/dble(nBfn)
  write(6,*) 'Sparsity analysis, iAO',iAO
  write(6,*) ' Total sparsity in %:',1.0d2*Total
  write(6,*) ' Complete Bfn sparsity in %:',.0d2*dble(lBfn)/dble(nBfn)
  write(6,*)
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(A1,mGrid,nBfn,Label='A1')
call mma_allocate(A2,mGrid,nBfn,Label='A2')
!                                                                      *
!***********************************************************************
!                                                                      *
select case (Functional_type)

  case (LDA_type)
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    call mma_Allocate(A_tri,nBfn*(nBfn+1)/2,Label='A_tri')
    AOInt(:,:,:) = Zero
    A2(1:mGrid,1:nBfn) = TabAO(1,1:mGrid,1:nBfn)
    do iD=1,nD
      A1(1:mGrid,1:nBfn) = Grid_AO(1,1:mGrid,1:nBfn,iD)
      call DGEMM_Tri('T','N',nBfn,nBfn,mGrid,One,A1,mGrid,A2,mGrid,Zero,A_Tri,nBfn)
      call Sym_Dist()
    end do
    call mma_deAllocate(A_tri)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  case (GGA_type)
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    A2(1:mGrid,1:nBfn) = TabAO(1,1:mGrid,1:nBfn)
    do iD=1,nD
      A1(1:mGrid,1:nBfn) = Grid_AO(1,1:mGrid,1:nBfn,iD)
      call DGEMM_('T','N',nBfn,nBfn,mGrid,One,A1,mGrid,A2,mGrid,Zero,AOInt(1,1,iD),nBfn)
      call Symmetrize()
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  case (meta_GGA_type1,meta_GGA_type2)
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    A2(1:mGrid,1:nBfn) = TabAO(1,1:mGrid,1:nBfn)
    do iD=1,nD
      A1(1:mGrid,1:nBfn) = Grid_AO(1,1:mGrid,1:nBfn,iD)
      call DGEMM_('T','N',nBfn,nBfn,mGrid,One,A1,mGrid,A2,mGrid,Zero,AOInt(1,1,iD),nBfn)
      call Symmetrize()
    end do

    call mma_allocate(A3,3,mGrid,nBfn,Label='A1')
    call mma_allocate(A4,3,mGrid,nBfn,Label='A2')

    call mma_Allocate(A_tri,nBfn*(nBfn+1)/2,Label='A_tri')
    A4(1:3,1:mGrid,1:nBfn) = TabAO(2:4,1:mGrid,1:nBfn)
    do iD=1,nD
      A3(1:3,1:mGrid,1:nBfn) = Grid_AO(2:4,1:mGrid,1:nBfn,iD)
      call DGEMM_Tri('T','N',nBfn,nBfn,3*mGrid,One,A3,3*mGrid,A4,3*mGrid,Zero,A_Tri,nBfn)
      call Sym_Dist()
    end do
    call mma_deallocate(A3)
    call mma_deallocate(A4)
    call mma_deAllocate(A_tri)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  case default
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    write(6,*) 'DFT_Int: Illegal functional type!'
    write(6,*) Functional_type
    call Abend()
    !                                                                  *
    !*******************************************************************
    !                                                                  *
end select
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(A1)
call mma_deallocate(A2)
!                                                                      *
!***********************************************************************
!                                                                      *
!#define _ANALYSIS_
#ifdef _ANALYSIS_
write(6,*)
write(6,*) ' Analysing AOInt'
Thr = 1.0D-14
do iD=1,nD
  lBfn = 0
  do iBfn=1,nBfn
    do jBfn=1,nBfn
      if (abs(AOInt(iBfn,jBfn,iD)) < Thr) lBfn = lBfn+1
    end do
  end do
  Total = dble(lBfn)/dble(nBfn**2)
  write(6,*) 'Sparsity analysis, iD',iD
  write(6,*) ' Total parcity in %:',1.0d2*Total
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
contains

subroutine Sym_Dist()

  integer iBfn, jBfn, ijBfn

  ijBfn = 0
  do iBfn=1,nBfn
    do jBfn=1,iBfn-1
      ijBfn = ijBfn+1
      AOInt_Sym = A_tri(ijBfn)
      AOInt(iBfn,jBfn,iD) = AOInt(iBfn,jBfn,iD)+AOInt_Sym
      AOInt(jBfn,iBfn,iD) = AOInt(jBfn,iBfn,iD)+AOInt_Sym
    end do
    ijBfn = ijBfn+1
    AOInt_Sym = A_tri(ijBfn)
    AOInt(iBfn,iBfn,iD) = AOInt(iBfn,iBfn,iD)+AOInt_Sym
  end do

end subroutine Sym_Dist

subroutine Symmetrize()

  integer iBfn, jBfn

  do iBfn=1,nBfn
    do jBfn=1,iBfn
      AOInt_Sym = AOInt(iBfn,jBfn,iD)+AOInt(jBfn,iBfn,iD)
      AOInt(iBfn,jBfn,iD) = AOInt_Sym
      AOInt(jBfn,iBfn,iD) = AOInt_Sym
    end do
  end do

end subroutine Symmetrize

end subroutine Do_NIntX
