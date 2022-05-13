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

use nq_Grid, only: Dens_AO, Grid_AO, iBfn_Index, TabAO
use nq_Info, only: Functional_type, GGA_type, LDA_type, meta_GGA_type1, meta_GGA_type2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: iD, mGrid, nBfn, nD
real(kind=wp) :: AOInt_Sym
real(kind=wp), allocatable :: A1(:,:), A2(:,:), A3(:,:,:), A4(:,:,:), A_tri(:)
!#define _ANALYSIS_
#ifdef _ANALYSIS_
real(kind=wp), parameter :: Thr = 1.0e-14_wp
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
mGrid = size(TabAO,2)
nBfn = size(iBfn_Index,2)
nD = size(Grid_AO,4)

#ifdef _ANALYSIS_
mAO = size(TabAO,1)
nFn = size(Grid_AO,1)
write(u6,*)
write(u6,*) ' Analysing Grid_AO'
do iD=1,nD
  do iFn=1,nFn
    lBfn = 0
    Total = Zero
    do iBfn=1,nBfn
      lGrid = 0
      do iGrid=1,mGrid
        if (abs(Grid_AO(iFn,iGrid,iBfn,iD)) < Thr) lGrid = lGrid+1
      end do
      if (lGrid == mGrid) lBfn = lBfn+1
      Total = Total+real(lGrid,kind=wp)/real(mGrid,kind=wp)
    end do
    Total = Total/real(nBfn,kind=wp)
    write(u6,*) 'Sparsity analysis, iD, iFn',iD,iFn
    write(u6,*) ' Total sparsity in %:',100.0_wp*Total
    write(u6,*) ' Complete Bfn sparsity in %:',100.0_wp*real(lBfn,kind=wp)/real(nBfn,kind=wp)
    write(u6,*)
  end do
end do
write(u6,*)
write(u6,*) ' Analysing TabAO'
do iAO=1,mAO
  lBfn = 0
  Total = Zero
  do iBfn=1,nBfn
    lGrid = 0
    do iGrid=1,mGrid
      if (abs(TabAO(iAO,iGrid,iBfn)) < Thr) lGrid = lGrid+1
    end do
    if (lGrid == mGrid) lBfn = lBfn+1
    Total = Total+real(lGrid,kind=wp)/real(mGrid,kind=wp)
  end do
  Total = Total/real(nBfn,kind=wp)
  write(u6,*) 'Sparsity analysis, iAO',iAO
  write(u6,*) ' Total sparsity in %:',100.0_wp*Total
  write(u6,*) ' Complete Bfn sparsity in %:',100.0_wp*real(lBfn,kind=wp)/real(nBfn,kind=wp)
  write(u6,*)
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
    Dens_AO(:,:,:) = Zero
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
      call DGEMM_('T','N',nBfn,nBfn,mGrid,One,A1,mGrid,A2,mGrid,Zero,Dens_AO(1,1,iD),nBfn)
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
      call DGEMM_('T','N',nBfn,nBfn,mGrid,One,A1,mGrid,A2,mGrid,Zero,Dens_AO(1,1,iD),nBfn)
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
    write(u6,*) 'DFT_Int: Illegal functional type!'
    write(u6,*) Functional_type
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
write(u6,*)
write(u6,*) ' Analysing Dens_AO'
do iD=1,nD
  lBfn = 0
  do iBfn=1,nBfn
    do jBfn=1,nBfn
      if (abs(Dens_AO(iBfn,jBfn,iD)) < Thr) lBfn = lBfn+1
    end do
  end do
  Total = real(lBfn,kind=wp)/real(nBfn**2,kind=wp)
  write(u6,*) 'Sparsity analysis, iD',iD
  write(u6,*) ' Total sparsity in %:',100.0_wp*Total
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
contains

subroutine Sym_Dist()

  integer(kind=iwp) :: iBfn, ijBfn, jBfn

  ijBfn = 0
  do iBfn=1,nBfn
    do jBfn=1,iBfn-1
      ijBfn = ijBfn+1
      AOInt_Sym = A_tri(ijBfn)
      Dens_AO(iBfn,jBfn,iD) = Dens_AO(iBfn,jBfn,iD)+AOInt_Sym
      Dens_AO(jBfn,iBfn,iD) = Dens_AO(jBfn,iBfn,iD)+AOInt_Sym
    end do
    ijBfn = ijBfn+1
    AOInt_Sym = A_tri(ijBfn)
    Dens_AO(iBfn,iBfn,iD) = Dens_AO(iBfn,iBfn,iD)+AOInt_Sym
  end do

end subroutine Sym_Dist

subroutine Symmetrize()

  integer(kind=iwp) :: iBfn, jBfn

  do iBfn=1,nBfn
    do jBfn=1,iBfn
      AOInt_Sym = Dens_AO(iBfn,jBfn,iD)+Dens_AO(jBfn,iBfn,iD)
      Dens_AO(iBfn,jBfn,iD) = AOInt_Sym
      Dens_AO(jBfn,iBfn,iD) = AOInt_Sym
    end do
  end do

end subroutine Symmetrize

end subroutine Do_NIntX
