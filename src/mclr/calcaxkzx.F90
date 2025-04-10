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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************

subroutine CalcAXkzx(AXkzx,GDMat,PUVX,NPUVX,IndPUVX,zx)

use Index_Functions, only: iTri, nTri_Elem
use MCLR_Data, only: nNA, nDens
use input_mclr, only: nRoots, ntBas, ntAsh, nSym, nAsh, nOrb
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero

implicit none
integer NPUVX
real*8, dimension(nTri_Elem(nRoots),nnA,nnA), intent(in) :: GDMat
real*8, dimension(NPUVX), intent(in) :: PUVX
integer, dimension(ntBas,ntAsh,ntAsh,ntAsh), intent(in) :: IndPUVX
real*8, dimension(nTri_Elem(nRoots-1)), intent(in) :: zx
real*8, dimension(nDens), intent(out) :: AXkzx
! Auxiliary Quantities
integer, dimension(nSym) :: Off_Act, Off_Orb
real*8, dimension(:), allocatable :: AXktmp
real*8, dimension(:,:), allocatable :: DKL1, DKL2
integer K, L, iKL, iKL2, iKK, iLL
integer p, iSym

Off_Act(1) = 0
Off_Orb(1) = 0
do ISym=2,nSym
  Off_Act(ISym) = Off_Act(ISym-1)+nAsh(iSym-1)
  Off_Orb(ISym) = Off_Orb(ISym-1)+nOrb(iSym-1)
end do

AXkzx(:) = Zero

call mma_allocate(DKL1,ntAsh,ntAsh)
call mma_allocate(DKL2,ntAsh,ntAsh)
call mma_allocate(AXktmp,nDens)

do K=2,nRoots
  do L=1,K-1
    iKL = iTri(K,L)
    iKK = iTri(K,K)
    iLL = iTri(L,L)
    iKL2 = nTri_Elem(K-2)+L
    do p=1,ntAsh
      DKL1(:,p) = GDMat(IKL,p,1:ntAsh)+GDMat(IKL,1:ntAsh,p)
      DKL2(:,p) = GDMat(IKK,p,1:ntAsh)-GDMat(ILL,p,1:ntAsh)
    end do
    AXktmp(:) = Zero
    call CalcAXk2(AXktmp,DKL1,DKL2,PUVX,NPUVX,IndPUVX,Off_Act,Off_Orb)
    call CalcAXk2(AXktmp,DKL2,DKL1,PUVX,NPUVX,IndPUVX,Off_Act,Off_Orb)
    AXkzx(:) = AXkzx(:)+zx(IKL2)*AXktmp(:)
  end do
end do

call mma_deallocate(DKL1)
call mma_deallocate(DKL2)
call mma_deallocate(Axktmp)

end subroutine CalcAXkzx
