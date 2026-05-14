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
use MCLR_Data, only: nDens, nNA
use input_mclr, only: nAsh, nOrb, nRoots, nSym, ntAsh, ntBas
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: AXkzx(nDens)
integer(kind=iwp), intent(in) :: NPUVX, IndPUVX(ntBas,ntAsh,ntAsh,ntAsh)
real(kind=wp), intent(in) :: GDMat(nTri_Elem(nRoots),nnA,nnA), PUVX(NPUVX), zx(nTri_Elem(nRoots-1))
integer(kind=iwp) :: iKK, iKL, iKL2, iLL, iSym, K, L, Off_Act(nSym), Off_Orb(nSym), p
real(kind=wp), allocatable :: AXktmp(:), DKL1(:,:), DKL2(:,:)

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
