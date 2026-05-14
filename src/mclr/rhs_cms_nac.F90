!**********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Paul B Calio                                     *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Based on rhs_cms.f from Jie J. Bao                             *
! ****************************************************************

subroutine RHS_CMS_NAC(Fock,CICSF)

use Index_Functions, only: nTri_Elem
use MCLR_Data, only: nAcPar, nAcPr2, nConf1, nDens, nNA
use input_mclr, only: nRoots, ntAsh, ntBas, ntBtri
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: Fock(nDens), CICSF(nconf1*nroots)
integer(kind=iwp) :: NPUVX
integer(kind=iwp), allocatable :: IndPUVX(:,:,:,:), IndTUVX(:,:,:,:)
real(kind=wp), allocatable :: AXkzx(:), AXPzx(:), AXX(:,:), bk(:), bP(:), bX(:), E_Final(:), FMO1t(:), FMO2t(:), GDMat(:,:,:), &
                              H(:,:), PUVX(:), R(:,:), W(:,:), zX(:)

! MEMORY ALLOCATION
call mma_allocate(AXkzx,nDens)
call mma_allocate(AXPzx,NConf1*nRoots)
call mma_allocate(AXX,nTri_Elem(nRoots-1),nTri_Elem(nRoots-1))
call mma_allocate(R,nRoots,nRoots)
call mma_allocate(H,nRoots,nRoots)
call mma_allocate(E_Final,nRoots)
call mma_allocate(bk,nDens)
call mma_allocate(bP,nConf1*nRoots)
call mma_allocate(bX,nTri_Elem(nRoots-1))
call mma_allocate(zX,nTri_Elem(nRoots-1))
call mma_allocate(IndPUVX,ntBas,ntAsh,ntAsh,ntAsh)
call mma_allocate(IndTUVX,ntAsh,ntAsh,ntAsh,ntAsh)
call mma_allocate(GDMat,nTri_Elem(nRoots),nnA,nnA)
call mma_allocate(W,nTri_Elem(nRoots),nTri_Elem(nRoots))
call Get_PUVXLen(NPUVX)
call mma_allocate(PUVX,NPUVX)
call mma_allocate(FMO1t,nRoots*ntBtri)
NACPAR = nTri_Elem(nnA)
NAcPr2 = nTri_Elem(nacpar)
call mma_allocate(FMO2t,nRoots*nacpr2)
! MAIN COURSE
! First, read results printed in MCPDFT module
call CMSRdMat(H,nRoots,nRoots,'ROT_HAM',7)
call Get_DArray('MS_FINAL_ROT',R,nRoots**2)

call Get_DArray('Last energies',E_Final,nRoots)
call Get_DArray('TwoEIntegral',PUVX,NPUVX)
call Get_Two_Ind(IndPUVX,IndTUVX)
call GetPDFTFocks(FMO1t,FMO2t,ntBtri)
! Calculate six additional terms in CMS Lagrangian equaiton
call CMSRHSGDMat(GDMat)
call CalcW(W,GDMAt,PUVX,NPUVX,IndTUVX)

call CalcAXX(AXX,W)

call CalcbXbP_CMSNAC(bX,bP,FMO1t,FMO2t,R,H,E_Final,ntBtri)

call SolveforzX(zX,AXX,bX)

call CalcAXkzx(AXkzx,GDMat,PUVX,NPUVX,IndPUVX,zx)

call CalcAXPzx(AXPzx,GDMat,PUVX,NPUVX,IndTUVX,W,zx)

call Calcbk_CMSNAC(bk,R,ntBtri,GDMat,zX)

call SolveforRHS(Fock,CICSF,AXkzx,AXPzx,bk,bP)

! MEMORY DEALLOCATION
call mma_deallocate(AXkzx)
call mma_deallocate(AXPzx)
call mma_deallocate(AXX)
call mma_deallocate(R)
call mma_deallocate(H)
call mma_deallocate(E_Final)
call mma_deallocate(bk)
call mma_deallocate(bP)
call mma_deallocate(bX)
call mma_deallocate(zX)
call mma_deallocate(IndPUVX)
call mma_deallocate(IndTUVX)
call mma_deallocate(GDMat)
call mma_deallocate(W)
call mma_deallocate(FMO1t)
call mma_deallocate(FMO2t)
call mma_deallocate(PUVX)

end subroutine RHS_CMS_NAC
