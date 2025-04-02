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
use stdalloc, only: mma_allocate, mma_deallocate
use MCLR_Data, only: nDens2, nConf1, nNA, nAcPar, nAcPr2
use input_mclr, only: nRoots, ntAsh, ntBas

implicit none
! Output
real*8, dimension(nDens2+6) :: Fock
real*8, dimension(nconf1*nroots) :: CICSF
! Auxiliary Quantities
real*8, dimension(nTri_Elem(nRoots),nnA,nnA) :: GDMat
real*8, dimension(nTri_Elem(nRoots),nTri_Elem(nRoots)) :: W
integer, dimension(ntBas,ntAsh,ntAsh,ntAsh) :: IndPUVX
integer, dimension(ntAsh,ntAsh,ntAsh,ntAsh) :: IndTUVX
real*8, dimension(:), allocatable :: PUVX, R, H, E_Final, AXkzx, AXPzx, AXX, bk, bP, bX, FMO1t, FMO2t, zX
integer NPUVX, NTri

! MEMORY ALLOCATION
call mma_allocate(AXkzx,nDens2)
call mma_allocate(AXPzx,NConf1*nRoots)
call mma_allocate(AXX,nTri_Elem(nRoots-1)**2)
call mma_allocate(R,nRoots**2)
call mma_allocate(H,nRoots**2)
call mma_allocate(E_Final,nRoots)
call mma_allocate(bk,nDens2)
call mma_allocate(bP,nConf1*nRoots)
call mma_allocate(bX,nTri_Elem(nRoots-1))
call mma_allocate(zX,nTri_Elem(nRoots-1))
call Get_PUVXLen(NPUVX)
call mma_allocate(PUVX,NPUVX)
call Get_Ntri(nTri)
call mma_allocate(FMO1t,nRoots*nTri)
NACPAR = nTri_Elem(nnA)
NAcPr2 = nTri_Elem(nacpar)
call mma_allocate(FMO2t,nRoots*nacpr2)
! MAIN COURSE
! First, read results printed in MCPDFT module
call CMSRdMat(H,nRoots,nRoots,'ROT_HAM',7)
call Get_DArray('MS_FINAL_ROT',R,nRoots**2)

call Get_DArray('Last energies',E_Final,nRoots)
call Read_PUVX(PUVX,NPUVX)
call Get_Two_Ind(IndPUVX,IndTUVX)
call GetPDFTFocks(FMO1t,FMO2t,nTri)
! Calculate six additional terms in CMS Lagrangian equaiton
call CMSRHSGDMat(GDMat)
call CalcW(W,GDMAt,PUVX,NPUVX,IndTUVX)

call CalcAXX(AXX,W)

call CalcbXbP_CMSNAC(bX,bP,FMO1t,FMO2t,R,H,E_Final,nTri)

call SolveforzX(zX,AXX,bX)

call CalcAXkzx(AXkzx,GDMat,PUVX,NPUVX,IndPUVX,zx)

call CalcAXPzx(AXPzx,GDMat,PUVX,NPUVX,IndTUVX,W,zx)

call Calcbk_CMSNAC(bk,R,nTri,GDMat,zX)

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
call mma_deallocate(FMO1t)
call mma_deallocate(FMO2t)
call mma_deallocate(PUVX)

end subroutine RHS_CMS_NAC
