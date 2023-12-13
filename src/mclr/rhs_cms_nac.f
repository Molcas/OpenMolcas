************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2021, Paul B Calio                                     *
************************************************************************
* ****************************************************************
* history:                                                       *
* Based on rhs_cms.f from Jie J. Bao                             *
* ****************************************************************
      subroutine RHS_CMS_NAC(Fock,CICSF)
      use stdalloc, only : mma_allocate, mma_deallocate
*#include "stdalloc.fh"
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "real.fh"
#include "sa.fh"

******Input
******Output
      Real*8,DIMENSION(nDens2+6)::Fock
      Real*8,DIMENSION(nconf1*nroots)::CICSF
******Auxiliary Quantities
      Real*8,DIMENSION((nRoots+1)*nRoots/2,nnA,nnA)::GDMat
      Real*8,DIMENSION((nRoots+1)*nRoots/2,(nRoots+1)*nRoots/2)::W
      INTEGER,DIMENSION(ntBas,ntAsh,ntAsh,ntAsh)::IndPUVX
      INTEGER,DIMENSION(ntAsh,ntAsh,ntAsh,ntAsh)::IndTUVX

      Real*8,DIMENSION(:),Allocatable::PUVX,R,H,E_Final,AXkzx,AXPzx,AXX,
     &bk,bP,bX,FMO1t,FMO2t,zX
      INTEGER NPUVX,NTri

******MEMORY ALLOCATION
      CALL mma_allocate(AXkzx,nDens2)
      CALL mma_allocate(AXPzx,NConf1*nRoots)
      CALL mma_allocate(AXX,((nRoots-1)*nRoots/2)**2)
      CALL mma_allocate(R,nRoots**2)
      CALL mma_allocate(H,nRoots**2)
      CALL mma_allocate(E_Final,nRoots)
      CALL mma_allocate(bk,nDens2)
      CALL mma_allocate(bP,nConf1*nRoots)
      CALL mma_allocate(bX,(nRoots-1)*nRoots/2)
      CALL mma_allocate(zX,(nRoots-1)*nRoots/2)
      Call Get_PUVXLen(NPUVX)
      CALL mma_allocate(PUVX,NPUVX)
      CALL Get_Ntri(nTri)
      CALL mma_allocate(FMO1t,nRoots*nTri)
      NACPAR=(nnA+1)*nnA/2
      NAcPr2=(nacpar+1)*nacpar/2
      CALL mma_allocate(FMO2t,nRoots*nacpr2)
******MAIN COURSE
******First, read results printed in MCPDFT module
      CALL CMSRdMat(H,nRoots,nRoots,'ROT_HAM',7)
      CALL Get_DArray('MS_FINAL_ROT    ',R,nRoots**2)

      CALL Get_DArray('Last energies',E_Final,nRoots)
      Call Read_PUVX(PUVX,NPUVX)
      CALL Get_Two_Ind(IndPUVX,IndTUVX)
      CALL GetPDFTFocks(FMO1t,FMO2t,nTri)
******Calculate six additional terms in CMS Lagrangian equaiton
      CALL CMSRHSGDMat(GDMat)
      CALL CalcW(W,GDMAt,PUVX,NPUVX,IndTUVX)

      CALL CalcAXX(AXX,W)

      CALL CalcbXbP_CMSNAC(bX,bP,FMO1t,FMO2t,R,H,E_Final,nTri)

      CALL SolveforzX(zX,AXX,bX)

      CALL CalcAXkzx(AXkzx,GDMat,PUVX,NPUVX,IndPUVX,zx)

      CALL CalcAXPzx(AXPzx,GDMat,PUVX,NPUVX,IndTUVX,W,zx)

      CALL Calcbk_CMSNAC(bk,R,nTri,GDMat,zX)

      CALL SolveforRHS(Fock,CICSF,AXkzx,AXPzx,bk,bP)

******MEMORY DEALLOCATION
      CALL mma_deallocate(AXkzx)
      CALL mma_deallocate(AXPzx)
      CALL mma_deallocate(AXX)
      CALL mma_deallocate(R)
      CALL mma_deallocate(H)
      CALL mma_deallocate(E_Final)
      CALL mma_deallocate(bk)
      CALL mma_deallocate(bP)
      CALL mma_deallocate(bX)
      CALL mma_deallocate(zX)
      CALL mma_deallocate(FMO1t)
      CALL mma_deallocate(FMO2t)
      CALL mma_deallocate(PUVX)
      RETURN
      end subroutine
