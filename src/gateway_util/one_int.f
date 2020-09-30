************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine One_Int(Kernel,Array,nArray,A,iAng,
     &                   iComp,nOrdOp,
     &                   Scr1,nScr1,Scr2,nScr2,naa,SAR,nSAR,
     &                   iShll_a,nPrim_a,Exp_a,nCntrc_a,Cff_a,
     &                   iCmp_a,
     &                   iShll_r,nPrim_r,Exp_r,nCntrc_r,Cff_r,
     &                   iCmp_r)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
      use Basis_Info
      use Her_RW
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
#include "stdalloc.fh"
      External Kernel
      Real*8 Array(nArray)
      Real*8, Intent(Out):: SAR(nSAR)
      Real*8, Intent(In):: A(3)
      Real*8, Intent(In):: Exp_a(nPrim_a), Exp_r(nPrim_r)
      Real*8, Intent(In):: Cff_a(nPrim_a,nCntrc_a)
      Real*8, Intent(In):: Cff_r(nPrim_r,nCntrc_r)
      Real*8 Scr1(nScr1), Scr2(nScr2)
      Real*8, Allocatable:: ZAR(:), ZIAR(:), KAR(:), PAR(:,:)
      Real*8, Allocatable:: pSAR(:)
*
      mArr = nArray/(nPrim_a*nPrim_r)
      Call mma_allocate(ZAR,nPrim_a*nPrim_r,Label='ZAR')
      Call mma_allocate(ZIAR,nPrim_a*nPrim_r,Label='ZIAR')
      Call mma_allocate(KAR,nPrim_a*nPrim_r,Label='KAR')
      Call mma_allocate(PAR,nPrim_a*nPrim_r,3,Label='PAR')
*
      mSAR = nPrim_a*nPrim_r * naa
      Call mma_allocate(pSAR,mSAR,Label='pSAR')
*
      Call ZXia(ZAR,ZIAR,nPrim_a,nPrim_r,Exp_a,Exp_r)
      Call SetUp1(Exp_a,nPrim_a,Exp_r,nPrim_r,A,A,KAR,PAR,ZIAR)
*
      nHer = (2*iAng+2+nOrdOp)/2
      Call Kernel(Exp_a,nPrim_a,Exp_r,nPrim_r,
     &            ZAR,ZIAR,KAR,PAR,
     &            pSAR,nPrim_a*nPrim_r,iComp,
     &            iAng,iAng,A,A,nHer,Array,mArr,A,nOrdOp)
*
      Call mma_deallocate(ZAR)
      Call mma_deallocate(ZIAR)
      Call mma_deallocate(KAR)
      Call mma_deallocate(PAR)
*
      Call DGEMM_('T','N',
     &            nPrim_r*naa,nCntrc_a,nPrim_a,
     &            1.0d0,pSAR,nPrim_a,
     &                  Cff_a,nPrim_a,
     &            0.0d0,Scr1,nPrim_r*naa)
      Call DGEMM_('T','N',
     &            naa*nCntrc_a,nCntrc_r,nPrim_r,
     &            1.0d0,Scr1,nPrim_r,
     &                  Cff_r,nPrim_r,
     &            0.0d0,Scr2,naa*nCntrc_a)
#ifdef _DEBUG_
      Call RecPrt('S_AR in Cartesian',' ',Scr2,
     &            naa,nCntrc_a*nCntrc_r)
#endif
*
*     Transform to spherical Gaussian if needed!
*
      If (Shells(iShll_a)%Transf.or.Shells(iShll_r)%Transf) Then
*
         Call CarSph(Scr2,naa,nCntrc_a*nCntrc_r,
     &               pSAR,nScr2,
     &               RSph(ipSph(iAng)),iAng,
     &               Shells(iShll_a)%Transf,
     &               Shells(iShll_a)%Prjct,
     &               RSph(ipSph(iAng)),iAng,
     &               Shells(iShll_r)%Transf,
     &               Shells(iShll_r)%Prjct,
     &               SAR,iCmp_a*iCmp_r)
      Else
         Call DGeTmO(Scr2,naa,naa,nCntrc_a*nCntrc_r,
     &               SAR,nCntrc_a*nCntrc_r)
      End If
      Call mma_deallocate(pSAR)
*define _DEBUG_
#ifdef _DEBUG_
      Call RecPrt('S_AR in Sphericals',' ',SAR,
     &                  iCmp_a*iCmp_r,nCntrc_a*nCntrc_r)
#endif
*
      Return
      End
