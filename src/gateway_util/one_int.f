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
      Subroutine One_Int(Kernel,DInf,nDInf,A,ip,Info,nInfo,jShll,iAng,
     &                   iComp,nOrdOp,
     &                   Scr1,nScr1,Scr2,nScr2,naa,ipSAR,nSAR,
     &                   iShll_a,nPrim_a,Exp_a,nCntrc_a,Cff_a,
     &                   iCmp_a,
     &                   iShll_r,nPrim_r,Exp_r,nCntrc_r,Cff_r,
     &                   iCmp_r)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
      use Her_RW
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
      External Kernel
      Real*8 DInf(nDInf)
      Real*8, Intent(In):: A(3)
      Real*8, Intent(In):: Exp_a(nPrim_a), Exp_r(nPrim_r)
      Real*8, Intent(In):: Cff_a(nPrim_a,nCntrc_a)
      Real*8, Intent(In):: Cff_r(nPrim_r,nCntrc_r)
      Real*8 Scr1(nScr1), Scr2(nScr2)
      Real*8, Allocatable:: ZAR(:), ZIAR(:), KAR(:), PAR(:,:)
*
      Call mma_allocate(ZAR,nPrim_a*nPrim_r,Label='ZAR')
      Call mma_allocate(ZIAR,nPrim_a*nPrim_r,Label='ZIAR')
      Call mma_allocate(KAR,nPrim_a*nPrim_r,Label='KAR')
      Call mma_allocate(PAR,nPrim_a*nPrim_r,3,Label='PAR')
*
      ipSAR = ip
      nSAR = nPrim_a*nPrim_r * naa
      ip = ip + nSAR
      nInfo = ipExp(jShll+1) - Info
      mArr = nDInf/(nPrim_a*nPrim_r) - nInfo
      If (mArr.lt.1) Then
         Call WarningMessage(2,'One_Int:  mArr < 1 .'
     &            //'Please, increase MOLCAS_MEM.')
         Call Abend()
      EndIf
*
      Call ZXia(ZAR,ZIAR,nPrim_a,nPrim_r,
     &          Exp_a,Exp_r)
      Call SetUp1(Exp_a,nPrim_a,
     &            Exp_r,nPrim_r,
     &            A,A,KAR,PAR,ZIAR)
*
      nHer = (2*iAng+2+nOrdOp)/2
      Call Kernel(Exp_a,nPrim_a,
     &            Exp_r,nPrim_r,
     &            ZAR,ZIAR,
     &            KAR,PAR,
     &            DInf(ipSAR),nPrim_a*nPrim_r,iComp,
     &            iAng,iAng,A,A,nHer,DInf(ip),mArr,A,nOrdOp)
*
      Call mma_deallocate(ZAR)
      Call mma_deallocate(ZIAR)
      Call mma_deallocate(KAR)
      Call mma_deallocate(PAR)
*
      Call DGEMM_('T','N',
     &            nPrim_r*naa,nCntrc_a,nPrim_a,
     &            1.0d0,DInf(ipSAR),nPrim_a,
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
      If (Transf(iShll_a).or.Transf(iShll_r)) Then
*
         Call CarSph(Scr2,naa,nCntrc_a*nCntrc_r,
     &               DInf(ipSAR),nScr2,
     &               RSph(ipSph(iAng)),iAng,
     &               Transf(iShll_a),Prjct(iShll_a),
     &               RSph(ipSph(iAng)),iAng,
     &               Transf(iShll_r),Prjct(iShll_r),
     &               Scr1,iCmp_a*iCmp_r)
         Call DCopy_(nCntrc_a*nCntrc_r*iCmp_a*iCmp_r,
     &               Scr1,1,DInf(ipSAR),1)
      Else
         Call DGeTmO(Scr2,naa,naa,nCntrc_a*nCntrc_r,
     &               DInf(ipSAR),nCntrc_a*nCntrc_r)
      End If
*define _DEBUG_
#ifdef _DEBUG_
      Call RecPrt('S_AR in Sphericals',' ',DInf(ipSAR),
     &                  iCmp_a*iCmp_r,nCntrc_a*nCntrc_r)
#endif
      nSAR = nCntrc_a*nCntrc_r * naa
*
      Return
      End
