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
     &                   iComp,nOrdOp,nScr1,nScr2,naa,ipSAR,nSAR,
     &                   iShll_a,nPrim_a,ipExp_a,nCntrc_a,ipCff_a,
     &                   iCmp_a,
     &                   iShll_r,nPrim_r,ipExp_r,nCntrc_r,ipCff_r,
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
      External Kernel
      Real*8 DInf(nDInf), A(3)
*
      ipSAR = ip
      nSAR = nPrim_a*nPrim_r * naa
      ip = ip + nSAR
      ipPAR = ip
      ip = ip + 3 * nPrim_a*nPrim_r
      ipZAR = ip
      ip = ip + nPrim_a*nPrim_r
      ipKAR = ip
      ip = ip + nPrim_a*nPrim_r
      ipZIAR = ip
      ip = ip + nPrim_a*nPrim_r
      nInfo = ipExp(jShll+1) - Info
      mArr = nDInf/(nPrim_a*nPrim_r) - nInfo
      If (mArr.lt.1) Then
         Call WarningMessage(2,'One_Int:  mArr < 1 .'
     &            //'Please, increase MOLCAS_MEM.')
         Call Abend()
      EndIf
*
      Call ZXia(DInf(ipZAR),DInf(ipZIAR),nPrim_a,nPrim_r,
     &          DInf(ipExp_a),DInf(ipExp_r))
      Call SetUp1(DInf(ipExp_a),nPrim_a,
     &            DInf(ipExp_r),nPrim_r,
     &            A,A,DInf(ipKAR),DInf(ipPAR),DInf(ipZIAR))
*
      nHer = (2*iAng+2+nOrdOp)/2
      Call Kernel(DInf(ipExp_a),nPrim_a,
     &            DInf(ipExp_r),nPrim_r,
     &            DInf(ipZAR),DInf(ipZIAR),
     &            DInf(ipKAR),DInf(ipPAR),
     &            DInf(ipSAR),nPrim_a*nPrim_r,iComp,
     &            iAng,iAng,A,A,nHer,DInf(ip),mArr,A,nOrdOp)
      ip = ip - 6 * nPrim_a * nPrim_r
*
      ipScrt1 = ip
      ip = ip + nScr1
      ipScrt2 = ip
      ip = ip + nScr2
      Call DGEMM_('T','N',
     &            nPrim_r*naa,nCntrc_a,nPrim_a,
     &            1.0d0,DInf(ipSAR),nPrim_a,
     &            DInf(ipCff_a),nPrim_a,
     &            0.0d0,DInf(ipScrt1),nPrim_r*naa)
      Call DGEMM_('T','N',
     &            naa*nCntrc_a,nCntrc_r,nPrim_r,
     &            1.0d0,DInf(ipScrt1),nPrim_r,
     &            DInf(ipCff_r),nPrim_r,
     &            0.0d0,DInf(ipScrt2),naa*nCntrc_a)
#ifdef _DEBUG_
      Call RecPrt('S_AR in Cartesian',' ',DInf(ipScrt2),
     &            naa,nCntrc_a*nCntrc_r)
#endif
*
*     Transform to spherical Gaussian if needed!
*
      If (Transf(iShll_a).or.Transf(iShll_r)) Then
*
         Call CarSph(DInf(ipScrt2),naa,nCntrc_a*nCntrc_r,
     &               DInf(ipSAR),nScr2,
     &               RSph(ipSph(iAng)),iAng,
     &               Transf(iShll_a),Prjct(iShll_a),
     &               RSph(ipSph(iAng)),iAng,
     &               Transf(iShll_r),Prjct(iShll_r),
     &               DInf(ipScrt1),iCmp_a*iCmp_r)
         Call DCopy_(nCntrc_a*nCntrc_r*iCmp_a*iCmp_r,
     &               DInf(ipScrt1),1,DInf(ipSAR),1)
      Else
         Call DGeTmO(DInf(ipScrt2),naa,naa,nCntrc_a*nCntrc_r,
     &               DInf(ipSAR),nCntrc_a*nCntrc_r)
      End If
*define _DEBUG_
#ifdef _DEBUG_
      Call RecPrt('S_AR in Sphericals',' ',DInf(ipSAR),
     &                  iCmp_a*iCmp_r,nCntrc_a*nCntrc_r)
#endif
      ip = ip - nScr2 - nScr1
      nSAR = nCntrc_a*nCntrc_r * naa
*
      Return
      End
