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
      Subroutine Do_Int_CASDFT2(Weights,nGrid,mAO,
     &                         TabSO,nTabSO,TabMO,nTabMO,
     &                         nD,FckInt,nFckInt,nFckDim,
     &                         RhoI,RhoA,mRho,dF_dRho,ndF_dRho,
     &                         dFdP2ontop,ndFdP2ontop,
     &                         TmpPUVX,nTmpPUVX)
      Implicit Real*8 (A-H,O-Z)
#include "nq_info.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "debug.fh"
#include "nsd.fh"
#include "setup.fh"
#include "nq_index.fh"
      Real*8 Weights(nGrid),
     &       TabSO(mAO,nGrid,nTabSO),
     &       TabMO(mAO,nGrid,nTabMO),
     &       FckInt(nFckInt,nFckDim),
     &       dF_dRho(ndF_dRho,nGrid),
     &       dFdP2ontop(ndFdP2ontop,nGrid),
     &       RhoI(mRho,nGrid),RhoA(mRho,nGrid),
     &       TmpPUVX(nTmpPUVX)
      Integer iOff_Bas(0:7)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      jOff_=0
      Do iIrrep = 0, mIrrep-1
         iOff_Bas(iIrrep)=jOff_
         jOff_=jOff_+mBas(iIrrep)
      End Do
*
      iOff=0
      Do iIrrep = 0, mIrrep - 1
         Do iBas_ = 1, mBas(iIrrep)
            iBas = iBas_ + iOff_Bas(iIrrep)
            Do jBas_ = 1, iBas_
               jBas = jBas_ + iOff_Bas(iIrrep)
               ijInd=iOff+iTri(iBas_,jBas_)
*
*              Sorry Sergey but I do not understand these two different
*              options. mAO is not properly set.
*
               If (mAO.eq.4) Then
                  Do iGrid=1,nGrid
                     OrbProd  = TabSO(1,iGrid,iBas)*TabSO(1,iGrid,jBas)
                     OrbProdx = TabSO(2,iGrid,iBas)*TabSO(1,iGrid,jBas)
     &                        + TabSO(1,iGrid,iBas)*TabSO(2,iGrid,jBas)
                     OrbPrody = TabSO(3,iGrid,iBas)*TabSO(1,iGrid,jBas)
     &                        + TabSO(1,iGrid,iBas)*TabSO(3,iGrid,jBas)
                     OrbProdz = TabSO(4,iGrid,iBas)*TabSO(1,iGrid,jBas)
     &                        + TabSO(1,iGrid,iBas)*TabSO(4,iGrid,jBas)
*
                     FckInt(ijInd,1) = FckInt(ijInd,1)
     &                               + Weights(iGrid)
     &                               *(OrbProd* dF_dRho(ipRa,iGrid)
     &                                +OrbProdx*dF_dRho(ipdRxa,iGrid)
     &                                +OrbPrody*dF_dRho(ipdRya,iGrid)
     &                                +OrbProdz*dF_dRho(ipdRza,iGrid))
                     FckInt(ijInd,2) = FckInt(ijInd,2)
     &                               + Weights(iGrid)
     &                               *(OrbProd* dF_dRho(ipRb,iGrid)
     &                                +OrbProdx*dF_dRho(ipdRxb,iGrid)
     &                                +OrbPrody*dF_dRho(ipdRyb,iGrid)
     &                                +OrbProdz*dF_dRho(ipdRzb,iGrid))
*
                  End Do ! iGrid
               Else
                  Do iGrid=1,nGrid
                     OrbProd  = TabSO(1,iGrid,iBas)*TabSO(1,iGrid,jBas)
                     FckInt(ijInd,1) = FckInt(ijInd,1)
     &                               + Weights(iGrid)
     &                               * OrbProd* dF_dRho(ipRa,iGrid)
                     FckInt(ijInd,2) = FckInt(ijInd,2)
     &                               + Weights(iGrid)
     &                               * OrbProd* dF_dRho(ipRb,iGrid)
*
                  End Do ! iGrid
               End If
*
            End Do  ! jBas
         End Do ! iBas
         iOff = iOff + (mBas(iIrrep)**2+mBas(iIrrep))/2
      End Do  ! iIrrep
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(TabMO)
         Call Unused_integer(nD)
         Call Unused_real_array(RhoI)
         Call Unused_real_array(RhoA)
         Call Unused_real_array(dFdP2ontop)
         Call Unused_real_array(TmpPUVX)
      End If
      End
