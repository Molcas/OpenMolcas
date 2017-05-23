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
************************************************************************
*                                                                      *
      Subroutine Do_NInt1_d(AOInt,nAOInt,ndF_dRho,dF_dRho,
     &                      Weights,mGrid,Rho,nRho,
     &                      Scr,TabAO1,iCmp,iBas,nGrid_Tot,iSpin,
     &                      mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 AOInt(nAOInt*nAOInt,iSpin), Rho(nRho,mGrid),
     &       dF_dRho(ndF_dRho,mGrid), Weights(mGrid),
     &       TabAO1(mAO,mGrid,iBas*iCmp),
     &       Scr(iSpin*nFn,mGrid,iBas*iCmp)

      nBB = iBas*iBas
      nGrid_Tot=nGrid_Tot+mGrid*iBas*iCmp*iBas*iCmp
*
      If (iSpin.ne.1) Go To 99
*
      Do iCB = 1, iBas*iCmp
         iC = (iCB-1)/iBas + 1
         iB = iCB - (iC-1)*iBas
         iiC = (iC-1)*iCmp + iC
         iiB = (iB-1)*iBas + iB
         ii = (iiC-1)*nBB + iiB
*                                                                      *
************************************************************************
*                                                                      *
         ToAdd1=Zero
         Do iGrid = 1, mGrid
            Tmp = TabAO1(1,iGrid,iCB)
     &                   * dF_dRho(ipR,iGrid) * Weights(iGrid)
            ToAdd1 = ToAdd1 + Tmp*TabAO1(1,iGrid,iCB)
            Scr(1,iGrid,iCB) = Tmp
         End Do
         AOInt(ii,1) = ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         Do jCB = iCB+1, iBas*iCmp
            jC = (jCB-1)/iBas + 1
            jB = jCB - (jC-1)*iBas
            ijC = (jC-1)*iCmp + iC
            ijB = (jB-1)*iBas + iB
            ij = (ijC-1)*nBB + ijB
            jiC= (iC-1)*iCmp + jC
            jiB= (iB-1)*iBas + jB
            ji = (jiC-1)*nBB + jiB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
            Do iGrid = 1, mGrid
               ToAdd1 = ToAdd1 + Scr(1,iGrid,iCB)*TabAO1(1,iGrid,jCB)
            End Do
            AOInt(ij,1) = ToAdd1
            AOInt(ji,1) = ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
*
      End Do
*
      Return
*
 99   Continue
*
      Do iCB = 1, iBas*iCmp
         iC = (iCB-1)/iBas + 1
         iB = iCB - (iC-1)*iBas
         iiC = (iC-1)*iCmp + iC
         iiB = (iB-1)*iBas + ib
         ii = (iiC-1)*nBB + iiB
*                                                                      *
************************************************************************
*                                                                      *
         ToAdd1=Zero
         ToAdd2=Zero
         Do iGrid = 1, mGrid
            Tmp1=TabAO1(1,iGrid,iCB)*dF_dRho(ipRa,iGrid)*Weights(iGrid)
            Tmp2=TabAO1(1,iGrid,iCB)*dF_dRho(ipRb,iGrid)*Weights(iGrid)
            ToAdd1 = ToAdd1 + Tmp1*TabAO1(1,iGrid,iCB)
            ToAdd2 = ToAdd2 + Tmp2*TabAO1(1,iGrid,iCB)
            Scr(1,iGrid,iCB) = Tmp1
            Scr(2,iGrid,iCB) = Tmp2
         End Do
         AOInt(ii,1) = ToAdd1
         AOInt(ii,2) = ToAdd2
*                                                                      *
************************************************************************
*                                                                      *
         Do jCB = iCB+1, iBas*iCmp
            jC = (jCB-1)/iBas + 1
            jB = jCB - (jC-1)*iBas
            ijC = (jC-1)*iCmp + iC
            ijB = (jB-1)*iBas + ib
            ij = (ijC-1)*nBB + ijB
            jiC= (iC-1)*iCmp + jC
            jiB= (iB-1)*iBas + jB
            ji = (jiC-1)*nBB + jiB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
            ToAdd2=Zero
            Do iGrid = 1, mGrid
               ToAdd1 = ToAdd1 + Scr(1,iGrid,iCB)*TabAO1(1,iGrid,jCB)
               ToAdd2 = ToAdd2 + Scr(2,iGrid,iCB)*TabAO1(1,iGrid,jCB)
            End Do
            AOInt(ij,1) = ToAdd1
            AOInt(ji,1) = ToAdd1
            AOInt(ij,2) = ToAdd2
            AOInt(ji,2) = ToAdd2
*                                                                      *
************************************************************************
*                                                                      *
         End Do
*
      End Do
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Rho)
      End
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Subroutine Do_NInt1(AOInt,nAOInt,mGrid,
     &                    TabAO1,iCmp,iBas,
     &                    TabAO2,jCmp,jBas,nGrid_Tot,iSpin,mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 AOInt(nAOInt*nAOInt,iSpin),
     &       TabAO1(iSpin*nFn,mGrid,iBas*iCmp),
     &       TabAO2(mAO,mGrid,jBas*jCmp)

      nBB = iBas*jBas
      nGrid_Tot=nGrid_Tot+mGrid*iBas*iCmp*jBas*jCmp
*
      If (iSpin.ne.1) Go To 99
*
      Do iCB = 1, iBas*iCmp
         iC = (iCB-1)/iBas + 1
         iB = iCB - (iC-1)*iBas
*                                                                      *
************************************************************************
*                                                                      *
         Do jCB = 1, jBas*jCmp
            jC = (jCB-1)/jBas + 1
            jB = jCB - (jC-1)*jBas
            ijC = (jC-1)*iCmp + iC
            ijB = (jB-1)*iBas + ib
            ij = (ijC-1)*nBB + ijB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
            Do iGrid = 1, mGrid
               ToAdd1 = ToAdd1 + TabAO1(1,iGrid,iCB)*TabAO2(1,iGrid,jCB)
            End Do
            AOInt(ij,1) = ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
*
      Return
*
 99   Continue
*
      Do iCB = 1, iBas*iCmp
         iC = (iCB-1)/iBas + 1
         iB = iCB - (iC-1)*iBas
*                                                                      *
************************************************************************
*                                                                      *
         Do jCB = 1, jBas*jCmp
            jC = (jCB-1)/jBas + 1
            jB = jCB - (jC-1)*jBas
            ijC = (jC-1)*iCmp + iC
            ijB = (jB-1)*iBas + ib
            ij = (ijC-1)*nBB + ijB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
            ToAdd2=Zero
            Do iGrid = 1, mGrid
               ToAdd1 = ToAdd1 + TabAO1(1,iGrid,iCB)*TabAO2(1,iGrid,jCB)
               ToAdd2 = ToAdd2 + TabAO1(2,iGrid,iCB)*TabAO2(1,iGrid,jCB)
            End Do
            AOInt(ij,1) = ToAdd1
            AOInt(ij,2) = ToAdd2
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
*
      Return
      End
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Subroutine Do_NInt2_d(AOInt,nAOInt,ndF_dRho, dF_dRho,
     &                      Weights,mGrid,Rho,nRho,
     &                      Scr,TabAO1,iCmp,iBas,nGrid_Tot,iSpin,
     &                      mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 AOInt(nAOInt*nAOInt,iSpin), dF_dRho(ndF_dRho,mGrid),
     &       Weights(mGrid), Rho(nRho,mGrid),
     &       TabAO1(mAO,mGrid,iBas*iCmp),
     &       Scr(iSpin*nFn,mGrid,iBas*iCmp)
*
      nBB = iBas*iBas
      nGrid_Tot=nGrid_Tot+mGrid*iBas*iCmp*iBas*iCmp
*
      If (iSpin.ne.1) Go To 99
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin = 1                                                        *
*                                                                      *
************************************************************************
*                                                                      *
      Do iCB = 1, iBas*iCmp
         iC = (iCB-1)/iBas + 1
         iB = iCB - (iC-1)*iBas
         iiC = (iC-1)*iCmp + iC
         iiB = (iB-1)*iBas + iB
         ii = (iiC-1)*nBB + iiB
*                                                                      *
************************************************************************
*                                                                      *
         ToAdd1=Zero
         Do iGrid = 1, mGrid

            gx=Rho(2,iGrid)*Weights(iGrid)
            gy=Rho(3,iGrid)*Weights(iGrid)
            gz=Rho(4,iGrid)*Weights(iGrid)

            Temp0=dF_dRho(ipR,iGrid)*Weights(iGrid)
            Temp1=gx*(2.0d0*dF_dRho(ipGxx,iGrid)+dF_dRho(ipGxy,iGrid))
            Temp2=gy*(2.0d0*dF_dRho(ipGxx,iGrid)+dF_dRho(ipGxy,iGrid))
            Temp3=gz*(2.0d0*dF_dRho(ipGxx,iGrid)+dF_dRho(ipGxy,iGrid))

            Scr(1,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp0
     &                       + TabAO1(2,iGrid,iCB) * Temp1
     &                       + TabAO1(3,iGrid,iCB) * Temp2
     &                       + TabAO1(4,iGrid,iCB) * Temp3
            Scr(2,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp1
            Scr(3,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp2
            Scr(4,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp3
            ToAdd1 = ToAdd1
     &             + Scr(1,iGrid,iCB) * TabAO1(1,iGrid,iCB)
     &             + Scr(2,iGrid,iCB) * TabAO1(2,iGrid,iCB)
     &             + Scr(3,iGrid,iCB) * TabAO1(3,iGrid,iCB)
     &             + Scr(4,iGrid,iCB) * TabAO1(4,iGrid,iCB)
         End Do
         AOInt(ii,1)=            ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         Do jCB = iCB+1, iBas*iCmp
            jC = (jCB-1)/iBas + 1
            jB = jCB - (jC-1)*iBas
            ijC = (jC-1)*iCmp + iC
            ijB = (jB-1)*iBas + iB
            ij = (ijC-1)*nBB + ijB
            jiC = (iC-1)*iCmp + jC
            jiB = (iB-1)*iBas + jB
            ji = (jiC-1)*nBB + jiB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + Scr(1,iGrid,iCB) * TabAO1(1,iGrid,jCB)
     &               + Scr(2,iGrid,iCB) * TabAO1(2,iGrid,jCB)
     &               + Scr(3,iGrid,iCB) * TabAO1(3,iGrid,jCB)
     &               + Scr(4,iGrid,iCB) * TabAO1(4,iGrid,jCB)
            End Do
            AOInt(ij,1)=            ToAdd1
            AOInt(ji,1)=            ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
      Return
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin =/= 1                                                      *
*                                                                      *
************************************************************************
*                                                                      *
 99   Continue
*                                                                      *
************************************************************************
*                                                                      *
      Do iCB = 1, iBas*iCmp
         iC = (iCB-1)/iBas + 1
         iB = iCB - (iC-1)*iBas
         iiC = (iC-1)*iCmp + iC
         iiB = (iB-1)*iBas + iB
         ii = (iiC-1)*nBB + iiB
*                                                                      *
************************************************************************
*                                                                      *
         ToAdd1=Zero
         ToAdd2=Zero
         Do iGrid = 1, mGrid

            gxa=rho(3,iGrid)*Weights(iGrid)
            gya=rho(4,iGrid)*Weights(iGrid)
            gza=rho(5,iGrid)*Weights(iGrid)
            gxb=rho(6,iGrid)*Weights(iGrid)
            gyb=rho(7,iGrid)*Weights(iGrid)
            gzb=rho(8,iGrid)*Weights(iGrid)

            Temp0a=dF_dRho(ipRa,iGrid) * Weights(iGrid)
            Temp0b=dF_dRho(ipRb,iGrid) * Weights(iGrid)
            Temp1a=2.0d0*dF_dRho(ipGaa,iGrid)*gxa
     &            +      dF_dRho(ipGab,iGrid)*gxb
            Temp1b=2.0d0*dF_dRho(ipGbb,iGrid)*gxb
     &            +      dF_dRho(ipGab,iGrid)*gxa
            Temp2a=2.0d0*dF_dRho(ipGaa,iGrid)*gya
     &            +      dF_dRho(ipGab,iGrid)*gyb
            Temp2b=2.0d0*dF_dRho(ipGbb,iGrid)*gyb
     &            +      dF_dRho(ipGab,iGrid)*gya
            Temp3a=2.0d0*dF_dRho(ipGaa,iGrid)*gza
     &            +      dF_dRho(ipGab,iGrid)*gzb
            Temp3b=2.0d0*dF_dRho(ipGbb,iGrid)*gzb
     &            +      dF_dRho(ipGab,iGrid)*gza
*
            Scr(1,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp0a
     &                       + TabAO1(2,iGrid,iCB) * Temp1a
     &                       + TabAO1(3,iGrid,iCB) * Temp2a
     &                       + TabAO1(4,iGrid,iCB) * Temp3a
            Scr(2,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp1a
            Scr(3,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp2a
            Scr(4,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp3a
            Scr(5,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp0b
     &                       + TabAO1(2,iGrid,iCB) * Temp1b
     &                       + TabAO1(3,iGrid,iCB) * Temp2b
     &                       + TabAO1(4,iGrid,iCB) * Temp3b
            Scr(6,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp1b
            Scr(7,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp2b
            Scr(8,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp3b
            ToAdd1= ToAdd1
     &            + Scr(1,iGrid,iCB) * TabAO1(1,iGrid,iCB)
     &            + Scr(2,iGrid,iCB) * TabAO1(2,iGrid,iCB)
     &            + Scr(3,iGrid,iCB) * TabAO1(3,iGrid,iCB)
     &            + Scr(4,iGrid,iCB) * TabAO1(4,iGrid,iCB)
            ToAdd2= ToAdd2
     &            + Scr(5,iGrid,iCB) * TabAO1(1,iGrid,iCB)
     &            + Scr(6,iGrid,iCB) * TabAO1(2,iGrid,iCB)
     &            + Scr(7,iGrid,iCB) * TabAO1(3,iGrid,iCB)
     &            + Scr(8,iGrid,iCB) * TabAO1(4,iGrid,iCB)
         End Do
         AOInt(ii,1)=            ToAdd1
         AOInt(ii,2)=            ToAdd2
*                                                                      *
************************************************************************
*                                                                      *
         Do jCB = iCB+1, iBas*iCmp
            jC = (jCB-1)/iBas + 1
            jB = jCB - (jC-1)*iBas
            ijC = (jC-1)*iCmp + iC
            ijB = (jB-1)*iBas + iB
            ij = (ijC-1)*nBB + ijB
            jiC = (iC-1)*iCmp + jC
            jiB = (iB-1)*iBas + jB
            ji = (jiC-1)*nBB + jiB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
            ToAdd2=Zero
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + Scr( 1,iGrid,iCB) * TabAO1(1,iGrid,jCB)
     &               + Scr( 2,iGrid,iCB) * TabAO1(2,iGrid,jCB)
     &               + Scr( 3,iGrid,iCB) * TabAO1(3,iGrid,jCB)
     &               + Scr( 4,iGrid,iCB) * TabAO1(4,iGrid,jCB)
               ToAdd2= ToAdd2
     &               + Scr( 5,iGrid,iCB) * TabAO1(1,iGrid,jCB)
     &               + Scr( 6,iGrid,iCB) * TabAO1(2,iGrid,jCB)
     &               + Scr( 7,iGrid,iCB) * TabAO1(3,iGrid,jCB)
     &               + Scr( 8,iGrid,iCB) * TabAO1(4,iGrid,jCB)
            End Do
            AOInt(ij,1)=            ToAdd1
            AOInt(ji,1)=            ToAdd1
            AOInt(ij,2)=            ToAdd2
            AOInt(ji,2)=            ToAdd2
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
*
      Return
      End
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Subroutine Do_NInt2(AOInt,nAOInt,mGrid,
     &                    TabAO1,iCmp,iBas,
     &                    TabAO2,jCmp,jBas,nGrid_Tot,iSpin,mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 AOInt(nAOInt*nAOInt,iSpin),
     &       TabAO1(nFn*iSpin,mGrid,iBas*iCmp),
     &       TabAO2(mAO,mGrid,jBas*jCmp)
*
      nBB = iBas*jBas
      nGrid_Tot=nGrid_Tot+mGrid*iBas*iCmp*jBas*jCmp
*
      If (iSpin.ne.1) Go To 99
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin = 1                                                        *
*                                                                      *
************************************************************************
*                                                                      *
      Do iCB = 1, iBas*iCmp
         iC = (iCB-1)/iBas + 1
         iB = iCB - (iC-1)*iBas
*                                                                      *
************************************************************************
*                                                                      *
*
         Do jCB = 1, jBas*jCmp
            jC = (jCB-1)/jBas + 1
            jB = jCB - (jC-1)*jBas
            ijC = (jC-1)*iCmp + iC
            ijB = (jB-1)*iBas + iB
            ij = (ijC-1)*nBB + ijB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
*
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + TabAO1(1,iGrid,iCB) * TabAO2(1,iGrid,jCB)
     &               + TabAO1(2,iGrid,iCB) * TabAO2(2,iGrid,jCB)
     &               + TabAO1(3,iGrid,iCB) * TabAO2(3,iGrid,jCB)
     &               + TabAO1(4,iGrid,iCB) * TabAO2(4,iGrid,jCB)
            End Do
            AOInt(ij,1)=            ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
*
      Return
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin =/= 1                                                      *
*                                                                      *
************************************************************************
*                                                                      *
 99   Continue
*                                                                      *
************************************************************************
*                                                                      *
      Do iCB = 1, iBas*iCmp
         iC = (iCB-1)/iBas + 1
         iB = iCB - (iC-1)*iBas
*                                                                      *
************************************************************************
*                                                                      *
         Do jCB = 1, jBas*jCmp
            jC = (jCB-1)/jBas + 1
            jB = jCB - (jC-1)*jBas
            ijC = (jC-1)*iCmp + iC
            ijB = (jB-1)*iBas + iB
            ij = (ijC-1)*nBB + ijB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
            ToAdd2=Zero
*
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + TabAO1( 1,iGrid,iCB) * TabAO2(1,iGrid,jCB)
     &               + TabAO1( 2,iGrid,iCB) * TabAO2(2,iGrid,jCB)
     &               + TabAO1( 3,iGrid,iCB) * TabAO2(3,iGrid,jCB)
     &               + TabAO1( 4,iGrid,iCB) * TabAO2(4,iGrid,jCB)
               ToAdd2= ToAdd2
     &               + TabAO1( 5,iGrid,iCB) * TabAO2(1,iGrid,jCB)
     &               + TabAO1( 6,iGrid,iCB) * TabAO2(2,iGrid,jCB)
     &               + TabAO1( 7,iGrid,iCB) * TabAO2(3,iGrid,jCB)
     &               + TabAO1( 8,iGrid,iCB) * TabAO2(4,iGrid,jCB)
            End Do
            AOInt(ij,1)=            ToAdd1
            AOInt(ij,2)=            ToAdd2
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Subroutine Do_NInt3_d(AOInt,nAOInt,ndF_dRho,
     &                      dF_dRho,Weights,mGrid,Rho,nRho,
     &                      Scr,TabAO1,iCmp,iBas,nGrid_Tot,iSpin,
     &                      mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 AOInt(nAOInt*nAOInt,iSpin), dF_dRho(ndF_dRho,mGrid),
     &       Weights(mGrid), Rho(nRho,mGrid),
     &       TabAO1(mAO,mGrid,iBas*iCmp),
     &       Scr(iSpin*nFn,mGrid,iBas*iCmp)
*
      nBB = iBas*iBas
      nGrid_Tot=nGrid_Tot+mGrid*iBas*iCmp*iBas*iCmp
*
      If (iSpin.ne.1) Go To 99
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin = 1                                                        *
*                                                                      *
************************************************************************
*                                                                      *
      Do iCB = 1, iBas*iCmp
         iC = (iCB-1)/iBas + 1
         iB = iCB - (iC-1)*iBas
         iiC = (iC-1)*iCmp + iC
         iiB = (iB-1)*iBas + iB
         ii = (iiC-1)*nBB + iiB
*                                                                      *
************************************************************************
*                                                                      *
         ToAdd1=Zero
         Do iGrid = 1, mGrid

            gx=Rho(2,iGrid)*Weights(iGrid)
            gy=Rho(3,iGrid)*Weights(iGrid)
            gz=Rho(4,iGrid)*Weights(iGrid)

            Temp0=dF_dRho(ipR,iGrid)*Weights(iGrid)
            Temp1=gx*(2.0d0*dF_dRho(ipGxx,iGrid)+dF_dRho(ipGxy,iGrid))
            Temp2=gy*(2.0d0*dF_dRho(ipGxx,iGrid)+dF_dRho(ipGxy,iGrid))
            Temp3=gz*(2.0d0*dF_dRho(ipGxx,iGrid)+dF_dRho(ipGxy,iGrid))

            Temp4=dF_dRho(ipT,iGrid)*Weights(iGrid)
            Temp5=dF_dRho(ipL,iGrid)*Weights(iGrid)

            Scr(1,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp0
     &                       + TabAO1(2,iGrid,iCB) * Temp1
     &                       + TabAO1(3,iGrid,iCB) * Temp2
     &                       + TabAO1(4,iGrid,iCB) * Temp3
            Scr(2,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp1
     &                       + TabAO1(2,iGrid,iCB) * Temp4
     &                       + TabAO1(2,iGrid,iCB) * Temp5
            Scr(3,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp2
     &                       + TabAO1(3,iGrid,iCB) * Temp4
     &                       + TabAO1(3,iGrid,iCB) * Temp5
            Scr(4,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp3
     &                       + TabAO1(4,iGrid,iCB) * Temp5
     &                       + TabAO1(4,iGrid,iCB) * Temp4
            Scr(5,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp5
            ToAdd1 = ToAdd1
     &             + Scr(1,iGrid,iCB) * TabAO1(1,iGrid,iCB)
     &             + Scr(2,iGrid,iCB) * TabAO1(2,iGrid,iCB)
     &             + Scr(3,iGrid,iCB) * TabAO1(3,iGrid,iCB)
     &             + Scr(4,iGrid,iCB) * TabAO1(4,iGrid,iCB)
     &             + Scr(5,iGrid,iCB) *(TabAO1(5,iGrid,iCB)
     &                                 +TabAO1(8,iGrid,iCB)
     &                                 +TabAO1(10,iGrid,iCB))
         End Do
         AOInt(ii,1)=            ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         Do jCB = iCB+1, iBas*iCmp
            jC = (jCB-1)/iBas + 1
            jB = jCB - (jC-1)*iBas
            ijC = (jC-1)*iCmp + iC
            ijB = (jB-1)*iBas + iB
            ij = (ijC-1)*nBB + ijB
            jiC = (iC-1)*iCmp + jC
            jiB = (iB-1)*iBas + jB
            ji = (jiC-1)*nBB + jiB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + Scr(1,iGrid,iCB) * TabAO1(1,iGrid,jCB)
     &               + Scr(2,iGrid,iCB) * TabAO1(2,iGrid,jCB)
     &               + Scr(3,iGrid,iCB) * TabAO1(3,iGrid,jCB)
     &               + Scr(4,iGrid,iCB) * TabAO1(4,iGrid,jCB)
     &               + Scr(5,iGrid,iCB) *(TabAO1(5,iGrid,jCB)
     &                                   +TabAO1(8,iGrid,jCB)
     &                                   +TabAO1(10,iGrid,jCB))
            End Do
            AOInt(ij,1)=            ToAdd1
            AOInt(ji,1)=            ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
*
      Return
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin =/= 1                                                      *
*                                                                      *
************************************************************************
*                                                                      *
 99   Continue
*                                                                      *
************************************************************************
*                                                                      *
      Do iCB = 1, iBas*iCmp
         iC = (iCB-1)/iBas + 1
         iB = iCB - (iC-1)*iBas
         iiC = (iC-1)*iCmp + iC
         iiB = (iB-1)*iBas + iB
         ii = (iiC-1)*nBB + iiB
*                                                                      *
************************************************************************
*                                                                      *
         ToAdd1=Zero
         ToAdd2=Zero
         Do iGrid = 1, mGrid

            gxa=rho(3,iGrid)*Weights(iGrid)
            gya=rho(4,iGrid)*Weights(iGrid)
            gza=rho(5,iGrid)*Weights(iGrid)
            gxb=rho(6,iGrid)*Weights(iGrid)
            gyb=rho(7,iGrid)*Weights(iGrid)
            gzb=rho(8,iGrid)*Weights(iGrid)

            Temp0a=dF_dRho(ipRa,iGrid) * Weights(iGrid)
            Temp0b=dF_dRho(ipRb,iGrid) * Weights(iGrid)
            Temp1a=2.0d0*dF_dRho(ipGaa,iGrid)*gxa
     &            +      dF_dRho(ipGab,iGrid)*gxb
            Temp1b=2.0d0*dF_dRho(ipGbb,iGrid)*gxb
     &            +      dF_dRho(ipGab,iGrid)*gxa
            Temp2a=2.0d0*dF_dRho(ipGaa,iGrid)*gya
     &            +      dF_dRho(ipGab,iGrid)*gyb
            Temp2b=2.0d0*dF_dRho(ipGbb,iGrid)*gyb
     &            +      dF_dRho(ipGab,iGrid)*gya
            Temp3a=2.0d0*dF_dRho(ipGaa,iGrid)*gza
     &            +      dF_dRho(ipGab,iGrid)*gzb
            Temp3b=2.0d0*dF_dRho(ipGbb,iGrid)*gzb
     &            +      dF_dRho(ipGab,iGrid)*gza
            Temp4a= dF_dRho(ipTa,iGrid)*Weights(iGrid)
            Temp4b= dF_dRho(ipTb,iGrid)*Weights(iGrid)
            Temp5a= dF_dRho(ipLa,iGrid)*Weights(iGrid)
            Temp5b= dF_dRho(ipLb,iGrid)*Weights(iGrid)
*
            Scr(1,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp0a
     &                       + TabAO1(2,iGrid,iCB) * Temp1a
     &                       + TabAO1(3,iGrid,iCB) * Temp2a
     &                       + TabAO1(4,iGrid,iCB) * Temp3a
            Scr(2,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp1a
     &                       + TabAO1(2,iGrid,iCB) * Temp4a
            Scr(3,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp2a
     &                       + TabAO1(3,iGrid,iCB) * Temp4a
            Scr(4,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp3a
     &                       + TabAO1(4,iGrid,iCB) * Temp4a
*
            Scr(5,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp0b
     &                       + TabAO1(2,iGrid,iCB) * Temp1b
     &                       + TabAO1(3,iGrid,iCB) * Temp2b
     &                       + TabAO1(4,iGrid,iCB) * Temp3b
            Scr( 6,iGrid,iCB)= TabAO1(1,iGrid,iCB) * Temp1b
     &                       + TabAO1(2,iGrid,iCB) * Temp4b
            Scr( 7,iGrid,iCB)= TabAO1(1,iGrid,iCB) * Temp2b
     &                       + TabAO1(3,iGrid,iCB) * Temp4b
            Scr( 8,iGrid,iCB)= TabAO1(1,iGrid,iCB) * Temp3b
     &                       + TabAO1(4,iGrid,iCB) * Temp4b
*
            Scr( 9,iGrid,iCB)= TabAO1(1,iGrid,iCB) * Temp5a
            Scr(10,iGrid,iCB)= TabAO1(1,iGrid,iCB) * Temp5b
            ToAdd1= ToAdd1
     &            + Scr( 1,iGrid,iCB) * TabAO1( 1,iGrid,iCB)
     &            + Scr( 2,iGrid,iCB) * TabAO1( 2,iGrid,iCB)
     &            + Scr( 3,iGrid,iCB) * TabAO1( 3,iGrid,iCB)
     &            + Scr( 4,iGrid,iCB) * TabAO1( 4,iGrid,iCB)
     &            + Scr( 9,iGrid,iCB) *(TabAO1( 5,iGrid,iCB)
     &                                + TabAO1( 8,iGrid,iCB)
     &                                + TabAO1(10,iGrid,iCB))
            ToAdd2= ToAdd2
     &            + Scr( 5,iGrid,iCB) * TabAO1( 1,iGrid,iCB)
     &            + Scr( 6,iGrid,iCB) * TabAO1( 2,iGrid,iCB)
     &            + Scr( 7,iGrid,iCB) * TabAO1( 3,iGrid,iCB)
     &            + Scr( 8,iGrid,iCB) * TabAO1( 4,iGrid,iCB)
     &            + Scr(10,iGrid,iCB) *(TabAO1( 5,iGrid,iCB)
     &                                + TabAO1( 8,iGrid,iCB)
     &                                + TabAO1(10,iGrid,iCB))
         End Do
         AOInt(ii,1)=            ToAdd1
         AOInt(ii,2)=            ToAdd2
*                                                                      *
************************************************************************
*                                                                      *
         Do jCB = iCB+1, iBas*iCmp
            jC = (jCB-1)/iBas + 1
            jB = jCB - (jC-1)*iBas
            ijC = (jC-1)*iCmp + iC
            ijB = (jB-1)*iBas + iB
            ij = (ijC-1)*nBB + ijB
            jiC = (iC-1)*iCmp + jC
            jiB = (iB-1)*iBas + jB
            ji = (jiC-1)*nBB + jiB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
            ToAdd2=Zero
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + Scr( 1,iGrid,iCB) * TabAO1(1,iGrid,jCB)
     &               + Scr( 2,iGrid,iCB) * TabAO1(2,iGrid,jCB)
     &               + Scr( 3,iGrid,iCB) * TabAO1(3,iGrid,jCB)
     &               + Scr( 4,iGrid,iCB) * TabAO1(4,iGrid,jCB)
     &               + Scr( 9,iGrid,iCB) *(TabAO1(5,iGrid,jCB)
     &                                   + TabAO1(8,iGrid,jCB)
     &                                   + TabAO1(10,iGrid,jCB))
               ToAdd2= ToAdd2
     &               + Scr( 5,iGrid,iCB) * TabAO1(1,iGrid,jCB)
     &               + Scr( 6,iGrid,iCB) * TabAO1(2,iGrid,jCB)
     &               + Scr( 7,iGrid,iCB) * TabAO1(3,iGrid,jCB)
     &               + Scr( 8,iGrid,iCB) * TabAO1(4,iGrid,jCB)
     &               + Scr(10,iGrid,iCB)*( TabAO1(5,iGrid,jCB)
     &                                   + TabAO1(8,iGrid,jCB)
     &                                   + TabAO1(10,iGrid,jCB))
            End Do
            AOInt(ij,1)=            ToAdd1
            AOInt(ji,1)=            ToAdd1
            AOInt(ij,2)=            ToAdd2
            AOInt(ji,2)=            ToAdd2
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
*
      Return
      End
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Subroutine Do_NInt3(AOInt,nAOInt,mGrid,
     &                    TabAO1,iCmp,iBas,
     &                    TabAO2,jCmp,jBas,nGrid_Tot,iSpin,mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 AOInt(nAOInt*nAOInt,iSpin),
     &       TabAO1(nFn*iSpin,mGrid,iBas*iCmp),
     &       TabAO2(mAO,mGrid,jBas*jCmp)
*
      nBB = iBas*jBas
      nGrid_Tot=nGrid_Tot+mGrid*iBas*iCmp*jBas*jCmp
*
      If (iSpin.ne.1) Go To 99
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin = 1                                                        *
*                                                                      *
************************************************************************
*                                                                      *
      Do iCB = 1, iBas*iCmp
         iC = (iCB-1)/iBas + 1
         iB = iCB - (iC-1)*iBas
*                                                                      *
************************************************************************
*                                                                      *
*
         Do jCB = 1, jBas*jCmp
            jC = (jCB-1)/jBas + 1
            jB = jCB - (jC-1)*jBas
            ijC = (jC-1)*iCmp + iC
            ijB = (jB-1)*iBas + iB
            ij = (ijC-1)*nBB + ijB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
*
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + TabAO1(1,iGrid,iCB) * TabAO2(1,iGrid,jCB)
     &               + TabAO1(2,iGrid,iCB) * TabAO2(2,iGrid,jCB)
     &               + TabAO1(3,iGrid,iCB) * TabAO2(3,iGrid,jCB)
     &               + TabAO1(4,iGrid,iCB) * TabAO2(4,iGrid,jCB)
     &               + TabAO1(5,iGrid,iCB) *(TabAO2(5,iGrid,jCB)
     &                                      +TabAO2(8,iGrid,jCB)
     &                                      +TabAO2(10,iGrid,jCB))
            End Do
            AOInt(ij,1)=            ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
*
      Return
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin =/= 1                                                      *
*                                                                      *
************************************************************************
*                                                                      *
 99   Continue
*                                                                      *
************************************************************************
*                                                                      *
      Do iCB = 1, iBas*iCmp
         iC = (iCB-1)/iBas + 1
         iB = iCB - (iC-1)*iBas
*                                                                      *
************************************************************************
*                                                                      *
         Do jCB = 1, jBas*jCmp
            jC = (jCB-1)/jBas + 1
            jB = jCB - (jC-1)*jBas
            ijC = (jC-1)*iCmp + iC
            ijB = (jB-1)*iBas + iB
            ij = (ijC-1)*nBB + ijB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
            ToAdd2=Zero
*
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + TabAO1( 1,iGrid,iCB) * TabAO2(1,iGrid,jCB)
     &               + TabAO1( 2,iGrid,iCB) * TabAO2(2,iGrid,jCB)
     &               + TabAO1( 3,iGrid,iCB) * TabAO2(3,iGrid,jCB)
     &               + TabAO1( 4,iGrid,iCB) * TabAO2(4,iGrid,jCB)
     &               + TabAO1( 9,iGrid,iCB) *(TabAO2(5,iGrid,jCB)
     &                                       +TabAO2(8,iGrid,jCB)
     &                                       +TabAO2(10,iGrid,jCB))
               ToAdd2= ToAdd2
     &               + TabAO1( 5,iGrid,iCB) * TabAO2(1,iGrid,jCB)
     &               + TabAO1( 6,iGrid,iCB) * TabAO2(2,iGrid,jCB)
     &               + TabAO1( 7,iGrid,iCB) * TabAO2(3,iGrid,jCB)
     &               + TabAO1( 8,iGrid,iCB) * TabAO2(4,iGrid,jCB)
     &               + TabAO1(10,iGrid,iCB) *(TabAO2(5,iGrid,jCB)
     &                                       +TabAO2(8,iGrid,jCB)
     &                                       +TabAO2(10,iGrid,jCB))
            End Do
            AOInt(ij,1)=            ToAdd1
            AOInt(ij,2)=            ToAdd2
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Subroutine Do_NInt4_d(AOInt,nAOInt,ndF_dRho,
     &                      dF_dRho,Weights,mGrid,Rho,nRho,
     &                      Scr,TabAO1,iCmp,iBas,nGrid_Tot,iSpin,
     &                      mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 AOInt(nAOInt*nAOInt,iSpin), dF_dRho(ndF_dRho,mGrid),
     &       Weights(mGrid), Rho(nRho,mGrid),
     &       TabAO1(mAO,mGrid,iBas*iCmp),
     &       Scr(iSpin*nFn,mGrid,iBas*iCmp)
*
      nBB = iBas*iBas
      nGrid_Tot=nGrid_Tot+mGrid*iBas*iCmp*iBas*iCmp
*
      If (iSpin.ne.1) Go To 99
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin = 1                                                        *
*                                                                      *
************************************************************************
*                                                                      *
      Do iCB = 1, iBas*iCmp
         iC = (iCB-1)/iBas + 1
         iB = iCB - (iC-1)*iBas
         iiC = (iC-1)*iCmp + iC
         iiB = (iB-1)*iBas + iB
         ii = (iiC-1)*nBB + iiB
*                                                                      *
************************************************************************
*                                                                      *
         ToAdd1=Zero
         Do iGrid = 1, mGrid
            gx=Rho(2,iGrid)*Weights(iGrid)
            gy=Rho(3,iGrid)*Weights(iGrid)
            gz=Rho(4,iGrid)*Weights(iGrid)

            Temp0=dF_dRho(ipR,iGrid)*Weights(iGrid)
            Temp1=gx*(2.0d0*dF_dRho(ipGxx,iGrid)+dF_dRho(ipGxy,iGrid))
            Temp2=gy*(2.0d0*dF_dRho(ipGxx,iGrid)+dF_dRho(ipGxy,iGrid))
            Temp3=gz*(2.0d0*dF_dRho(ipGxx,iGrid)+dF_dRho(ipGxy,iGrid))

            Temp4=dF_dRho(ipT,iGrid)*Weights(iGrid)

            Scr(1,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp0
     &                       + TabAO1(2,iGrid,iCB) * Temp1
     &                       + TabAO1(3,iGrid,iCB) * Temp2
     &                       + TabAO1(4,iGrid,iCB) * Temp3
            Scr(2,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp1
     &                       + TabAO1(2,iGrid,iCB) * Temp4
            Scr(3,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp2
     &                       + TabAO1(3,iGrid,iCB) * Temp4
            Scr(4,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp3
     &                       + TabAO1(4,iGrid,iCB) * Temp4
            ToAdd1 = ToAdd1
     &             + Scr(1,iGrid,iCB) * TabAO1(1,iGrid,iCB)
     &             + Scr(2,iGrid,iCB) * TabAO1(2,iGrid,iCB)
     &             + Scr(3,iGrid,iCB) * TabAO1(3,iGrid,iCB)
     &             + Scr(4,iGrid,iCB) * TabAO1(4,iGrid,iCB)
         End Do
         AOInt(ii,1)=            ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         Do jCB = iCB+1, iBas*iCmp
            jC = (jCB-1)/iBas + 1
            jB = jCB - (jC-1)*iBas
            ijC = (jC-1)*iCmp + iC
            ijB = (jB-1)*iBas + iB
            ij = (ijC-1)*nBB + ijB
            jiC = (iC-1)*iCmp + jC
            jiB = (iB-1)*iBas + jB
            ji = (jiC-1)*nBB + jiB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + Scr(1,iGrid,iCB) * TabAO1(1,iGrid,jCB)
     &               + Scr(2,iGrid,iCB) * TabAO1(2,iGrid,jCB)
     &               + Scr(3,iGrid,iCB) * TabAO1(3,iGrid,jCB)
     &               + Scr(4,iGrid,iCB) * TabAO1(4,iGrid,jCB)
            End Do
            AOInt(ij,1)=            ToAdd1
            AOInt(ji,1)=            ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
*
      Return
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin =/= 1                                                      *
*                                                                      *
************************************************************************
*                                                                      *
 99   Continue
*                                                                      *
************************************************************************
*                                                                      *
      Do iCB = 1, iBas*iCmp
         iC = (iCB-1)/iBas + 1
         iB = iCB - (iC-1)*iBas
         iiC = (iC-1)*iCmp + iC
         iiB = (iB-1)*iBas + iB
         ii = (iiC-1)*nBB + iiB
*                                                                      *
************************************************************************
*                                                                      *
         ToAdd1=Zero
         ToAdd2=Zero
         Do iGrid = 1, mGrid

            gxa=rho(3,iGrid)*Weights(iGrid)
            gya=rho(4,iGrid)*Weights(iGrid)
            gza=rho(5,iGrid)*Weights(iGrid)
            gxb=rho(6,iGrid)*Weights(iGrid)
            gyb=rho(7,iGrid)*Weights(iGrid)
            gzb=rho(8,iGrid)*Weights(iGrid)

            Temp0a=dF_dRho(ipRa,iGrid) * Weights(iGrid)
            Temp0b=dF_dRho(ipRb,iGrid) * Weights(iGrid)
            Temp1a=2.0d0*dF_dRho(ipGaa,iGrid)*gxa
     &            +      dF_dRho(ipGab,iGrid)*gxb
            Temp1b=2.0d0*dF_dRho(ipGbb,iGrid)*gxb
     &            +      dF_dRho(ipGab,iGrid)*gxa
            Temp2a=2.0d0*dF_dRho(ipGaa,iGrid)*gya
     &            +      dF_dRho(ipGab,iGrid)*gyb
            Temp2b=2.0d0*dF_dRho(ipGbb,iGrid)*gyb
     &            +      dF_dRho(ipGab,iGrid)*gya
            Temp3a=2.0d0*dF_dRho(ipGaa,iGrid)*gza
     &            +      dF_dRho(ipGab,iGrid)*gzb
            Temp3b=2.0d0*dF_dRho(ipGbb,iGrid)*gzb
     &            +      dF_dRho(ipGab,iGrid)*gza
            Temp4a= dF_dRho(ipTa,iGrid)*Weights(iGrid)
            Temp4b= dF_dRho(ipTb,iGrid)*Weights(iGrid)
*
            Scr(1,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp0a
     &                       + TabAO1(2,iGrid,iCB) * Temp1a
     &                       + TabAO1(3,iGrid,iCB) * Temp2a
     &                       + TabAO1(4,iGrid,iCB) * Temp3a
            Scr(2,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp1a
     &                       + TabAO1(2,iGrid,iCB) * Temp4a
            Scr(3,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp2a
     &                       + TabAO1(3,iGrid,iCB) * Temp4a
            Scr(4,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp3a
     &                       + TabAO1(4,iGrid,iCB) * Temp4a
*
            Scr(5,iGrid,iCB) = TabAO1(1,iGrid,iCB) * Temp0b
     &                       + TabAO1(2,iGrid,iCB) * Temp1b
     &                       + TabAO1(3,iGrid,iCB) * Temp2b
     &                       + TabAO1(4,iGrid,iCB) * Temp3b
            Scr( 6,iGrid,iCB)= TabAO1(1,iGrid,iCB) * Temp1b
     &                       + TabAO1(2,iGrid,iCB) * Temp4b
            Scr( 7,iGrid,iCB)= TabAO1(1,iGrid,iCB) * Temp2b
     &                       + TabAO1(3,iGrid,iCB) * Temp4b
            Scr( 8,iGrid,iCB)= TabAO1(1,iGrid,iCB) * Temp3b
     &                       + TabAO1(4,iGrid,iCB) * Temp4b
*
            ToAdd1= ToAdd1
     &            + Scr( 1,iGrid,iCB) * TabAO1( 1,iGrid,iCB)
     &            + Scr( 2,iGrid,iCB) * TabAO1( 2,iGrid,iCB)
     &            + Scr( 3,iGrid,iCB) * TabAO1( 3,iGrid,iCB)
     &            + Scr( 4,iGrid,iCB) * TabAO1( 4,iGrid,iCB)
            ToAdd2= ToAdd2
     &            + Scr( 5,iGrid,iCB) * TabAO1( 1,iGrid,iCB)
     &            + Scr( 6,iGrid,iCB) * TabAO1( 2,iGrid,iCB)
     &            + Scr( 7,iGrid,iCB) * TabAO1( 3,iGrid,iCB)
     &            + Scr( 8,iGrid,iCB) * TabAO1( 4,iGrid,iCB)
         End Do
         AOInt(ii,1)=            ToAdd1
         AOInt(ii,2)=            ToAdd2
*                                                                      *
************************************************************************
*                                                                      *
         Do jCB = iCB+1, iBas*iCmp
            jC = (jCB-1)/iBas + 1
            jB = jCB - (jC-1)*iBas
            ijC = (jC-1)*iCmp + iC
            ijB = (jB-1)*iBas + iB
            ij = (ijC-1)*nBB + ijB
            jiC = (iC-1)*iCmp + jC
            jiB = (iB-1)*iBas + jB
            ji = (jiC-1)*nBB + jiB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
            ToAdd2=Zero
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + Scr( 1,iGrid,iCB) * TabAO1(1,iGrid,jCB)
     &               + Scr( 2,iGrid,iCB) * TabAO1(2,iGrid,jCB)
     &               + Scr( 3,iGrid,iCB) * TabAO1(3,iGrid,jCB)
     &               + Scr( 4,iGrid,iCB) * TabAO1(4,iGrid,jCB)
               ToAdd2= ToAdd2
     &               + Scr( 5,iGrid,iCB) * TabAO1(1,iGrid,jCB)
     &               + Scr( 6,iGrid,iCB) * TabAO1(2,iGrid,jCB)
     &               + Scr( 7,iGrid,iCB) * TabAO1(3,iGrid,jCB)
     &               + Scr( 8,iGrid,iCB) * TabAO1(4,iGrid,jCB)
            End Do
            AOInt(ij,1)=            ToAdd1
            AOInt(ji,1)=            ToAdd1
            AOInt(ij,2)=            ToAdd2
            AOInt(ji,2)=            ToAdd2
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
*
      Return
      End
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Subroutine Do_NInt4(AOInt,nAOInt,mGrid,
     &                    TabAO1,iCmp,iBas,
     &                    TabAO2,jCmp,jBas,nGrid_Tot,iSpin,mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 AOInt(nAOInt*nAOInt,iSpin),
     &       TabAO1(nFn*iSpin,mGrid,iBas*iCmp),
     &       TabAO2(mAO,mGrid,jBas*jCmp)
*
      nBB = iBas*jBas
      nGrid_Tot=nGrid_Tot+mGrid*iBas*iCmp*jBas*jCmp
*
      If (iSpin.ne.1) Go To 99
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin = 1                                                        *
*                                                                      *
************************************************************************
*                                                                      *
      Do iCB = 1, iBas*iCmp
         iC = (iCB-1)/iBas + 1
         iB = iCB - (iC-1)*iBas
*                                                                      *
************************************************************************
*                                                                      *
*
         Do jCB = 1, jBas*jCmp
            jC = (jCB-1)/jBas + 1
            jB = jCB - (jC-1)*jBas
            ijC = (jC-1)*iCmp + iC
            ijB = (jB-1)*iBas + iB
            ij = (ijC-1)*nBB + ijB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
*
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + TabAO1(1,iGrid,iCB) * TabAO2(1,iGrid,jCB)
     &               + TabAO1(2,iGrid,iCB) * TabAO2(2,iGrid,jCB)
     &               + TabAO1(3,iGrid,iCB) * TabAO2(3,iGrid,jCB)
     &               + TabAO1(4,iGrid,iCB) * TabAO2(4,iGrid,jCB)
            End Do
            AOInt(ij,1)=            ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
*
      Return
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin =/= 1                                                      *
*                                                                      *
************************************************************************
*                                                                      *
 99   Continue
*                                                                      *
************************************************************************
*                                                                      *
      Do iCB = 1, iBas*iCmp
         iC = (iCB-1)/iBas + 1
         iB = iCB - (iC-1)*iBas
*                                                                      *
************************************************************************
*                                                                      *
         Do jCB = 1, jBas*jCmp
            jC = (jCB-1)/jBas + 1
            jB = jCB - (jC-1)*jBas
            ijC = (jC-1)*iCmp + iC
            ijB = (jB-1)*iBas + iB
            ij = (ijC-1)*nBB + ijB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
            ToAdd2=Zero
*
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + TabAO1( 1,iGrid,iCB) * TabAO2(1,iGrid,jCB)
     &               + TabAO1( 2,iGrid,iCB) * TabAO2(2,iGrid,jCB)
     &               + TabAO1( 3,iGrid,iCB) * TabAO2(3,iGrid,jCB)
     &               + TabAO1( 4,iGrid,iCB) * TabAO2(4,iGrid,jCB)
               ToAdd2= ToAdd2
     &               + TabAO1( 5,iGrid,iCB) * TabAO2(1,iGrid,jCB)
     &               + TabAO1( 6,iGrid,iCB) * TabAO2(2,iGrid,jCB)
     &               + TabAO1( 7,iGrid,iCB) * TabAO2(3,iGrid,jCB)
     &               + TabAO1( 8,iGrid,iCB) * TabAO2(4,iGrid,jCB)
            End Do
            AOInt(ij,1)=            ToAdd1
            AOInt(ij,2)=            ToAdd2
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
