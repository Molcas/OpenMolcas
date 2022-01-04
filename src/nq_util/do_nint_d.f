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
      Subroutine Do_NInt_d(ndF_dRho,dF_dRho,mGrid,
     &                      Grid_AO,TabAO1,nBfn,iSpin,mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      use nq_Grid, only: GradRho, Weights
      use nq_Grid, only: vRho, vTau, vLapl
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
#include "nq_info.fh"
      Real*8 dF_dRho(ndF_dRho,mGrid),
     &       TabAO1(mAO,mGrid,nBfn), Grid_AO(nFn,mGrid,nBfn,iSpin)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      If (Functional_type.eq.LDA_type) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      If (iSpin.ne.1) Go To 99
*
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
         Do iGrid = 1, mGrid
            Tmp = TabAO1(1,iGrid,iCB)
     &                   * vRho(1,iGrid) * Weights(iGrid)
            Grid_AO(1,iGrid,iCB,1) = Tmp
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End Do
      Return
*
 99   Continue
*
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
         Do iGrid = 1, mGrid
            Tmp1=TabAO1(1,iGrid,iCB)*vRho(1,iGrid)*Weights(iGrid)
            Tmp2=TabAO1(1,iGrid,iCB)*vRho(2,iGrid)*Weights(iGrid)
            Grid_AO(1,iGrid,iCB,1) = Tmp1
            Grid_AO(1,iGrid,iCB,2) = Tmp2
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End Do
      Return
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.GGA_type) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*
      If (iSpin.ne.1) Go To 98
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin = 1                                                        *
*                                                                      *
************************************************************************
*                                                                      *
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
         Do iGrid = 1, mGrid

            gx=GradRho(1,iGrid)*Weights(iGrid)
            gy=GradRho(2,iGrid)*Weights(iGrid)
            gz=GradRho(3,iGrid)*Weights(iGrid)

            Temp0=vRho(1,iGrid)*Weights(iGrid)
            Temp1=gx*2.0d0*dF_dRho(ipGxx,iGrid)
            Temp2=gy*2.0d0*dF_dRho(ipGxx,iGrid)
            Temp3=gz*2.0d0*dF_dRho(ipGxx,iGrid)

            Grid_AO(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp0
     &                         + TabAO1(2,iGrid,iCB) * Temp1
     &                         + TabAO1(3,iGrid,iCB) * Temp2
     &                         + TabAO1(4,iGrid,iCB) * Temp3
            Grid_AO(2,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp1
            Grid_AO(3,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp2
            Grid_AO(4,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp3
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End Do
      Return
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin =/= 1                                                      *
*                                                                      *
************************************************************************
*                                                                      *
 98   Continue
*                                                                      *
************************************************************************
*                                                                      *
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
         Do iGrid = 1, mGrid

            gxa=Gradrho(1,iGrid)*Weights(iGrid)
            gya=Gradrho(2,iGrid)*Weights(iGrid)
            gza=Gradrho(3,iGrid)*Weights(iGrid)
            gxb=Gradrho(4,iGrid)*Weights(iGrid)
            gyb=Gradrho(5,iGrid)*Weights(iGrid)
            gzb=Gradrho(6,iGrid)*Weights(iGrid)

            Temp0a=vRho(1,iGrid) * Weights(iGrid)
            Temp0b=vRho(2,iGrid) * Weights(iGrid)
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
            Grid_AO(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp0a
     &                         + TabAO1(2,iGrid,iCB) * Temp1a
     &                         + TabAO1(3,iGrid,iCB) * Temp2a
     &                         + TabAO1(4,iGrid,iCB) * Temp3a
            Grid_AO(2,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp1a
            Grid_AO(3,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp2a
            Grid_AO(4,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp3a
            Grid_AO(1,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp0b
     &                         + TabAO1(2,iGrid,iCB) * Temp1b
     &                         + TabAO1(3,iGrid,iCB) * Temp2b
     &                         + TabAO1(4,iGrid,iCB) * Temp3b
            Grid_AO(2,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp1b
            Grid_AO(3,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp2b
            Grid_AO(4,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp3b
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End Do
      Return
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.meta_GGA_type2) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      If (iSpin.ne.1) Go To 97
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin = 1                                                        *
*                                                                      *
************************************************************************
*                                                                      *
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
         Do iGrid = 1, mGrid

            gx=GradRho(1,iGrid)*Weights(iGrid)
            gy=GradRho(2,iGrid)*Weights(iGrid)
            gz=GradRho(3,iGrid)*Weights(iGrid)

            Temp0=vRho(1,iGrid)*Weights(iGrid)
            Temp1=gx*2.0d0*dF_dRho(ipGxx,iGrid)
            Temp2=gy*2.0d0*dF_dRho(ipGxx,iGrid)
            Temp3=gz*2.0d0*dF_dRho(ipGxx,iGrid)

            Temp4=vTau(1,iGrid)*Weights(iGrid)
            Temp5=vLapl(1,iGrid)*Weights(iGrid)

            Grid_AO(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp0
     &                         + TabAO1(2,iGrid,iCB) * Temp1
     &                         + TabAO1(3,iGrid,iCB) * Temp2
     &                         + TabAO1(4,iGrid,iCB) * Temp3
            Grid_AO(2,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp1
     &                         + TabAO1(2,iGrid,iCB) * Temp4
     &                         + TabAO1(2,iGrid,iCB) * Temp5
            Grid_AO(3,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp2
     &                         + TabAO1(3,iGrid,iCB) * Temp4
     &                         + TabAO1(3,iGrid,iCB) * Temp5
            Grid_AO(4,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp3
     &                         + TabAO1(4,iGrid,iCB) * Temp5
     &                         + TabAO1(4,iGrid,iCB) * Temp4
            Grid_AO(5,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp5
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End Do
      Return
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin =/= 1                                                      *
*                                                                      *
************************************************************************
*                                                                      *
 97   Continue
*                                                                      *
************************************************************************
*                                                                      *
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
         Do iGrid = 1, mGrid

            gxa=Gradrho(1,iGrid)*Weights(iGrid)
            gya=Gradrho(2,iGrid)*Weights(iGrid)
            gza=Gradrho(3,iGrid)*Weights(iGrid)
            gxb=Gradrho(4,iGrid)*Weights(iGrid)
            gyb=Gradrho(5,iGrid)*Weights(iGrid)
            gzb=Gradrho(6,iGrid)*Weights(iGrid)

            Temp0a=vRho(1,iGrid) * Weights(iGrid)
            Temp0b=vRho(2,iGrid) * Weights(iGrid)
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
            Temp4a= vTau(1,iGrid)*Weights(iGrid)
            Temp4b= vTau(2,iGrid)*Weights(iGrid)
            Temp5a= vLapl(1,iGrid)*Weights(iGrid)
            Temp5b= vLapl(2,iGrid)*Weights(iGrid)
*
            Grid_AO(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp0a
     &                       + TabAO1(2,iGrid,iCB) * Temp1a
     &                       + TabAO1(3,iGrid,iCB) * Temp2a
     &                       + TabAO1(4,iGrid,iCB) * Temp3a
            Grid_AO(2,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp1a
     &                       + TabAO1(2,iGrid,iCB) * Temp4a
            Grid_AO(3,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp2a
     &                       + TabAO1(3,iGrid,iCB) * Temp4a
            Grid_AO(4,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp3a
     &                       + TabAO1(4,iGrid,iCB) * Temp4a
            Grid_AO(5,iGrid,iCB,1)= TabAO1(1,iGrid,iCB) * Temp5a
*
            Grid_AO(1,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp0b
     &                       + TabAO1(2,iGrid,iCB) * Temp1b
     &                       + TabAO1(3,iGrid,iCB) * Temp2b
     &                       + TabAO1(4,iGrid,iCB) * Temp3b
            Grid_AO(2,iGrid,iCB,2)= TabAO1(1,iGrid,iCB) * Temp1b
     &                       + TabAO1(2,iGrid,iCB) * Temp4b
            Grid_AO(3,iGrid,iCB,2)= TabAO1(1,iGrid,iCB) * Temp2b
     &                       + TabAO1(3,iGrid,iCB) * Temp4b
            Grid_AO(4,iGrid,iCB,2)= TabAO1(1,iGrid,iCB) * Temp3b
     &                       + TabAO1(4,iGrid,iCB) * Temp4b
            Grid_AO(5,iGrid,iCB,2)= TabAO1(1,iGrid,iCB) * Temp5b
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End Do
      Return
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.meta_GGA_type1) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*
      If (iSpin.ne.1) Go To 96
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin = 1                                                        *
*                                                                      *
************************************************************************
*                                                                      *
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
         Do iGrid = 1, mGrid
            gx=GradRho(1,iGrid)*Weights(iGrid)
            gy=GradRho(2,iGrid)*Weights(iGrid)
            gz=GradRho(3,iGrid)*Weights(iGrid)

            Temp0=vRho(1,iGrid)*Weights(iGrid)
            Temp1=gx*2.0d0*dF_dRho(ipGxx,iGrid)
            Temp2=gy*2.0d0*dF_dRho(ipGxx,iGrid)
            Temp3=gz*2.0d0*dF_dRho(ipGxx,iGrid)

            Temp4=vTau(1,iGrid)*Weights(iGrid)

            Grid_AO(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp0
     &                         + TabAO1(2,iGrid,iCB) * Temp1
     &                         + TabAO1(3,iGrid,iCB) * Temp2
     &                         + TabAO1(4,iGrid,iCB) * Temp3
            Grid_AO(2,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp1
     &                         + TabAO1(2,iGrid,iCB) * Temp4
            Grid_AO(3,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp2
     &                         + TabAO1(3,iGrid,iCB) * Temp4
            Grid_AO(4,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp3
     &                         + TabAO1(4,iGrid,iCB) * Temp4
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End Do
      Return
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin =/= 1                                                      *
*                                                                      *
************************************************************************
*                                                                      *
 96   Continue
*                                                                      *
************************************************************************
*                                                                      *
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
         Do iGrid = 1, mGrid

            gxa=Gradrho(1,iGrid)*Weights(iGrid)
            gya=Gradrho(2,iGrid)*Weights(iGrid)
            gza=Gradrho(3,iGrid)*Weights(iGrid)
            gxb=Gradrho(4,iGrid)*Weights(iGrid)
            gyb=Gradrho(5,iGrid)*Weights(iGrid)
            gzb=Gradrho(6,iGrid)*Weights(iGrid)

            Temp0a=vRho(1,iGrid) * Weights(iGrid)
            Temp0b=vRho(2,iGrid) * Weights(iGrid)
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
            Temp4a= vTau(1,iGrid)*Weights(iGrid)
            Temp4b= vTau(2,iGrid)*Weights(iGrid)
*
            Grid_AO(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp0a
     &                         + TabAO1(2,iGrid,iCB) * Temp1a
     &                         + TabAO1(3,iGrid,iCB) * Temp2a
     &                         + TabAO1(4,iGrid,iCB) * Temp3a
            Grid_AO(2,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp1a
     &                         + TabAO1(2,iGrid,iCB) * Temp4a
            Grid_AO(3,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp2a
     &                         + TabAO1(3,iGrid,iCB) * Temp4a
            Grid_AO(4,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp3a
     &                         + TabAO1(4,iGrid,iCB) * Temp4a
*
            Grid_AO(1,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp0b
     &                         + TabAO1(2,iGrid,iCB) * Temp1b
     &                         + TabAO1(3,iGrid,iCB) * Temp2b
     &                         + TabAO1(4,iGrid,iCB) * Temp3b
            Grid_AO(2,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp1b
     &                         + TabAO1(2,iGrid,iCB) * Temp4b
            Grid_AO(3,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp2b
     &                         + TabAO1(3,iGrid,iCB) * Temp4b
            Grid_AO(4,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp3b
     &                         + TabAO1(4,iGrid,iCB) * Temp4b
*
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End Do
      Return
*                                                                      *
************************************************************************
************************************************************************
*
      Else
*                                                                      *
************************************************************************
************************************************************************
*
         Write (6,*) 'DFT_Int: Illegal functional type!'
         Call Abend()

*                                                                      *
************************************************************************
************************************************************************
*
      End If
*                                                                      *
************************************************************************
************************************************************************
      Return
      End
