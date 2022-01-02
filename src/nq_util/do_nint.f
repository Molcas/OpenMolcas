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
      Subroutine Do_NInt1_d(ndF_dRho,dF_dRho,
     &                      Weights,mGrid,
     &                      Scr,TabAO1,nBfn,nGrid_Tot,iSpin,
     &                      mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 dF_dRho(ndF_dRho,mGrid), Weights(mGrid),
     &       TabAO1(mAO,mGrid,nBfn), Scr(nFn,mGrid,nBfn,iSpin)

      nGrid_Tot=nGrid_Tot+mGrid*nBfn**2
*
      If (iSpin.ne.1) Go To 99
*
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
         Do iGrid = 1, mGrid
            Tmp = TabAO1(1,iGrid,iCB)
     &                   * dF_dRho(ipR,iGrid) * Weights(iGrid)
            Scr(1,iGrid,iCB,1) = Tmp
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End Do
*
      Return
*
 99   Continue
*
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
         Do iGrid = 1, mGrid
            Tmp1=TabAO1(1,iGrid,iCB)*dF_dRho(ipRa,iGrid)*Weights(iGrid)
            Tmp2=TabAO1(1,iGrid,iCB)*dF_dRho(ipRb,iGrid)*Weights(iGrid)
            Scr(1,iGrid,iCB,1) = Tmp1
            Scr(1,iGrid,iCB,2) = Tmp2
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End Do
*
      Return
      End
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Subroutine Do_NInt2_d(ndF_dRho, dF_dRho,
     &                      Weights,mGrid,
     &                      Scr,TabAO1,nBfn,nGrid_Tot,iSpin,
     &                      mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      use nq_Grid, only: GradRho
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 dF_dRho(ndF_dRho,mGrid),
     &       Weights(mGrid),
     &       TabAO1(mAO,mGrid,nBfn),
     &       Scr(nFn,mGrid,nBfn,iSpin)
*
      nGrid_Tot=nGrid_Tot+mGrid*nBfn**2
*
      If (iSpin.ne.1) Go To 99
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

            Temp0=dF_dRho(ipR,iGrid)*Weights(iGrid)
            Temp1=gx*(2.0d0*dF_dRho(ipGxx,iGrid)+dF_dRho(ipGxy,iGrid))
            Temp2=gy*(2.0d0*dF_dRho(ipGxx,iGrid)+dF_dRho(ipGxy,iGrid))
            Temp3=gz*(2.0d0*dF_dRho(ipGxx,iGrid)+dF_dRho(ipGxy,iGrid))

            Scr(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp0
     &                         + TabAO1(2,iGrid,iCB) * Temp1
     &                         + TabAO1(3,iGrid,iCB) * Temp2
     &                         + TabAO1(4,iGrid,iCB) * Temp3
            Scr(2,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp1
            Scr(3,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp2
            Scr(4,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp3
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
 99   Continue
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
            Scr(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp0a
     &                         + TabAO1(2,iGrid,iCB) * Temp1a
     &                         + TabAO1(3,iGrid,iCB) * Temp2a
     &                         + TabAO1(4,iGrid,iCB) * Temp3a
            Scr(2,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp1a
            Scr(3,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp2a
            Scr(4,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp3a
            Scr(1,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp0b
     &                         + TabAO1(2,iGrid,iCB) * Temp1b
     &                         + TabAO1(3,iGrid,iCB) * Temp2b
     &                         + TabAO1(4,iGrid,iCB) * Temp3b
            Scr(2,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp1b
            Scr(3,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp2b
            Scr(4,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp3b
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End Do
*
      Return
      End
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Subroutine Do_NInt3_d(ndF_dRho,
     &                      dF_dRho,Weights,mGrid,
     &                      Scr,TabAO1,nBfn,nGrid_Tot,iSpin,
     &                      mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      use nq_Grid, only: GradRho
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 dF_dRho(ndF_dRho,mGrid),
     &       Weights(mGrid),
     &       TabAO1(mAO,mGrid,nBfn),
     &       Scr(nFn,mGrid,nBfn,iSpin)
*
      nGrid_Tot=nGrid_Tot+mGrid*nBfn**2
*
      If (iSpin.ne.1) Go To 99
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

            Temp0=dF_dRho(ipR,iGrid)*Weights(iGrid)
            Temp1=gx*(2.0d0*dF_dRho(ipGxx,iGrid)+dF_dRho(ipGxy,iGrid))
            Temp2=gy*(2.0d0*dF_dRho(ipGxx,iGrid)+dF_dRho(ipGxy,iGrid))
            Temp3=gz*(2.0d0*dF_dRho(ipGxx,iGrid)+dF_dRho(ipGxy,iGrid))

            Temp4=dF_dRho(ipT,iGrid)*Weights(iGrid)
            Temp5=dF_dRho(ipL,iGrid)*Weights(iGrid)

            Scr(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp0
     &                         + TabAO1(2,iGrid,iCB) * Temp1
     &                         + TabAO1(3,iGrid,iCB) * Temp2
     &                         + TabAO1(4,iGrid,iCB) * Temp3
            Scr(2,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp1
     &                         + TabAO1(2,iGrid,iCB) * Temp4
     &                         + TabAO1(2,iGrid,iCB) * Temp5
            Scr(3,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp2
     &                         + TabAO1(3,iGrid,iCB) * Temp4
     &                         + TabAO1(3,iGrid,iCB) * Temp5
            Scr(4,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp3
     &                         + TabAO1(4,iGrid,iCB) * Temp5
     &                         + TabAO1(4,iGrid,iCB) * Temp4
            Scr(5,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp5
         End Do
*                                                                      *
************************************************************************
*                                                                      *
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
            Scr(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp0a
     &                       + TabAO1(2,iGrid,iCB) * Temp1a
     &                       + TabAO1(3,iGrid,iCB) * Temp2a
     &                       + TabAO1(4,iGrid,iCB) * Temp3a
            Scr(2,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp1a
     &                       + TabAO1(2,iGrid,iCB) * Temp4a
            Scr(3,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp2a
     &                       + TabAO1(3,iGrid,iCB) * Temp4a
            Scr(4,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp3a
     &                       + TabAO1(4,iGrid,iCB) * Temp4a
            Scr(5,iGrid,iCB,1)= TabAO1(1,iGrid,iCB) * Temp5a
*
            Scr(1,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp0b
     &                       + TabAO1(2,iGrid,iCB) * Temp1b
     &                       + TabAO1(3,iGrid,iCB) * Temp2b
     &                       + TabAO1(4,iGrid,iCB) * Temp3b
            Scr(2,iGrid,iCB,2)= TabAO1(1,iGrid,iCB) * Temp1b
     &                       + TabAO1(2,iGrid,iCB) * Temp4b
            Scr(3,iGrid,iCB,2)= TabAO1(1,iGrid,iCB) * Temp2b
     &                       + TabAO1(3,iGrid,iCB) * Temp4b
            Scr(4,iGrid,iCB,2)= TabAO1(1,iGrid,iCB) * Temp3b
     &                       + TabAO1(4,iGrid,iCB) * Temp4b
            Scr(5,iGrid,iCB,2)= TabAO1(1,iGrid,iCB) * Temp5b
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End Do
*
      Return
      End
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Subroutine Do_NInt4_d(ndF_dRho,
     &                      dF_dRho,Weights,mGrid,
     &                      Scr,TabAO1,nBfn,nGrid_Tot,iSpin,
     &                      mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      use nq_Grid, only: GradRho
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 dF_dRho(ndF_dRho,mGrid),
     &       Weights(mGrid),
     &       TabAO1(mAO,mGrid,nBfn),
     &       Scr(nFn,mGrid,nBfn,iSpin)
*
      nGrid_Tot=nGrid_Tot+mGrid*nBfn**2
*
      If (iSpin.ne.1) Go To 99
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

            Temp0=dF_dRho(ipR,iGrid)*Weights(iGrid)
            Temp1=gx*(2.0d0*dF_dRho(ipGxx,iGrid)+dF_dRho(ipGxy,iGrid))
            Temp2=gy*(2.0d0*dF_dRho(ipGxx,iGrid)+dF_dRho(ipGxy,iGrid))
            Temp3=gz*(2.0d0*dF_dRho(ipGxx,iGrid)+dF_dRho(ipGxy,iGrid))

            Temp4=dF_dRho(ipT,iGrid)*Weights(iGrid)

            Scr(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp0
     &                         + TabAO1(2,iGrid,iCB) * Temp1
     &                         + TabAO1(3,iGrid,iCB) * Temp2
     &                         + TabAO1(4,iGrid,iCB) * Temp3
            Scr(2,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp1
     &                         + TabAO1(2,iGrid,iCB) * Temp4
            Scr(3,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp2
     &                         + TabAO1(3,iGrid,iCB) * Temp4
            Scr(4,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp3
     &                         + TabAO1(4,iGrid,iCB) * Temp4
         End Do
*                                                                      *
************************************************************************
*                                                                      *
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
            Scr(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp0a
     &                         + TabAO1(2,iGrid,iCB) * Temp1a
     &                         + TabAO1(3,iGrid,iCB) * Temp2a
     &                         + TabAO1(4,iGrid,iCB) * Temp3a
            Scr(2,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp1a
     &                         + TabAO1(2,iGrid,iCB) * Temp4a
            Scr(3,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp2a
     &                         + TabAO1(3,iGrid,iCB) * Temp4a
            Scr(4,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp3a
     &                         + TabAO1(4,iGrid,iCB) * Temp4a
*
            Scr(1,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp0b
     &                         + TabAO1(2,iGrid,iCB) * Temp1b
     &                         + TabAO1(3,iGrid,iCB) * Temp2b
     &                         + TabAO1(4,iGrid,iCB) * Temp3b
            Scr(2,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp1b
     &                         + TabAO1(2,iGrid,iCB) * Temp4b
            Scr(3,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp2b
     &                         + TabAO1(3,iGrid,iCB) * Temp4b
            Scr(4,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp3b
     &                         + TabAO1(4,iGrid,iCB) * Temp4b
*
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End Do
*
      Return
      End
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
