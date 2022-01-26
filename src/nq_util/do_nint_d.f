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
      Subroutine Do_NInt_d(mGrid,Grid_AO,TabAO1,nBfn,nD,mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
!     use nq_Grid, only: Rho
      use nq_Grid, only: GradRho, Weights
      use nq_Grid, only: vRho, vSigma, vTau, vLapl
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_info.fh"
      Real*8 TabAO1(mAO,mGrid,nBfn), Grid_AO(nFn,mGrid,nBfn,nD)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
!     Thr=1.0D-14
!     If (nD.eq.1) Thr=Thr/Two
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Select Case(Functional_type)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Case (LDA_type)
!
!        F(Rho)
!
!        for the integrals we need:
!
!        phi_i dF/dRho phi_j
!
!        Grid_AO contains
!        1: phi_i dF/dRho
!
!        Final integral assembled as, done in do_nIntx.
!
!        Grid_AO(1)_i phi_j
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Select Case (nD)
*                                                                      *
************************************************************************
*                                                                      *
      Case(1)
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, mGrid

!        If (Rho(1,iGrid)<Thr) Cycle
         Tmp =  vRho(1,iGrid) * Weights(iGrid)


         Do iCB = 1, nBfn
            Grid_AO(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Tmp
         End Do

      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Case(2)
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, mGrid

!        If (Rho(1,iGrid)+Rho(2,iGrid)<Thr) Cycle
         Tmp1=  vRho(1,iGrid) * Weights(iGrid)
         Tmp2=  vRho(2,iGrid) * Weights(iGrid)

         Do iCB = 1, nBfn
            Grid_AO(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Tmp1
            Grid_AO(1,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Tmp2
         End Do

      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Case Default
        Write (6,*) 'Invalid nD value:', nD
        Call Abend()
      End Select  ! nD
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Case (GGA_type)
!
!        F(Rho,Sigma) : Sigma=GradRho*GradRho
!
!        dF/dGradRho = dF/dSigma dSigma/dGradRho = 2 dF/dSigma GradRho
!
!        for the integrals we need:
!
!        phi_i dF/dRho phi_j
!     +  phi_i 2 (dF/dSigma) {GradRho Grad(phi_j)}
!     +  {Grad(phi_i) GradRho} 2 (dF/dSigma) phi_j
!
!        Grid_AO contains
!        1:  0.5 * phi_i dF/dRho
!          + {Grad(phi_i GradRho} 2 (dF/dSigma)
!
!        Final integral assembled as, done in do_nIntx.
!
!        Grid_AO(1)_i phi_j
!      + phi_i Grid_AO(1)_j
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*
      Select Case(nD)
*                                                                      *
************************************************************************
*                                                                      *
      Case(1)
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, mGrid

!        If (Rho(1,iGrid)<Thr) Cycle
         gx=GradRho(1,iGrid)*Weights(iGrid)
         gy=GradRho(2,iGrid)*Weights(iGrid)
         gz=GradRho(3,iGrid)*Weights(iGrid)

         Temp0=0.5D0*vRho(1,iGrid)*Weights(iGrid)
         Temp1=gx*2.0d0*vSigma(1,iGrid)
         Temp2=gy*2.0d0*vSigma(1,iGrid)
         Temp3=gz*2.0d0*vSigma(1,iGrid)

         Do iCB = 1, nBfn

            Grid_AO(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp0
     &                             + TabAO1(2,iGrid,iCB) * Temp1
     &                             + TabAO1(3,iGrid,iCB) * Temp2
     &                             + TabAO1(4,iGrid,iCB) * Temp3
         End Do

      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Case(2)
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, mGrid

!        If (Rho(1,iGrid)+Rho(2,iGrid)<Thr) Cycle
         gxa=Gradrho(1,iGrid)*Weights(iGrid)
         gya=Gradrho(2,iGrid)*Weights(iGrid)
         gza=Gradrho(3,iGrid)*Weights(iGrid)
         gxb=Gradrho(4,iGrid)*Weights(iGrid)
         gyb=Gradrho(5,iGrid)*Weights(iGrid)
         gzb=Gradrho(6,iGrid)*Weights(iGrid)

         Temp0a=0.5D0*vRho(1,iGrid) * Weights(iGrid)
         Temp0b=0.5D0*vRho(2,iGrid) * Weights(iGrid)
         Temp1a=2.0d0*vSigma(1,iGrid)*gxa + vSigma(2,iGrid)*gxb
         Temp1b=2.0d0*vSigma(3,iGrid)*gxb + vSigma(2,iGrid)*gxa
         Temp2a=2.0d0*vSigma(1,iGrid)*gya + vSigma(2,iGrid)*gyb
         Temp2b=2.0d0*vSigma(3,iGrid)*gyb + vSigma(2,iGrid)*gya
         Temp3a=2.0d0*vSigma(1,iGrid)*gza + vSigma(2,iGrid)*gzb
         Temp3b=2.0d0*vSigma(3,iGrid)*gzb + vSigma(2,iGrid)*gza

         Do iCB = 1, nBfn
*
            Grid_AO(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp0a
     &                             + TabAO1(2,iGrid,iCB) * Temp1a
     &                             + TabAO1(3,iGrid,iCB) * Temp2a
     &                             + TabAO1(4,iGrid,iCB) * Temp3a
            Grid_AO(1,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp0b
     &                             + TabAO1(2,iGrid,iCB) * Temp1b
     &                             + TabAO1(3,iGrid,iCB) * Temp2b
     &                             + TabAO1(4,iGrid,iCB) * Temp3b
         End Do

      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Case Default
        Write (6,*) 'Invalid nD value:', nD
        Call Abend()
      End Select  ! nD
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Case (meta_GGA_type1)
!
!        F(Rho,Sigma,Tau) : Sigma=GradRho*GradRho
!
!        dF/dGradRho = dF/dSigma dSigma/dGradRho = 2 dF/dSigma GradRho
!
!        for the integrals we need:
!
!        phi_i dF/dRho phi_j
!     +  phi_i 2 (dF/dSigma) {GradRho Grad(phi_j)}
!     +  {Grad(phi_i) GradRho} 2 (dF/dSigma) phi_j
!     +  dF/dTau {Grad(phi_i) Grad(phi_j)}
!
!        Grid_AO contains
!        1: 0.5 * phi_i dF/dRho + {Grad(phi_i) GradRho} 2 (dF/dSigma)
!        2: Grad(phi_i)_x dF/dTau
!        3: Grad(phi_i)_y dF/dTau
!        4: Grad(phi_i)_z dF/dTau
!
!        Final integral assembled as, done in do_nIntx.
!
!        Grid_AO(1)_i phi_j
!      + Phi_i * Grid_AO(1)_j
!      + Grid_AO(2)_i Grad(phi_j)_x
!      + Grid_AO(3)_i Grad(phi_j)_y
!      + Grid_AO(4)_i Grad(phi_j)_z
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*
      Select Case(nD)
*                                                                      *
************************************************************************
*                                                                      *
      Case(1)
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, mGrid

!        If (Rho(1,iGrid)<Thr) Cycle
         gx=GradRho(1,iGrid)*Weights(iGrid)
         gy=GradRho(2,iGrid)*Weights(iGrid)
         gz=GradRho(3,iGrid)*Weights(iGrid)

         Temp0=0.5D0*vRho(1,iGrid)*Weights(iGrid)
         Temp1=gx*2.0d0*vSigma(1,iGrid)
         Temp2=gy*2.0d0*vSigma(1,iGrid)
         Temp3=gz*2.0d0*vSigma(1,iGrid)

         Temp4=vTau(1,iGrid)*Weights(iGrid)

         Do iCB = 1, nBfn

            Grid_AO(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp0
     &                             + TabAO1(2,iGrid,iCB) * Temp1
     &                             + TabAO1(3,iGrid,iCB) * Temp2
     &                             + TabAO1(4,iGrid,iCB) * Temp3
            Grid_AO(2,iGrid,iCB,1) = TabAO1(2,iGrid,iCB) * Temp4
            Grid_AO(3,iGrid,iCB,1) = TabAO1(3,iGrid,iCB) * Temp4
            Grid_AO(4,iGrid,iCB,1) = TabAO1(4,iGrid,iCB) * Temp4
         End Do

      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Case(2)
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, mGrid

!        If (Rho(1,iGrid)+Rho(2,iGrid)<Thr) Cycle
         gxa=Gradrho(1,iGrid)*Weights(iGrid)
         gya=Gradrho(2,iGrid)*Weights(iGrid)
         gza=Gradrho(3,iGrid)*Weights(iGrid)
         gxb=Gradrho(4,iGrid)*Weights(iGrid)
         gyb=Gradrho(5,iGrid)*Weights(iGrid)
         gzb=Gradrho(6,iGrid)*Weights(iGrid)

         Temp0a=0.5D0*vRho(1,iGrid) * Weights(iGrid)
         Temp0b=0.5D0*vRho(2,iGrid) * Weights(iGrid)
         Temp1a=2.0d0*vSigma(1,iGrid)*gxa + vSigma(2,iGrid)*gxb
         Temp1b=2.0d0*vSigma(3,iGrid)*gxb + vSigma(2,iGrid)*gxa
         Temp2a=2.0d0*vSigma(1,iGrid)*gya + vSigma(2,iGrid)*gyb
         Temp2b=2.0d0*vSigma(3,iGrid)*gyb + vSigma(2,iGrid)*gya
         Temp3a=2.0d0*vSigma(1,iGrid)*gza + vSigma(2,iGrid)*gzb
         Temp3b=2.0d0*vSigma(3,iGrid)*gzb + vSigma(2,iGrid)*gza
         Temp4a= vTau(1,iGrid)*Weights(iGrid)
         Temp4b= vTau(2,iGrid)*Weights(iGrid)

         Do iCB = 1, nBfn
*
            Grid_AO(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp0a
     &                             + TabAO1(2,iGrid,iCB) * Temp1a
     &                             + TabAO1(3,iGrid,iCB) * Temp2a
     &                             + TabAO1(4,iGrid,iCB) * Temp3a
            Grid_AO(2,iGrid,iCB,1) = TabAO1(2,iGrid,iCB) * Temp4a
            Grid_AO(3,iGrid,iCB,1) = TabAO1(3,iGrid,iCB) * Temp4a
            Grid_AO(4,iGrid,iCB,1) = TabAO1(4,iGrid,iCB) * Temp4a
*
            Grid_AO(1,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp0b
     &                             + TabAO1(2,iGrid,iCB) * Temp1b
     &                             + TabAO1(3,iGrid,iCB) * Temp2b
     &                             + TabAO1(4,iGrid,iCB) * Temp3b
            Grid_AO(2,iGrid,iCB,2) = TabAO1(2,iGrid,iCB) * Temp4b
            Grid_AO(3,iGrid,iCB,2) = TabAO1(3,iGrid,iCB) * Temp4b
            Grid_AO(4,iGrid,iCB,2) = TabAO1(4,iGrid,iCB) * Temp4b
*
         End Do

      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Case Default
        Write (6,*) 'Invalid nD value:', nD
        Call Abend()
      End Select  ! nD
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Case (meta_GGA_type2)
!
!        F(Rho,Sigma,Tau,Lapl) : Sigma=GradRho*GradRho
!
!        dF/dGradRho = dF/dSigma dSigma/dGradRho = 2 dF/dSigma GradRho
!
!        for the integrals we need:
!
!        phi_i dF/dRho phi_j
!     +  phi_i 2 dF/dSigma {GradRho Grad(phi_j)}
!     +  {Grad(phi_i) GradRho} 2 (dF/dSigma) phi_j
!     +  dF/dTau {Grad(phi_i) Grad(phi_j)}
!     +  dF/dLapl Lapl(phi_i) phi_j
!     +  dF/dLapl 2 {Grad(phi_i) Grad(phi_j)}
!     +  dF/dLapl phi_i Lapl(phi_j)
!
!        Grid_AO contains
!        1: 0.5 phi_i dF/dRho + {Grad(phi_i) GradRho} 2 (dF/dSigma)
!          +Lapl(phi_i) dF/dLapl
!        2: Grad(phi_i)_x (dF/dTau + 2 dF/dLapl)
!        3: Grad(phi_i)_y (dF/dTau + 2 dF/dLapl)
!        4: Grad(phi_i)_z (dF/dTau + 2 dF/dLapl)
!
!        Final integral assembled as, done in do_nIntx.
!
!        Grid_AO(1)_i phi_j
!      + phi_i Grid_AO(1)_j
!      + Grid_AO(2)_i Grad(phi_j)_x
!      + Grid_AO(3)_i Grad(phi_j)_y
!      + Grid_AO(4)_i Grad(phi_j)_z
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Select Case(nD)
*                                                                      *
************************************************************************
*                                                                      *
      Case(1)
*                                                                      *
************************************************************************
*                                                                      *

      Do iGrid = 1, mGrid

!        If (Rho(1,iGrid)<Thr) Cycle
         gx=GradRho(1,iGrid)*Weights(iGrid)
         gy=GradRho(2,iGrid)*Weights(iGrid)
         gz=GradRho(3,iGrid)*Weights(iGrid)

         Temp0=0.5D0*vRho(1,iGrid)*Weights(iGrid)
         Temp1=gx*2.0d0*vSigma(1,iGrid)
         Temp2=gy*2.0d0*vSigma(1,iGrid)
         Temp3=gz*2.0d0*vSigma(1,iGrid)

         Temp4=vTau(1,iGrid)*Weights(iGrid)
         Temp5=vLapl(1,iGrid)*Weights(iGrid)
         Temp45=Temp4 + Two * Temp5

         Do iCB = 1, nBfn

            Grid_AO(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp0
     &                             + TabAO1(2,iGrid,iCB) * Temp1
     &                             + TabAO1(3,iGrid,iCB) * Temp2
     &                             + TabAO1(4,iGrid,iCB) * Temp3
     &                             +(TabAO1(5,iGrid,iCB)
     &                              +TabAO1(8,iGrid,iCB)
     &                              +TabAO1(10,iGrid,iCB))*Temp5
            Grid_AO(2,iGrid,iCB,1) = TabAO1(2,iGrid,iCB) * Temp45
            Grid_AO(3,iGrid,iCB,1) = TabAO1(3,iGrid,iCB) * Temp45
            Grid_AO(4,iGrid,iCB,1) = TabAO1(4,iGrid,iCB) * Temp45
         End Do

      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Case(2)
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, mGrid

!        If (Rho(1,iGrid)+Rho(2,iGrid)<Thr) Cycle
         gxa=Gradrho(1,iGrid)*Weights(iGrid)
         gya=Gradrho(2,iGrid)*Weights(iGrid)
         gza=Gradrho(3,iGrid)*Weights(iGrid)
         gxb=Gradrho(4,iGrid)*Weights(iGrid)
         gyb=Gradrho(5,iGrid)*Weights(iGrid)
         gzb=Gradrho(6,iGrid)*Weights(iGrid)

         Temp0a=0.5D0*vRho(1,iGrid) * Weights(iGrid)
         Temp0b=0.5D0*vRho(2,iGrid) * Weights(iGrid)
         Temp1a=2.0d0*vSigma(1,iGrid)*gxa + vSigma(2,iGrid)*gxb
         Temp1b=2.0d0*vSigma(3,iGrid)*gxb + vSigma(2,iGrid)*gxa
         Temp2a=2.0d0*vSigma(1,iGrid)*gya + vSigma(2,iGrid)*gyb
         Temp2b=2.0d0*vSigma(3,iGrid)*gyb + vSigma(2,iGrid)*gya
         Temp3a=2.0d0*vSigma(1,iGrid)*gza + vSigma(2,iGrid)*gzb
         Temp3b=2.0d0*vSigma(3,iGrid)*gzb + vSigma(2,iGrid)*gza
         Temp4a= vTau(1,iGrid)*Weights(iGrid)
         Temp4b= vTau(2,iGrid)*Weights(iGrid)
         Temp5a= vLapl(1,iGrid)*Weights(iGrid)
         Temp5b= vLapl(2,iGrid)*Weights(iGrid)
         Temp45a = Temp4a + Two * Temp5a
         Temp45b = Temp4b + Two * Temp5b

         Do iCB = 1, nBfn
*
            Grid_AO(1,iGrid,iCB,1) = TabAO1(1,iGrid,iCB) * Temp0a
     &                             + TabAO1(2,iGrid,iCB) * Temp1a
     &                             + TabAO1(3,iGrid,iCB) * Temp2a
     &                             + TabAO1(4,iGrid,iCB) * Temp3a
     &                             +(TabAO1(5,iGrid,iCB)
     &                              +TabAO1(8,iGrid,iCB)
     &                              +TabAO1(10,iGrid,iCB))*Temp5a
            Grid_AO(2,iGrid,iCB,1) = TabAO1(2,iGrid,iCB) * Temp45a
            Grid_AO(3,iGrid,iCB,1) = TabAO1(3,iGrid,iCB) * Temp45a
            Grid_AO(4,iGrid,iCB,1) = TabAO1(4,iGrid,iCB) * Temp45a
*
            Grid_AO(1,iGrid,iCB,2) = TabAO1(1,iGrid,iCB) * Temp0b
     &                             + TabAO1(2,iGrid,iCB) * Temp1b
     &                             + TabAO1(3,iGrid,iCB) * Temp2b
     &                             + TabAO1(4,iGrid,iCB) * Temp3b
     &                             +(TabAO1(5,iGrid,iCB)
     &                              +TabAO1(8,iGrid,iCB)
     &                              +TabAO1(10,iGrid,iCB))*Temp5b
            Grid_AO(2,iGrid,iCB,2) = TabAO1(2,iGrid,iCB) * Temp45b
            Grid_AO(3,iGrid,iCB,2) = TabAO1(3,iGrid,iCB) * Temp45b
            Grid_AO(4,iGrid,iCB,2) = TabAO1(4,iGrid,iCB) * Temp45b
         End Do

      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Case Default
        Write (6,*) 'Invalid nD value:', nD
        Call Abend()
      End Select  ! nD
*                                                                      *
************************************************************************
************************************************************************
*
      Case Default
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
      End Select
*                                                                      *
************************************************************************
************************************************************************
      Return
      End
