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
* Copyright (C) 2000, Roland Lindh                                     *
************************************************************************
      Subroutine DiracX(mGrid,iSpin,F_xc,Coeff,T_X)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. November 2000                           *
************************************************************************
C-Ajitha Modifying the kernel output structure
      use KSDFT_Info, only: F_xca, F_xcb
      use nq_Grid, only: Rho, l_casdft
      use nq_Grid, only: vRho
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
#include "ksdft.fh"
      Real*8 F_xc(mGrid)
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
      THIRD =One/Three
      FTHIRD=Four/Three
      CVX   =(two**THIRD)*((three/Pi)**THIRD)
      Rho_min=T_X*1.0D-2
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute value of energy and integrad on the grid
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin=1
*
      If (iSpin.eq.1) Then
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, mGrid
         d_alpha =Rho(1,iGrid)
         DTot=Two*d_alpha
         If (DTot.lt.T_X) Go To 100
*
*------- Exchange contributions to energy
*
         functional=-Three/Two*CVX*d_alpha**FTHIRD
         F_xc(iGrid)=F_xc(iGrid)+Coeff*functional
*
*------- Exchange contributions to the AO integrals
*
         func_d_rho_alpha=-CVX*d_alpha**THIRD
*

         vRho(1,iGrid) = vRho(1,iGrid)
     &                               + Coeff*func_d_rho_alpha
*
 100     Continue
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin=/=1
*
      Else
*                                                                      *
************************************************************************
*                                                                      *
      If (l_casdft) Then
      Do iGrid = 1, mGrid
         d_alpha =Max(Rho_Min,Rho(1,iGrid))
         d_beta  =Max(Rho_Min,Rho(2,iGrid))
         DTot=d_alpha+d_beta
         If (DTot.lt.T_X) Cycle
*------- Exchange contributions to energy
*
         functional =-Three/Four*CVX*(d_alpha**FTHIRD+d_beta**FTHIRD)
         functionala=-Three/Four*CVX*(d_alpha**FTHIRD)
         functionalb=-Three/Four*CVX*(                d_beta**FTHIRD)
         F_xc(iGrid) =F_xc(iGrid) +Coeff*functional
         F_xca(iGrid)=F_xca(iGrid)+Coeff*functionala
         F_xcb(iGrid)=F_xcb(iGrid)+Coeff*functionalb
*
*------- Exchange contributions to the AO integrals
*
         func_d_rho_alpha=-CVX*d_alpha**THIRD
         func_d_rho_beta =-CVX*d_beta **THIRD
*
         vRho(1,iGrid) = vRho(1,iGrid)
     &                               + Coeff*func_d_rho_alpha
         vRho(2,iGrid) = vRho(2,iGrid)
     &                               + Coeff*func_d_rho_beta
*
      End Do
      Else
      Do iGrid = 1, mGrid
         d_alpha =Max(Rho_Min,Rho(1,iGrid))
         d_beta  =Max(Rho_Min,Rho(2,iGrid))
         DTot=d_alpha+d_beta
         If (DTot.lt.T_X) Cycle
*------- Exchange contributions to energy
*
         functional =-Three/Four*CVX*(d_alpha**FTHIRD+d_beta**FTHIRD)
         F_xc(iGrid) =F_xc(iGrid) +Coeff*functional
*
*------- Exchange contributions to the AO integrals
*
         func_d_rho_alpha=-CVX*d_alpha**THIRD
         func_d_rho_beta =-CVX*d_beta **THIRD
*
         vRho(1,iGrid) = vRho(1,iGrid)
     &                               + Coeff*func_d_rho_alpha
         vRho(2,iGrid) = vRho(2,iGrid)
     &                               + Coeff*func_d_rho_beta
*
      End Do
      End If
*
      End If
*
      Return
      End
