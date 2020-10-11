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
      Subroutine DiracX_OFE(mGrid,Rho,nRho,iSpin,F_xc,dF_dRho,
     &                      ndF_dRho,Coeff,T_X)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. November 2000                           *
************************************************************************
C-Ajitha Modifying the kernel output structure
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),F_xc(mGrid)
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
         d_alpha =Rho(ipR,iGrid)
         DTot=Two*d_alpha
         If (DTot.lt.T_X) Go To 100
*
*------- Exchange contributions to energy
*
         functional=-Three/Two*CVX*d_alpha**FTHIRD
         F_xc(iGrid)=F_xc(iGrid)+functional
*
*------- Exchange contributions to the AO intergrals
*
         func_d_rho_alpha=-CVX*d_alpha**THIRD
*

        dF_dRho(ipR,iGrid) = dF_dRho(ipR,iGrid)
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
      Do iGrid = 1, mGrid
         d_alpha =Max(Rho_Min,Rho(ipRa,iGrid))
         d_beta  =Max(Rho_Min,Rho(ipRb,iGrid))
         DTot=d_alpha+d_beta
         If (DTot.lt.T_X) Go To 200
*------- Exchange contributions to energy
*
         functional=-Three/Four*CVX*(d_alpha**FTHIRD+d_beta**FTHIRD)
         F_xc(iGrid)=F_xc(iGrid)+functional
*
*------- Exchange contributions to the AO intergrals
*
         func_d_rho_alpha=-CVX*d_alpha**THIRD
         func_d_rho_beta =-CVX*d_beta **THIRD
*
         dF_dRho(ipRa,iGrid) = dF_dRho(ipRa,iGrid)
     &                               + Coeff*func_d_rho_alpha
         dF_dRho(ipRb,iGrid) = dF_dRho(ipRb,iGrid)
     &                               + Coeff*func_d_rho_beta
*
 200     Continue
*
      End Do
*
      End If
*
      Return
      End
