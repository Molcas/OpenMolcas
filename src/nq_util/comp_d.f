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
* Copyright (C) Giovanni Li Manni                                      *
************************************************************************
      Real*8 Function Comp_d(Weights,mGrid,Rho,nRho,iSpin,T_X,iSwitch)
************************************************************************
*                                                                      *
* Object: integrate densities (alpha, beta, total, gradients....)      *
*         the object integrated is dictaded by iSwitch value:          *
*         iSwitch = 0  .... total density                              *
*         iSwitch = 1  .... alpha density                              *
*         iSwitch = 2  .... beta density                               *
*                                                                      *
* Author: G. Li Manni... taking Sir R. Lindh as model                  *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 Weights(mGrid), Rho(nRho,mGrid)
*
      Comp_d=Zero
      Rho_min=T_X*1.0D-2
************************************************************************
*     iSpin=1
************************************************************************
      If (iSpin.eq.1) Then
         Do iGrid = 1, mGrid
           d_alpha=Rho(1,iGrid)
           if(iSwitch.eq.1) then
             DTot=d_alpha
           else if(iSwitch.eq.2) then
             DTot=d_alpha
           else !if(iSwitch.eq.0) then
             DTot=Two*d_alpha
           end if
           If (DTot.lt.T_X) Go To 199
           Comp_d = Comp_d + DTot*Weights(iGrid)
 199       Continue
         End Do
************************************************************************
*     iSpin=/=1
************************************************************************
      Else
        Do iGrid = 1, mGrid
          d_alpha=Max(Rho_min,Rho(1,iGrid))
          d_beta =Max(Rho_min,Rho(2,iGrid))
          if(iSwitch.eq.1) then
            DTot=d_alpha
          else if(iSwitch.eq.2) then
            DTot=d_beta
          else !if(iSwitch.eq.0) then
            DTot=d_alpha+d_beta
          end if
          If (DTot.lt.T_X) Go To 299
          Comp_d = Comp_d + DTot*Weights(iGrid)
 299      Continue
        End Do
      End If
************************************************************************
      Return
      End
