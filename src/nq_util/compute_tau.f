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
      Real*8 Function Compute_Tau(Weights,mGrid,iSpin)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. November 2000                           *
************************************************************************
      use nq_Grid, only: Rho, Tau
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 Weights(mGrid)
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*
      Compute_Tau=Zero
*
*     iSpin=1
*
      If (iSpin.eq.1) Then
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, mGrid
*
         d_alpha=Rho(1,iGrid)
         DTot=d_alpha
         TauA=Two*Tau(1,iGrid)
*
*------- Accumulate contributions to the integrated Tau
*
         Compute_Tau = Compute_Tau + TauA*Weights(iGrid)
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
*
         d_alpha=Rho(1,iGrid)
         d_beta =Rho(2,iGrid)
         DTot=d_alpha+d_beta
         TauA=(Tau(1,iGrid)+Tau(2,iGrid))
*
*------- Accumulate contributions to the integrated density
*
         Compute_Tau = Compute_Tau + TauA*Weights(iGrid)
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
      Return
      End
