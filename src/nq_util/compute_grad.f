!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2000, Roland Lindh                                     *
!***********************************************************************
      Real*8 Function Compute_Grad(Weights,mGrid,iSpin)
!***********************************************************************
!      Author:Roland Lindh, Department of Chemical Physics, University *
!             of Lund, SWEDEN. November 2000                           *
!***********************************************************************
      use nq_Grid, only: Sigma
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 Weights(mGrid)
!                                                                      *
!***********************************************************************
!                                                                      *
!
      Compute_Grad=Zero
!
!     iSpin=1
!
      If (iSpin.eq.1) Then
!                                                                      *
!***********************************************************************
!                                                                      *
      Do iGrid = 1, mGrid
!
         Gamma=Sqrt(Sigma(1,iGrid))
!
!------- Accumulate contributions to the integrated Tau
!
         Compute_Grad = Compute_Grad + Two*Gamma*Weights(iGrid)
!
      End Do
!                                                                      *
!***********************************************************************
!                                                                      *
!     iSpin=/=1
!
      Else
!                                                                      *
!***********************************************************************
!                                                                      *
      Do iGrid = 1, mGrid
!
         Gamma=Sqrt(Sigma(1,iGrid)+Two*Sigma(2,iGrid)+Sigma(3,iGrid))
!
!------- Accumulate contributions to the integrated density
!
         Compute_Grad = Compute_Grad + Gamma*Weights(iGrid)
!
      End Do
!                                                                      *
!***********************************************************************
!                                                                      *
      End If
!                                                                      *
!***********************************************************************
!                                                                      *
!
      Return
      End
