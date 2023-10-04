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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************
      SubRoutine CloseP()
!***********************************************************************
!                                                                      *
! Object: to close the handling of the 2nd order density matrix.       *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************
      use pso_stuff, only: Case_MP2, LuGam, lPSO, G1, G2, CMO, DSVar,
     &                     DS, DVar, D0, Gamma_On, LuGamma, Bin, G_Toc,
     &                     SO2cI
      use stdalloc, only: mma_deallocate
      Implicit None
      Logical DoCholesky
!
      If(case_mp2) then
         Call DecideOnCholesky(DoCholesky)
         If(.not. DoCholesky) Then
            Call DaClos(LuGam)
         End If
      End If
      If (Gamma_On) Then
!*********** columbus interface ****************************************
         Call DaClos(LuGamma)
         Call mma_deallocate(Bin)
         Call mma_deallocate(G_Toc)
         Call mma_deallocate(SO2cI)
      End If
!
      If (lPSO) Then
         Call mma_deallocate(G2)
         Call mma_deallocate(G1)
      End If
      call mma_deallocate(CMO)
      call mma_deallocate(DSVar)
      call mma_deallocate(DS)
      call mma_deallocate(DVar)
      call mma_deallocate(D0)
!
      Return
      End SubRoutine CloseP
