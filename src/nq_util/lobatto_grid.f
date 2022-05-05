!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine Lobatto_Grid(L_Max)
!***********************************************************************
!                                                                      *
!     Computes datas useful for the angular quadrature.                *
!                                                                      *
!***********************************************************************
      use nq_Structure, only: Info_Ang
      use nq_Info
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
!                                                                      *
!***********************************************************************
!                                                                      *
      Interface
         Subroutine Do_Lobatto(L_Eff,nPoints,R)
         Implicit None
         Integer L_Eff, nPoints
         Real*8, Allocatable:: R(:,:)
         End Subroutine Do_Lobatto
      End Interface
!                                                                      *
!***********************************************************************
!                                                                      *
!     Observe that we use standard GGL for orders 1 and 2.
      Call GGL_Grid(2)
!                                                                      *
!***********************************************************************
!                                                                      *
!---- Generate angular grid a la Lobatto
!
      Do L_Eff = 3, L_Max
         nAngularGrids=nAngularGrids+1
!
         Info_Ang(nAngularGrids)%L_Eff=L_Eff
         Call Do_Lobatto(L_Eff,                                         &
     &                   Info_Ang(nAngularGrids)%nPoints,               &
     &                   Info_Ang(nAngularGrids)%R)
!
      End Do        ! L_Eff
!                                                                      *
!***********************************************************************
!                                                                      *
      Return
      End
