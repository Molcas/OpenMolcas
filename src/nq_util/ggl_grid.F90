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
      Subroutine GGL_Grid(L_Max)
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
         Subroutine Do_GGL(L_Eff,nPoints,R)
         Implicit None
         Integer L_Eff, nPoints
         Real*8, Allocatable:: R(:,:)
         End Subroutine Do_GGL
      End Interface
!                                                                      *
!***********************************************************************
!                                                                      *
!     Generate angular grid from Gauss and Gauss-Legendre quadrature
!
!---- Theta (polar angle): 0 =< theta =< pi
!     Gauss-Legendre Quadrature (L_Quad+1)/2 points
!---- Phi (azimuthal angle): 0=< phi =< 2*pi
!     Gauss-Quadrature (L_Quad+1) points
!
      Do L_Eff = 1, L_Max
         nAngularGrids=nAngularGrids+1
!
         Info_Ang(nAngularGrids)%L_eff=L_eff
         Call Do_GGL(L_Eff,                                             &
     &               Info_Ang(nAngularGrids)%nPoints,                   &
     &               Info_Ang(nAngularGrids)%R)
!
      End Do          ! L_Eff
!                                                                      *
!***********************************************************************
!                                                                      *
      Return
      End
