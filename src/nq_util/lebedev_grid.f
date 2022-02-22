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
      Subroutine Lebedev_Grid(L_Max)
************************************************************************
*                                                                      *
*     Computes datas useful for the angular quadrature.                *
*                                                                      *
************************************************************************
      use nq_Structure, only: Info_Ang
      use nq_Info
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Parameter (nSet=11)
      Integer Lebedev_order(nSet)
      Data Lebedev_order/5,7,11,17,23,29,35,41,47,53,59/
*                                                                      *
************************************************************************
*                                                                      *
      Interface
         Subroutine Do_GGL(L_Eff,nPoints,R)
         Implicit None
         Integer L_Eff, nPoints
         Real*8, Allocatable:: R(:,:)
         End Subroutine Do_GGL
         Subroutine Do_Lebedev(L_Eff,nPoints,R)
         Implicit None
         Integer L_Eff, nPoints
         Real*8, Allocatable:: R(:,:)
         End Subroutine Do_Lebedev
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
*     Use the GGL grids to minimize the number of grid points.
*
      If (L_Max.lt.3) Return
      nAngularGrids=nAngularGrids+1
      Info_Ang(nAngularGrids)%L_Eff=3
      Call Do_GGL(3,
     &            Info_Ang(nAngularGrids)%nPoints,
     &            Info_Ang(nAngularGrids)%R)
*                                                                      *
************************************************************************
*                                                                      *
*---- Generate angular grid a la Lebedev
*
      Do iSet = 1, nSet
         If (Lebedev_order(iSet).le.L_Max) Then
            nAngularGrids=nAngularGrids+1
            L_Eff=Lebedev_order(iSet)
*
            Info_Ang(nAngularGrids)%L_Eff=L_Eff
            Call Do_Lebedev(L_Eff,
     &                      Info_Ang(nAngularGrids)%nPoints,
     &                      Info_Ang(nAngularGrids)%R)
         Else
*
            Return
*
         End If
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
