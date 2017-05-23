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
      Implicit Real*8 (a-h,o-z)
#include "nq_info.fh"
#include "real.fh"
#include "WrkSpc.fh"
      Parameter (nSet=11)
      Integer Lebedev_order(nSet)
      Data Lebedev_order/5,7,11,17,23,29,35,41,47,53,59/
*                                                                      *
************************************************************************
*                                                                      *
*     Use the GGL grids to minimize the number of grid points.
*
      If (L_Max.lt.3) Go To 99
      Call Do_GGL(3,nPoints,ipR)
      nAngularGrids=nAngularGrids+1
      Info_Ang(1,nAngularGrids)=3
      Info_Ang(2,nAngularGrids)=nPoints
      Info_Ang(3,nAngularGrids)=ipR
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
            Call Do_Lebedev(L_Eff,nPoints,ipR)
*
            Info_Ang(1,nAngularGrids)=L_Eff
            Info_Ang(2,nAngularGrids)=nPoints
            Info_Ang(3,nAngularGrids)=ipR
         Else
*
            Go To 99
*
         End If
      End Do
*                                                                      *
************************************************************************
*                                                                      *
 99   Continue
      Return
      End
