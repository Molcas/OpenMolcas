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
      Subroutine GGL_Grid(L_Max)
************************************************************************
*                                                                      *
*     Computes datas useful for the angular quadrature.                *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "nq_info.fh"
#include "real.fh"
#include "WrkSpc.fh"
*                                                                      *
************************************************************************
*                                                                      *
*     Generate angular grid from Gauss and Gauss-Legendre quadrature
*
*---- Theta (polar angle): 0 =< theta =< pi
*     Gauss-Legendre Quadrature (L_Quad+1)/2 points
*---- Phi (azimuthal angle): 0=< phi =< 2*pi
*     Gauss-Quadrature (L_Quad+1) points
*
      Do L_Eff = 1, L_Max
         nAngularGrids=nAngularGrids+1
*
         Call Do_GGL(L_Eff,nPoints,ipR)
*
         Info_Ang(1,nAngularGrids)=L_Eff
         Info_Ang(2,nAngularGrids)=nPoints
         Info_Ang(3,nAngularGrids)=ipR
*
      End Do          ! L_Eff
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
