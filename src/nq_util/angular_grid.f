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
      Subroutine Angular_Grid()
************************************************************************
*                                                                      *
*     Computes datas useful for the angular quadrature.                *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "nq_info.fh"
#include "real.fh"
#include "debug.fh"
#include "WrkSpc.fh"
      Logical Check
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions                                              *
*                                                                      *
      Check(i,j)=iAnd(i,2**(j-1)).ne.0
*                                                                      *
************************************************************************
*                                                                      *
      nAngularGrids=0
      If (Check(iOpt_Angular,3)) Then
*                                                                      *
************************************************************************
*                                                                      *
*------- Generate angular grid a la Lebedev
*
         Call Lebedev_Grid(L_Quad)
*                                                                      *
************************************************************************
*                                                                      *
      Else If (Check(iOpt_Angular,1)) Then
*                                                                      *
************************************************************************
*                                                                      *
*------- Generate angular grid a la Lobatto
*
         Call Lobatto_Grid(L_Quad)
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
*------- Generate angular grid from Gauss and Gauss-Legendre quadrature
*
         Call GGL_Grid(L_Quad)
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
      If (Debug) Then
         Do iSet = 1, nAngularGrids
            iR =Info_Ang(3,iSet)
            nGP=Info_Ang(2,iSet)
            l  =Info_Ang(1,iSet)
            Write (6,*) 'l=',l
            Call RecPrt('Angular grid',' ',Work(iR),4,nGP)
         End Do
      End If
*
      Return
      End
