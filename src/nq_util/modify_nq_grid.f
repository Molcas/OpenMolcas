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
      Subroutine Modify_NQ_grid
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "nq_info.fh"
#include "grid_on_disk.fh"
      Parameter (L_Quad_Low=23, Threshold_High=1.0D-7, nR_Low=50)
*                                                                      *
************************************************************************
*                                                                      *
      Call qEnter('Modify')
*                                                                      *
************************************************************************
*                                                                      *
*     Reduce the size and the accuracy of the grid temporarily
*
      L_Quad_Save=L_Quad
      Threshold_save=Threshold
      nR_Save=nR
      ThrC = Crowding
*
      L_Quad=Min(L_Quad,L_Quad_Low)
      If (Quadrature(1:3).eq.'LMG') Then
         Threshold     =Max(Threshold_High,Threshold)
      Else
         nR     =Min(nR_Low,nR)
      End If
      Crowding=Max(ThrC-Two,One)
*
      Write (6,*)
      Write (6,*) 'Modify the NQ grid!'
      Write (6,*)
      Call Funi_Print()
*                                                                      *
************************************************************************
*                                                                      *
*     Change the Grid set index.
*
      iGrid_Set=Intermediate
*                                                                      *
************************************************************************
*                                                                      *
      Call qExit('Modify')
      Return
      End
