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
      Subroutine Reset_NQ_Grid
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "nq_info.fh"
#include "grid_on_disk.fh"
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*     Reset the size and the accuracy of the grid to the requested
*     values.
*
      L_Quad=L_Quad_Save
      If (Quadrature(1:3).eq.'LMG') Then
         Threshold=Threshold_Save
      Else
         nR=nR_save
      End If
*
      Crowding     =ThrC
*
      Write (6,*)
      Write (6,*) 'Reset the NQ grid!'
      Write (6,*)
      Call Funi_Print()
*                                                                      *
************************************************************************
*                                                                      *
*     Change the Grid set index
*
      iGrid_Set=Final
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
