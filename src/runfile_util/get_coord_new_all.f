************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Roland Lindh                                           *
************************************************************************
      Subroutine Get_Coord_New_All(Coord_All,nAtoms_All)
************************************************************
*
*   <DOC>
*     <Name>Get\_Coord\_All</Name>
*     <Syntax>Call Get\_Coord\_New\_All(Coord\_All,nAtoms\_All)</Syntax>
*     <Arguments>
*       \Argument{Coord\_All}{Array of coordinates}{Real*8 (3,nAtoms\_All)}{out}
*       \Argument{nAtoms\_All}{Number of atoms}{Integer}{in}
*     </Arguments>
*     <Purpose>Get the new Coordinates from RUNFILE</Purpose>
*     <Dependencies></Dependencies>
*     <Author>R. Lindh</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*       Place cartesian coordinates (in a.u. into array Coord\_All(3,*)
*     </Description>
*    </DOC>
*
************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Real*8 Coord_All(3,nAtoms_All)
      Real*8, Dimension (:,:), Allocatable :: CU
      Interface
        Subroutine Get_Coord_New(CU,nAtoms)
        Real*8, Dimension (:,:), Allocatable :: CU
        Integer nAtoms
        End Subroutine
      End Interface
*
      Call Get_nAtoms_All(nAtoms_Allx)
      If (nAtoms_All.ne.nAtoms_Allx) Then
         Write (6,*) 'Get_Coord_All: nAtoms_All.ne.nAtoms_Allx'
         Write (6,*) 'nAtoms_All=',nAtoms_All
         Write (6,*) 'nAtoms_Allx=',nAtoms_Allx
         Call QTrace
         Call Abend
      End If
      Call Get_Coord_New(CU,nAtoms)
      Call Get_Coord_All_(CU,nAtoms,Coord_All,nAtoms_All)
      Call mma_deallocate(CU)
*
      Return
      End
