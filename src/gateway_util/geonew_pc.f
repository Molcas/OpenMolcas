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
* Copyright (C) 1991,2003, Roland Lindh                                *
************************************************************************
      SubRoutine GeoNew_PC()
************************************************************************
*                                                                      *
* Object: to pick up the geometry from a special file. This will only  *
*         make any difference of there exist a file otherwise SEWARD   *
*         will use the geometry as specified by the standard input     *
*         file.                                                        *
*                                                                      *
*     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
*             University of Lund, SWEDEN                               *
*             March '91                                                *
*                                                                      *
*     Modified to work with point charges. RL 20030507                 *
************************************************************************
      use external_centers
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
      Real*8, Dimension(:), Allocatable :: CN
      Interface
        Subroutine Get_PC_Coord_New(CN,lBuf)
        Real*8, Dimension(:), Allocatable :: CN
        Integer lBuf
        End Subroutine
      End Interface
*
*     Check if there is a data field called 'GeoNewPC'
*
      Call Get_PC_Coord_New(CN,lBuf)
      nAtoms=lbuf/nData_XF
*
*     Quit if the datadfield 'NewGeom' is not available
*
      If ( lBuf.eq.0 ) then
         Return
      End If
*
*     Replace coodinates read in subroutine input
*
      call dcopy_(nAtoms*nData_XF,CN,1,XF,1)
      Write (6,*)
      Write (6,'(A)') '    Point Charge data read from RUNFILE'
      Write (6,*)
*
*     Epilogue, end
*
      Call mma_deallocate(CN)
      Return
      End
