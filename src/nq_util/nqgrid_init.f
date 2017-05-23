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
      Subroutine NQGrid_Init()
      Implicit Real*8 (A-H,O-Z)
#include "grid_on_disk.fh"
      Integer Other_Type
      Parameter (Other_Type=2)
*                                                                      *
************************************************************************
*                                                                      *
*     Make the grid file dirty
*
*---- Open the file.
      Lu_Grid=77
      Call DaName_MF_WA(Lu_Grid,'NQGRID')
*
*---- Write the status flag and disk addresses fo the sets.
*
      iDisk_Set(Final)=-1
      iDisk_Set(Intermediate)=-1
      G_S(Final)=Regenerate
      G_S(Intermediate)=Regenerate
      Old_Functional_Type=Other_Type
*
      iDisk_Grid=0
      Call iDaFile(Lu_Grid,1,G_S,5,iDisk_Grid)
*
      iDisk_Set(Final)=iDisk_Grid
      iDisk_Set(Intermediate)=iDisk_Grid
*
      iDisk_Grid=0
      Call iDaFile(Lu_Grid,1,G_S,5,iDisk_Grid)
*
      Call DaClos(Lu_Grid)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
