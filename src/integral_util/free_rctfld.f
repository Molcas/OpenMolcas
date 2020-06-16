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
       Subroutine Free_RctFld(iXPolType)
       use PCM_arrays
       Implicit Real*8 (a-h,o-z)
#include "rctfld.fh"
#include "status.fh"
#include "stdalloc.fh"
*
       If (RctFld_Status.eq.InActive) Return
       Call GetMem('MM','Free','Real',ipMM,nMM)
       If (lLangevin .or. (iXPolType.gt.0)) Then
          If(iXPolType.eq.2) Then
             nPolComp = 6
          Else
             nPolComp = 1
          EndIf
          Call GetMem('Field ','Free','Real',ipField ,nGrid*4)
          Call GetMem('dField','Free','Real',ipdField,nGrid*4)
          Call GetMem('Dip   ','Free','Real',ipDip   ,nGrid*3)
          Call GetMem('PolEf ','Free','Real',ipPolEf ,nGrid*nPolComp)
          Call GetMem('DipEf ','Free','Real',ipDipEf ,nGrid  )
          Call GetMem('Grid  ','Free','Real',ipGrid  ,nGrid*3)
*
          Call GetMem('favxyz','Free','Real',ipfavxyz,nCavxyz)
          Call GetMem('davxyz','Free','Real',ipdavxyz,nCavxyz)
          Call GetMem('cavxyz','Free','Real',ipcavxyz,nCavxyz)
          Call GetMem('ravxyz','Free','Real',ipravxyz,nCavxyz)
       End If
       If (PCM) Then
          Call GetMem('PCMSph','Free','Real',ip_Sph,nPCM_info)
          Call mma_deallocate(Vert)
          Call mma_deallocate(PCMTess)
          Call mma_deallocate(PCMSph)
*
*---- Free the space for geometric derivatives
*
          If (DoDeriv) Then
             Call mma_deallocate(dTes)
             Call mma_deallocate(dPnt)
             Call mma_deallocate(dRad)
             Call mma_deallocate(dCntr)
             Call mma_deallocate(PCM_SQ)
          End If
*
       End If
       RctFld_Status=InActive
*
       Return
       End
