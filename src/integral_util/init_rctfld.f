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
       Subroutine Init_RctFld(NonEq,iCharge)
       Implicit Real*8 (a-h,o-z)
#include "rctfld.fh"
#include "status.fh"
#include "WrkSpc.fh"
#include "itmax.fh"
#include "info.fh"
       Logical NonEq
*
       If (RctFld_Status.eq.Active) Return
       mMM = (lMax+1)*(lMax+2)*(lMax+3)/6
       nMM = 2 * mMM
       Call GetMem('MM','Allo','Real',ipMM,nMM)
       If (iXPolType.gt.0) nGrid = nXF
       If (lLangevin .or. (iXPolType.gt.0)) Then
          If(lLangevin) Then
             maxa = INT(radlat/scala)
             maxb = INT(radlat/scalb)
             maxc = INT(radlat/scalc)
             nabc=(2*maxa+2)*(2*maxb+2)*(2*maxc+2)
             nGrid=nGrid+nabc*latato
          EndIf
          If(iXPolType.eq.2) Then
             nPolComp = 6
          Else
             nPolComp = 1
          EndIf
          Call GetMem('Field ','Allo','Real',ipField ,nGrid*4)
          Call GetMem('dField','Allo','Real',ipdField,nGrid*4)
          Call GetMem('Dip   ','Allo','Real',ipDip   ,nGrid*3)
          Call GetMem('PolEf ','Allo','Real',ipPolEf ,nGrid*nPolComp)
          Call GetMem('DipEf ','Allo','Real',ipDipEf ,nGrid  )
          Call GetMem('Grid  ','Allo','Real',ipGrid  ,nGrid*3)
*
          nCavxyz = (lMax+1)*(lMax+2)*(lMax+3)/6
          Call GetMem('favxyz','Allo','Real',ipfavxyz,nCavxyz)
          Call GetMem('davxyz','Allo','Real',ipdavxyz,nCavxyz)
          Call GetMem('cavxyz','Allo','Real',ipCavxyz,nCavxyz)
          Call GetMem('ravxyz','Allo','Real',ipravxyz,nCavxyz)
       End If
       If (.Not.PCM) NonEq_Ref=NonEq
       Call Init_PCM(NonEq,iCharge)
       RctFld_Status=Active
*
       Return
       End
