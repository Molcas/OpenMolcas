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
       use Langevin_arrays
       use PCM_arrays, only: MM
       use external_centers, only: nXF, iXPolType
       Implicit Real*8 (a-h,o-z)
#include "rctfld.fh"
#include "status.fh"
#include "stdalloc.fh"
#include "itmax.fh"
#include "info.fh"
       Logical NonEq
*
       tK=1.0D-99 ! Boltzman factor, initial set to 0 K
       If (RctFld_Status.eq.Active) Return
       mMM = (lMax+1)*(lMax+2)*(lMax+3)/6
       nMM = 2 * mMM
       Call mma_allocate(MM,mMM,2,Label='MM')
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
          Call mma_allocate(Field,4,nGrid,Label='Field')
          Call mma_allocate(dField,4,nGrid,Label='dField')
          Call mma_allocate(Dip,3,nGrid,Label='Dip')
          Call mma_allocate(PolEf,nPolComp,nGrid,Label='PolEf')
          Call mma_allocate(DipEf,nGrid,Label='DipEf')
          Call mma_allocate(Grid,3,nGrid,Label='Grid')
*
          nCavxyz = (lMax+1)*(lMax+2)*(lMax+3)/6
          Call mma_allocate(davxyz,nCavxyz,Label='davxyz')
          Call mma_allocate(cavxyz,nCavxyz,Label='cavxyz')
          Call mma_allocate(ravxyz,nCavxyz,Label='ravxyz')
       End If
       If (.Not.PCM) NonEq_Ref=NonEq
       Call Init_PCM(NonEq,iCharge)
       RctFld_Status=Active
*
       Return
       End
