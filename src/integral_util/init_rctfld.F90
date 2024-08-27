!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Init_RctFld(NonEq,iCharge)

use Langevin_arrays, only: Cavxyz, Davxyz, Ravxyz, DField, Dip, DipEF, Field, Grid, PolEF
use external_centers, only: nXF, iXPolType
use stdalloc, only: mma_allocate
use rctfld_module, only: TK, lMax, nMM, nGrid, lLangevin, MaxA, RadLat, Scala, MaxB, Scalb, MaxC, Scalc, nABC, lAtAto, PCM, &
                         NonEQ_Ref, nCavxyz, MM

implicit none
logical NonEq
integer iCharge
integer MMM, nPolComp

tK = 1.0D-99 ! Boltzman factor, initial set to 0 K
if (allocated(MM)) return
mMM = (lMax+1)*(lMax+2)*(lMax+3)/6
nMM = 2*mMM
call mma_allocate(MM,mMM,2,Label='MM')
MM(:,:) = 0.0d0
if (iXPolType > 0) nGrid = nXF
if (lLangevin .or. (iXPolType > 0)) then
  if (lLangevin) then
    maxa = int(radlat/scala)
    maxb = int(radlat/scalb)
    maxc = int(radlat/scalc)
    nabc = (2*maxa+2)*(2*maxb+2)*(2*maxc+2)
    nGrid = nGrid+nabc*latato
  end if
  if (iXPolType == 2) then
    nPolComp = 6
  else
    nPolComp = 1
  end if
  call mma_allocate(Field,4,nGrid,Label='Field')
  call mma_allocate(dField,4,nGrid,Label='dField')
  call mma_allocate(Dip,3,nGrid,Label='Dip')
  call mma_allocate(PolEf,nPolComp,nGrid,Label='PolEf')
  call mma_allocate(DipEf,nGrid,Label='DipEf')
  call mma_allocate(Grid,3,nGrid,Label='Grid')

  nCavxyz = (lMax+1)*(lMax+2)*(lMax+3)/6
  call mma_allocate(davxyz,nCavxyz,Label='davxyz')
  call mma_allocate(cavxyz,nCavxyz,Label='cavxyz')
  call mma_allocate(ravxyz,nCavxyz,Label='ravxyz')
end if
if (.not. PCM) NonEq_Ref = NonEq
call Init_PCM(NonEq,iCharge)

return

end subroutine Init_RctFld
