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

subroutine Free_RctFld()

use PCM_arrays, only: Centr, dCntr, dPnt, dRad, dTes, IntSph, NewSph, nVert, PCM_N, PCM_SQ, PCMDm, PCMiSph, PCMSph, PCMTess, SSph, &
                      Vert
use Langevin_arrays, only: Cavxyz, Davxyz, dField, Dip, DipEF, Field, Grid, PolEF, Ravxyz
use rctfld_module, only: DoDeriv, lLangevin, MM, PCM
use External_Centers, only: iXPolType
use stdalloc, only: mma_deallocate

implicit none

if (allocated(MM)) then

  call mma_deallocate(MM)

  if (lLangevin .or. (iXPolType > 0)) then
    call mma_deallocate(Field)
    call mma_deallocate(dField)
    call mma_deallocate(Dip)
    call mma_deallocate(PolEf)
    call mma_deallocate(DipEf)
    call mma_deallocate(Grid)
    call mma_deallocate(davxyz)
    call mma_deallocate(cavxyz)
    call mma_deallocate(ravxyz)
  end if

  if (PCM) then
    call mma_deallocate(NewSph)
    call mma_deallocate(IntSph)
    call mma_deallocate(NVert)
    call mma_deallocate(PCMiSph)
    call mma_deallocate(PCM_N)
    call mma_deallocate(PCMDM)
    call mma_deallocate(SSph)
    call mma_deallocate(Centr)
    call mma_deallocate(Vert)
    call mma_deallocate(PCMTess)
    call mma_deallocate(PCMSph)

    ! Free the space for geometric derivatives

    if (DoDeriv) then
      call mma_deallocate(dTes)
      call mma_deallocate(dPnt)
      call mma_deallocate(dRad)
      call mma_deallocate(dCntr)
      call mma_deallocate(PCM_SQ)
    end if

  end if

end if

return

end subroutine Free_RctFld
