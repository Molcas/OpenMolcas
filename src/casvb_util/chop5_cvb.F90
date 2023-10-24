!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine chop5_cvb()

use casvb_global, only: corth, iconstruc, idelstr, ifxorb, ifxstr, iorbrel, iorts, ipermzeta, irels, irots, izeta, ndimrel, ndrot, &
                        norb, norbrel, nort, north, nprorb, nsyme, nvb, nzeta, nzrvb, orbfr_is_unit, release, relorb, symelm, &
                        tconstr, trprm
use stdalloc, only: mma_allocate, mma_deallocate

implicit none

if (release(5)) then
  call mma_deallocate(symelm)
  call mma_deallocate(iorbrel)
  call mma_deallocate(north)
  call mma_deallocate(corth)
  call mma_deallocate(irels)
  call mma_deallocate(relorb)
  call mma_deallocate(ifxorb)
  call mma_deallocate(ifxstr)
  call mma_deallocate(idelstr)
  call mma_deallocate(iorts)
  call mma_deallocate(irots)
  call mma_deallocate(izeta)
  call mma_deallocate(trprm)
  call mma_deallocate(tconstr)
  call mma_deallocate(ipermzeta)
end if
release(5) = .true.
release(6) = .false.

!FIXME: These deallocations should not be needed
if (allocated(symelm)) call mma_deallocate(symelm)
if (allocated(iorbrel)) call mma_deallocate(iorbrel)
if (allocated(north)) call mma_deallocate(north)
if (allocated(corth)) call mma_deallocate(corth)
if (allocated(irels)) call mma_deallocate(irels)
if (allocated(relorb)) call mma_deallocate(relorb)
if (allocated(ifxorb)) call mma_deallocate(ifxorb)
if (allocated(ifxstr)) call mma_deallocate(ifxstr)
if (allocated(idelstr)) call mma_deallocate(idelstr)
if (allocated(iorts)) call mma_deallocate(iorts)
if (allocated(irots)) call mma_deallocate(irots)
if (allocated(izeta)) call mma_deallocate(izeta)
if (allocated(trprm)) call mma_deallocate(trprm)
if (allocated(tconstr)) call mma_deallocate(tconstr)
if (allocated(ipermzeta)) call mma_deallocate(ipermzeta)
call mma_allocate(symelm,norb,norb,nsyme,label='symelm')
call mma_allocate(iorbrel,ndimrel,label='iorbrel')
call mma_allocate(north,norb,label='north')
call mma_allocate(corth,norb,norb**2,label='corth')
call mma_allocate(irels,2,norb-1,label='irels')
call mma_allocate(relorb,norb,norb,min(norb-1,norbrel),label='relorb')
call mma_allocate(ifxorb,norb,label='ifxorb')
call mma_allocate(ifxstr,nvb,label='ifxstr')
call mma_allocate(idelstr,nzrvb,label='idelstr')
call mma_allocate(iorts,2,nort,label='iorts')
call mma_allocate(irots,2,ndrot,label='irots')
call mma_allocate(izeta,nsyme,label='izeta')
if (.not. orbfr_is_unit) then
  call mma_allocate(trprm,nprorb,nprorb,label='trprm')
else
  call mma_allocate(trprm,0,0,label='trprm')
end if
if (iconstruc == 2) then
  call mma_allocate(tconstr,nvb,nvb,label='tconstr')
else
  call mma_allocate(tconstr,0,0,label='tconstr')
end if
call mma_allocate(ipermzeta,norb,nzeta,label='ipermzeta')

return

end subroutine chop5_cvb
