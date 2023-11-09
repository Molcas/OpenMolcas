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

subroutine change5_cvb()

use casvb_global, only: iconstruc, imethod, lfxvb, lzrvb, ndimrel, ndrot, nfxorb, nfxvb, norbrel, nort, nsyme, nvb, nzeta, nzrvb, &
                        orbfr_is_unit, plc_const, strucopt
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nfxvbr, nzrvbr
logical(kind=iwp) :: changed, construc
logical(kind=iwp), external :: chpcmp_cvb ! ... Change of dimensioning variables ...

! Dimensioning for symmetry handling:
changed = .false.
if (chpcmp_cvb(nsyme)) changed = .true.
if (chpcmp_cvb(ndimrel)) changed = .true.
if (chpcmp_cvb(norbrel)) changed = .true.
if (chpcmp_cvb(nvb)) changed = .true.
if (chpcmp_cvb(nzrvb)) changed = .true.
if (chpcmp_cvb(nort)) changed = .true.
if (chpcmp_cvb(ndrot)) changed = .true.

orbfr_is_unit = ((ndimrel == 0) .and. (nfxorb == 0) .and. (nort == 0) .and. (.not. plc_const))
! Set ORBFR_IS_UNIT if optimization method is 'NONE':
if (imethod == 11) orbfr_is_unit = .true.
if (chpcmp_cvb(merge(1,0,orbfr_is_unit))) changed = .true.
nfxvbr = nfxvb
if (lfxvb == 1) nfxvbr = nvb-nfxvb
nzrvbr = nzrvb
if (lzrvb == 1) nzrvbr = nvb-nzrvb
! NPR arrays -- depend on nature of opt procedure
construc = ((nzrvbr > 0) .or. ((nfxvbr > 0) .and. (nfxvbr < nvb)) .or. (nzeta > 0))
if (construc) then
  if ((nvb > 20) .or. (.not. strucopt)) then
    iconstruc = 1
  else
    iconstruc = 2
  end if
else
  iconstruc = 0
end if
if (chpcmp_cvb(iconstruc)) changed = .true.
if (changed) call touch_cvb('MEM5')

return

end subroutine change5_cvb
