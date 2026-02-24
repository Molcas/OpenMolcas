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

subroutine cpinp(LUnit,iRc)

use UnixInfo, only: ProgName
use spool, only: Close_LuSpool, Disable_Spool, SpoolInp
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: LUnit, iRC
integer(kind=iwp) :: istatus, LUSpool
character(len=180) :: line
character :: ch
#ifdef _DMRG_
character(len=180) :: line2
#endif
integer(kind=iwp), external :: IsFreeUnit
#include "warnings.h"

iRc = _RC_ALL_IS_WELL_
! The following code will open, and return the unit number LUSpool,
! of an ASCII file with a copy of the presently used input.
! The records are strings of 180 characters, conforming to the
! standards set e.g. in src/util/inputil.f. This may change in the
! future, so look out!

call SpoolInp(LUSpool)
call Disable_Spool()
rewind(LUSpool)
! Now open a new file, and copy only the input between the '&RASSCF'
! and the 'End of Input' markers, inclusive, and skipping any commented
! lines. (The latter is put there by sbin/auto.plx, so it is safe to
! assume it is not abbreviated). The copied lines are left adjusted.
! Positioning LUSpool after the '&RASSCF' marker.
if (ProgName(1:5) == 'dmrgs') then
  call RdNLst(LuSpool,'DMRGSCF')
  call setpos(luspool,'OOPT',line,irc)
else
  call RdNLst(LuSpool,'RASSCF')
end if
! Opening a new file:
LUnit = 99
LUnit = IsFreeUnit(LUnit)
call Molcas_Open(LUnit,'CleanInput')
! Copy only the relevant lines of input:
line = ' '
line(1:7) = '&RASSCF'
write(LUnit,'(A180)') line
do
  read(luspool,'(A180)',iostat=istatus) line
  if (istatus /= 0) then
    ! Something went wrong...Let the caller handle it:
    iRc = _RC_INPUT_ERROR_
    return
  end if
  line = adjustl(line)
# ifdef _DMRG_
  if (ProgName(1:5) == 'dmrgs') then
    line2 = line
    call upcase(line2(1:4))
    if (line2(1:4) == 'ENDO') then
      line = ' '
      line(1:12) = 'End of Input'
      write(LUnit,'(A180)') line
      exit
    end if
  end if
# endif
  ch = line(1:1)
  if ((ch /= ' ') .and. (ch /= '*') .and. (ch /= '!')) write(LUnit,'(A180)') line
  call upcase(line(1:12))
  if (line(1:12) == 'END OF INPUT') exit
end do
call close_luspool(LUSpool)

return

end subroutine cpinp
