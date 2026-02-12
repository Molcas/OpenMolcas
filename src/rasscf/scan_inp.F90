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

subroutine Scan_Inp(iRc)
! ------------------------------------------------------------
! Scan input lines after the '&RASSCF' marker and until
! finding keyword 'END ' or the end of file.
! Keywords are identified according to file 'input_ras.F90'
! Logical flags in 'input_ras.F90' are set according to input.
! Return codes are _RC_ALL_IS_WELL_ or _RC_INPUT_ERROR_
! ------------------------------------------------------------

#ifdef _DMRG_
use UnixInfo, only: ProgName
#endif
use input_ras, only: nKeys, CMD, KeyFlags, KeyEND, LuInput
use PrintLevel, only: DEBUG
use output_ras, only: IPRLOC

implicit none
integer iRC
integer I, iCMD
character(len=4) Command
character(len=180) Line
#ifdef _DMRG_
logical qcmaquis_input
#endif
#include "warnings.h"

#ifdef _DMRG_
qcmaquis_input = .false.
#endif

! If the return code is already set to indicate an error, there will
! be an error trace written out.
! Also at very high print level, there will be an error trace.
if ((IPRLOC(1) >= DEBUG) .or. (iRc /= _RC_ALL_IS_WELL_)) goto 200

! Find keywords in input and set keyword flags
do I=0,NKeys
  KeyFlags(I) = .false.
end do
rewind(LuInput)
10 continue
read(LuInput,'(A)',end=9910,err=9920) Line
Command = Line(1:4)
call UpCase(Command)
do iCmd=1,NKeys

# ifdef _DMRG_
  if ((ProgName(1:5) == 'rassc') .and. (command == 'ENDR')) then
    qcmaquis_input = .false.
  else if ((ProgName(1:5) == 'dmrgs') .and. (command == 'ENDO')) then
    goto 9990
  end if
# endif

  if (Command == Cmd(iCmd)) then

#   ifdef _DMRG_
    !> check for QCMaquis input section
    if ((ProgName(1:5) == 'rassc') .and. (command == 'RGIN')) then
      qcmaquis_input = .true.
      KeyFlags(iCmd) = .true.
    end if
    if (qcmaquis_input) goto 20
#   endif

    KeyFlags(iCmd) = .true.
    ! Special case: Skip title line.
    if (Command == 'TITL') read(LuInput,'(A)',end=9910,err=9920) Line
    ! SVC (ugly hack) Special case: Skip fileorb line. FIXME: we need a more
    ! robust input scanning method so that these things are not necessary
    if (Command == 'FILE') read(LuInput,'(A)',end=9910,err=9920) Line
    goto 20

  end if

end do
20 continue
if (.not. KeyEND) goto 10
Go To 9990

200 continue
! Similar functionality, but with written trace:
do I=0,NKeys
  KeyFlags(I) = .false.
end do
write(6,*) ' Scanning the input for keywords:'
write(6,*) ' Rewinding LUInput=',LUInput
rewind(LuInput)
write(6,*) ' OK after rewind.'
210 continue
write(6,*) ' Reading a line...'
read(LuInput,'(A)',end=9910,err=9920) Line
write(6,*) ' '''//line(1:64)//' ...'''
Command = Line(1:4)
call UpCase(Command)
do iCmd=1,NKeys

# ifdef _DMRG_
  if ((ProgName(1:5) == 'rassc') .and. (command == 'ENDR')) then
    qcmaquis_input = .false.
  else if ((ProgName(1:5) == 'dmrgs') .and. (command == 'ENDO')) then
    goto 9990
  end if
# endif

  if (Command == Cmd(iCmd)) then

#   ifdef _DMRG_
    !> check for QCMaquis input section
    if ((ProgName(1:5) == 'rassc') .and. (command == 'RGIN')) then
      qcmaquis_input = .true.
      KeyFlags(iCmd) = .true.
    end if
    if (qcmaquis_input) goto 220
#   endif

    write(6,*) ' Understood keyword '''//Cmd(iCmd)//''''
    KeyFlags(iCmd) = .true.
    ! Special case: Skip title line.
    if (Command == 'TITL') then
      write(6,*) ' Dummy read title line.'
      read(LuInput,'(A)',end=9910,err=9920) Line
    end if
    goto 220
  end if

end do
220 continue
if (.not. KeyEND) goto 210
Go To 9990

! Error exits ---------------------------------------
9910 continue
write(6,*) ' Tried to read a new line. Hit End of record.'
write(6,*) ' Last word was ',Command
irc = _RC_INPUT_ERROR_
goto 9990
!----------------------------------------------------
9920 continue
write(6,*) ' Tried, and failed, to read a new line.'
write(6,*) ' Last word was ',Command
irc = _RC_INPUT_ERROR_
goto 9990
!----------------------------------------------------
9990 continue

end subroutine Scan_Inp
