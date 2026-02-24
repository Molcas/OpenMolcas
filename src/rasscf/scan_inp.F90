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
use input_ras, only: CMD, Key, KeyFlags, LuInput, nKeys
use PrintLevel, only: DEBUG
use output_ras, only: IPRLOC
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: iRC
integer(kind=iwp) :: iCMD, istatus
character(len=180) :: Line
character(len=4) :: Command
#ifdef _DMRG_
logical(kind=iwp) :: qcmaquis_input
#endif
#include "warnings.h"

#ifdef _DMRG_
qcmaquis_input = .false.
#endif

! If the return code is already set to indicate an error, there will
! be an error trace written out.
! Also at very high print level, there will be an error trace.
if ((IPRLOC(1) < DEBUG) .and. (iRc == _RC_ALL_IS_WELL_)) then

  ! Find keywords in input and set keyword flags
  KeyFlags(:) = .false.
  rewind(LuInput)
  outer1: do
    read(LuInput,'(A)',iostat=istatus) Line
    if (istatus /= 0) then
      call Error(istatus)
      return
    end if
    Command = Line(1:4)
    call UpCase(Command)
    do iCmd=1,NKeys

#     ifdef _DMRG_
      if ((ProgName(1:5) == 'rassc') .and. (command == 'ENDR')) then
        qcmaquis_input = .false.
      else if ((ProgName(1:5) == 'dmrgs') .and. (command == 'ENDO')) then
        exit outer1
      end if
#     endif

      if (Command == Cmd(iCmd)) then

#       ifdef _DMRG_
        !> check for QCMaquis input section
        if ((ProgName(1:5) == 'rassc') .and. (command == 'RGIN')) then
          qcmaquis_input = .true.
          KeyFlags(iCmd) = .true.
        end if
        if (qcmaquis_input) exit
#       endif

        KeyFlags(iCmd) = .true.
        ! Special case: Skip title line.
        if (Command == 'TITL') read(LuInput,'(A)',iostat=istatus) Line
        if (istatus /= 0) then
          call Error(istatus)
          return
        end if
        ! SVC (ugly hack) Special case: Skip fileorb line. FIXME: we need a more
        ! robust input scanning method so that these things are not necessary
        if (Command == 'FILE') read(LuInput,'(A)',iostat=istatus) Line
        if (istatus /= 0) then
          call Error(istatus)
          return
        end if
        exit

      end if

    end do
    if (Key('END')) exit outer1
  end do outer1

else

  ! Similar functionality, but with written trace:
  KeyFlags(:) = .false.
  write(u6,*) ' Scanning the input for keywords:'
  write(u6,*) ' Rewinding LUInput=',LUInput
  rewind(LuInput)
  write(u6,*) ' OK after rewind.'
  outer2: do
    write(u6,*) ' Reading a line...'
    read(LuInput,'(A)',iostat=istatus) Line
    if (istatus /= 0) then
      call Error(istatus)
      return
    end if
    write(u6,*) ' '''//line(1:64)//' ...'''
    Command = Line(1:4)
    call UpCase(Command)
    do iCmd=1,NKeys

#     ifdef _DMRG_
      if ((ProgName(1:5) == 'rassc') .and. (command == 'ENDR')) then
        qcmaquis_input = .false.
      else if ((ProgName(1:5) == 'dmrgs') .and. (command == 'ENDO')) then
        exit outer2
      end if
#     endif

      if (Command == Cmd(iCmd)) then

#       ifdef _DMRG_
        !> check for QCMaquis input section
        if ((ProgName(1:5) == 'rassc') .and. (command == 'RGIN')) then
          qcmaquis_input = .true.
          KeyFlags(iCmd) = .true.
        end if
        if (qcmaquis_input) exit
#       endif

        write(u6,*) ' Understood keyword '''//Cmd(iCmd)//''''
        KeyFlags(iCmd) = .true.
        ! Special case: Skip title line.
        if (Command == 'TITL') then
          write(u6,*) ' Dummy read title line.'
          read(LuInput,'(A)',iostat=istatus) Line
          if (istatus /= 0) then
            call Error(istatus)
            return
          end if
        end if
        exit
      end if

    end do
    if (Key('END')) exit outer2
  end do outer2

end if

contains

subroutine Error(code)

  integer(kind=iwp), intent(in) :: code

  if (code < 0) then
    write(u6,*) ' Tried to read a new line. Hit End of record.'
  else
    write(u6,*) ' Tried, and failed, to read a new line.'
  end if
  write(u6,*) ' Last word was ',Command
  irc = _RC_INPUT_ERROR_

end subroutine Error

end subroutine Scan_Inp
