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
!  SetMem
!
!> @brief
!>   Initialize and change status for memory control (common / MOLCAS_GetMem / MemCtl)
!>
!> @details
!> String is a string of any size and is not case sensitive.
!> The string contains keyword and the status of keyword in a form:
!> ``KEYWORD=STATUS``. The ``STATUS`` of the keyword can be ``ON`` or ``OFF``.
!> The ::SetMem subroutine can recognize only keywords:
!>
!> - ``TRACE``:  traces memory.
!> - ``SYSOUT``: unit of file which will be used as output for all kind of prints connected with the memory control.
!> - ``CLEAR``:  sets memory block for extra checking for memory allocation.
!>               A  ``_GARBLE_`` preprocessor option can be defined during compilation process to add an additional checking.
!> - ``QUERY``:  prints status of the Molcas_query.
!> - ``CHECK``:  check the internal state of the MA.
!>
!> @param[in] String ``TRACE=ON`` / ``TRACE=OFF`` / ``SYSOUT=ON`` / ``SYSOUT=OFF`` / ``CLEAR=ON`` / ``CLEAR=OFF`` / ``QUERY=ON`` / ``QUERY=OFF`` / ``CHECK=ON`` / ``CHECK=OFF``
!***********************************************************************

subroutine SetMem(String)

use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: String
#include "SysCtl.fh"
#include "mama.fh"
integer(kind=iwp) :: iW, lToken
character(len=20) :: Token

!----------------------------------------------------------------------*
!     Initialize the Common / MemCtl / the first time it is referenced *
!----------------------------------------------------------------------*
if (MemCtl(ipStat) /= ON) call IniMem()
!----------------------------------------------------------------------*
!     read default parameters from Common / MemCtl /                   *
!----------------------------------------------------------------------*
iW = MemCtl(ipSysOut)
if (MemCtl(ipTrace) == ON) write(iW,*) ' <<< Entering SetMem >>>'
!----------------------------------------------------------------------*
!     extract the first token and convert it into standard format      *
!----------------------------------------------------------------------*
call StdFmt(String,Token)
if (Token == ' ') return
!----------------------------------------------------------------------*
!     replace default values                                           *
!----------------------------------------------------------------------*
lToken = len(Token)
if (Token(1:6) == 'TRACE=') then
  if (Token(7:8) == 'ON') then
    MemCtl(ipTrace) = ON
    return
  else if (Token(7:9) == 'OFF') then
    MemCtl(ipTrace) = OFF
    return
  end if
else if (Token(1:7) == 'SYSOUT=') then
  read(Token(8:lToken),*) MemCtl(ipSysOut)
  return
else if (Token(1:6) == 'CLEAR=') then
  if (Token(7:8) == 'ON') then
    MemCtl(ipClear) = ON
    return
  else if (Token(7:9) == 'OFF') then
    MemCtl(ipClear) = OFF
    return
  end if
else if (Token(1:6) == 'QUERY=') then
  if (Token(7:8) == 'ON') then
    MemCtl(ipQuery) = ON
    return
  else if (Token(7:9) == 'OFF') then
    MemCtl(ipQuery) = OFF
    return
  end if
else if (Token(1:6) == 'CHECK=') then
  if (Token(7:8) == 'ON') then
    MemCtl(ipCheck) = ON
    return
  else if (Token(7:9) == 'OFF') then
    MemCtl(ipCheck) = OFF
    return
  end if
else
  write(u6,*) 'SetMem: illegal option'
  write(u6,'(2A)') 'Option:',Token
  call Abend()
end if
!----------------------------------------------------------------------*
!     exit                                                             *
!----------------------------------------------------------------------*
if (MemCtl(ipTrace) == ON) write(iW,*) ' <<< Exiting SetMem >>>'

end subroutine SetMem
