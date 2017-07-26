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
*  SetMem
*
*> @brief
*>   Initialize and change status for memory control (common / MOLCAS_GetMem / MemCtl)
*>
*> @details
*> String is a string of any size and is not case sensitive.
*> The string contains keyword and the status of keyword in a form:
*> ``KEYWORD=STATUS``. The ``STATUS`` of the keyword can be ``ON`` or ``OFF``.
*> The ::SetMem subroutine can recognize only keywords:
*>
*> - ``TRACE``:  traces memory.
*> - ``SYSOUT``: unit of file which will be used as output for all kind of prints connected with the memory control.
*> - ``CLEAR``:  sets memory block for extra checking for memory allocation.
*>               A  ``_GARBLE_`` preprocessor option can be defined during compilation process to add an additional checking.
*> - ``QUERY``:  prints status of the Molcas_query.
*> - ``CHECK``:  check the internal state of the MA.
*>
*> @param[in] String ``TRACE=ON`` / ``TRACE=OFF`` / ``SYSOUT=ON`` / ``SYSOUT=OFF`` / ``CLEAR=ON`` / ``CLEAR=OFF`` / ``QUERY=ON`` / ``QUERY=OFF`` / ``CHECK=ON`` / ``CHECK=OFF``
************************************************************************
      Subroutine SetMem (String)
*

#include "SysCtl.fh"
#include "mama.fh"
*
*
      Character*(*) String
      Character*20  Token
*----------------------------------------------------------------------*
*     Initialize the Common / MemCtl / the first time it is referenced *
*----------------------------------------------------------------------*
      If ( MemCtl(ipStat).ne.ON ) then
         Call IniMem
      End if
*----------------------------------------------------------------------*
*     read default parameters from Common / MemCtl /                   *
*----------------------------------------------------------------------*
      iW=MemCtl(ipSysOut)
      If ( MemCtl(ipTrace).eq.ON ) then
         Write(iW,*) ' <<< Entering SetMem >>>'
      End If
*----------------------------------------------------------------------*
*     extract the first token and convert it into standard format      *
*----------------------------------------------------------------------*
      Call StdFmt(String,Token)
      If ( Token.eq.' ' ) Return
*----------------------------------------------------------------------*
*     replace default values                                           *
*----------------------------------------------------------------------*
      lToken=LEN(Token)
      If ( Token(1:6).eq.'TRACE=') then
         If ( Token(7:8).eq.'ON' ) then
            MemCtl(ipTrace) = ON
            Return
         Else If ( Token(7:9).eq.'OFF' ) then
            MemCtl(ipTrace) = OFF
            Return
         End If
      Else if ( Token(1:7).eq.'SYSOUT=' ) then
         Read(Token(8:lToken),*) MemCtl(ipSysOut)
         Return
      Else if ( Token(1:6).eq.'CLEAR=') then
         If ( Token(7:8).eq.'ON' ) then
            MemCtl(ipClear) = ON
            Return
         Else If ( Token(7:9).eq.'OFF' ) then
            MemCtl(ipClear) = OFF
            Return
         End If
      Else if ( Token(1:6).eq.'QUERY=') then
         If ( Token(7:8).eq.'ON' ) then
            MemCtl(ipQuery) = ON
            Return
         Else If ( Token(7:9).eq.'OFF' ) then
            MemCtl(ipQuery) = OFF
            Return
         End If
      Else if ( Token(1:6).eq.'CHECK=') then
         If ( Token(7:8).eq.'ON' ) then
            MemCtl(ipCheck) = ON
            Return
         Else If ( Token(7:9).eq.'OFF' ) then
            MemCtl(ipCheck) = OFF
            Return
         End If
      Else
         Write (6,*) 'SetMem: illegal option'
         Write (6,'(2A)') 'Option:',Token
         Call QTrace()
         Call Abend()
      End if
*----------------------------------------------------------------------*
*     exit                                                             *
*----------------------------------------------------------------------*
      If ( MemCtl(ipTrace).eq.ON ) then
         Write(iW,*) ' <<< Exiting SetMem >>>'
      End If
      Return
      End
