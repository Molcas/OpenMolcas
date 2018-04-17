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
      Subroutine SetPos(LUnit,KeyIn,Line,iRc)
      Implicit Real*8 (A-H,O-Z)
      Character*(*) KeyIn
      Character*(*) Line
      Character*16 Command
      Character*16 Key
#include "warnings.fh"
#include "rasdim.fh"
#include "input_ras.fh"
#include "output_ras.fh"
      Parameter(ROUTINE='SetPos  ')
      Intrinsic len, min

      Call QENTER(ROUTINE)

* Read until, and including, a line beginning with a particular
* string in an ASCII file, assumed already opened, with unit
* number LUnit. That line is returned.
* Key lengths up to 16 bytes can be used, it is determined by
* the size of the input variable.
*
      IPRLEV=IPRLOC(1)
      iRc=_RC_ALL_IS_WELL_
      KLen=MIN(16,LEN(KeyIn))
      Key=' '
      Command=' '
      Rewind(LUnit)

      Key(1:KLen)=KeyIn(1:KLen)
      call upcase(Key)
10    Continue
      Read(LUnit,'(A)',End=9910,Err=9920) Line
      Command(1:KLen)=Line(1:KLen)
      call upcase(Command)
      If (Command.ne.Key) GoTo 10
      Return

*---  Error exits ----------------------
9910  CONTINUE
      If(IPRLEV.ge.TERSE) Then
       write(6,*)' SETPOS: Attempt to find an input line beginning'
       write(6,*)' with the keyword ''',KeyIn,''' failed.'
      End If
*      Call Quit(_RC_INPUT_ERROR_)
      iRc=_RC_INPUT_ERROR_
      Return
9920  CONTINUE
      If(IPRLEV.ge.TERSE) Then
       write(6,*)' SETPOS: Attempt to find an input line beginning'
       write(6,*)' with the keyword ''',KeyIn,''' failed.'
      End If
*      Call Quit(_RC_INPUT_ERROR_)
      iRc=_RC_INPUT_ERROR_
      Return
      End
