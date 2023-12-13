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
      Subroutine Scan_Inp_m(iRc)
* ------------------------------------------------------------
* Scan input lines after the '&MCPDFT' marker and until
* finding keyword 'END ' or the end of file.
* Keywords are identified according to file 'input_ras_mcpdft.fh'
* Logical flags in 'input_ras_mcpdft.fh' are set according to input.
* Return codes are _RC_ALL_IS_WELL_ or _RC_INPUT_ERROR_
* ------------------------------------------------------------

      use mcpdft_output, only: debug, lf, iPrLoc

      Implicit Real*8 (A-H,O-Z)

#include "rasdim.fh"
#include "warnings.h"
#include "WrkSpc.fh"
#include "input_ras_mcpdft.fh"
*
      Character*4 Command
      Character*180  Line
*

* If the return code is already set to indicate an error, there will
* be an error trace written out.
* Also at very high print level, there will be an error trace.
!        goto 200
      if(IPRLOC(1).ge.DEBUG  .or. iRc.ne._RC_ALL_IS_WELL_) goto 200

*--Find keywords in input and set keyword flags
      Do I=0,NKeys
       KeyFlags(I)=.False.
      End Do
      Rewind(LuInput)
  10  Continue
      Read(LuInput,'(A)',End=9910,Err=9920) Line
      Command=Line(1:4)
      Call UpCase(Command)
      Do iCmd=1,NKeys
        If ( Command.eq.Cmd(iCmd) ) Then
         KeyFlags(iCmd)=.TRUE.
* Special case: Skip title line.
         If ( Command.eq.'TITL') Then
          Read(LuInput,'(A)',End=9910,Err=9920) Line
         End If
* SVC (ugly hack) Special case: Skip fileorb line. FIXME: we need a more
* robust input scanning method so that these things are not necessary
         If ( Command.eq.'FILE') Then
          Read(LuInput,'(A)',End=9910,Err=9920) Line
         End If
         GoTo 20
        End If
      End Do
  20  Continue
      If ( .not.KeyEND ) GoTo 10
      Go To 9990

 200  continue
* Similar functionality, but with written trace:
      Do I=0,NKeys
       KeyFlags(I)=.False.
      End Do
      write(lf,*)' Scanning the input for keywords:'
      write(lf,*)' Rewinding LUInput=',LUInput
      Rewind(LuInput)
      write(lf,*)' OK after rewind.'
 210  Continue
      write(lf,*)' Reading a line...'
      Read(LuInput,'(A)',End=9910,Err=9920) Line
      write(lf,*)' '''//line(1:64)//' ...'''
      Command=Line(1:4)
      Call UpCase(Command)
      Do iCmd=1,NKeys
        If ( Command.eq.Cmd(iCmd) ) Then
          write(lf,*)' Understood keyword '''//Cmd(iCmd)//''''
          KeyFlags(iCmd)=.TRUE.
* Special case: Skip title line.
          If ( Command.eq.'TITL') Then
           write(lf,*)' Dummy read title line.'
           Read(LuInput,'(A)',End=9910,Err=9920) Line
          End If
          GoTo 220
        End If
      End Do
 220  Continue
      If ( .not.KeyEND ) GoTo 210
      Go To 9990

* Error exits ---------------------------------------
9910  CONTINUE
      write(lf,*)' Tried to read a new line. Hit End of record.'
      write(lf,*)' Last word was ',Command
      irc=_RC_INPUT_ERROR_
      GOTO 9990
*----------------------------------------------------
9920  CONTINUE
      write(lf,*)' Tried, and failed, to read a new line.'
      write(lf,*)' Last word was ',Command
      irc=_RC_INPUT_ERROR_
      GOTO 9990
*----------------------------------------------------
9990  CONTINUE
      Return
      End
