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
      Subroutine Scan_Inp(iRc)
      Implicit Real*8 (A-H,O-Z)
* ------------------------------------------------------------
* Scan input lines after the '&RASSCF' marker and until
* finding keyword 'END ' or the end of file.
* Keywords are identified according to file 'input_ras.fh'
* Logical flags in 'input_ras.fh' are set according to input.
* Return codes are _RC_ALL_IS_WELL_ or _RC_INPUT_ERROR_
* ------------------------------------------------------------
#include "rasdim.fh"
#include "warnings.fh"
#include "WrkSpc.fh"
#include "input_ras.fh"
#include "output_ras.fh"
      Parameter(ROUTINE='Scan_Inp')
*
      Character*4 Command
      Character*180  Line
      Character*180 Get_LN
      External Get_LN

#ifdef _DMRG_
      External Get_ProgName
      Character*100 Get_ProgName
      Character*100 ProgName
      logical qcmaquis_input
#endif
*
      Call qEnter('Scan_Inp')

#ifdef _DMRG_
      qcmaquis_input = .false.
      ProgName       = Get_ProgName()
#endif

* If the return code is already set to indicate an error, there will
* be an error trace written out.
* Also at very high print level, there will be an error trace.
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

#ifdef _DMRG_
        if(ProgName(1:5).eq.'rassc' .and. command.eq.'ENDR')then
          qcmaquis_input = .false.
        else if(ProgName(1:5).eq.'dmrgs' .and. command.eq.'ENDO')then
          goto 9990
        end if
#endif

        If ( Command.eq.Cmd(iCmd) ) Then

#ifdef _DMRG_
          !> check for QCMaquis input section
          if(ProgName(1:5).eq.'rassc' .and. command.eq.'RGIN')then
            qcmaquis_input = .true.
            KeyFlags(iCmd) = .true.
          end if
          if(qcmaquis_input) goto 20
#endif

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
      write(6,*)' Scanning the input for keywords:'
      write(6,*)' Rewinding LUInput=',LUInput
      Rewind(LuInput)
      write(6,*)' OK after rewind.'
 210  Continue
      write(6,*)' Reading a line...'
      Read(LuInput,'(A)',End=9910,Err=9920) Line
      write(6,*)' '''//line(1:64)//' ...'''
      Command=Line(1:4)
      Call UpCase(Command)
      Do iCmd=1,NKeys

#ifdef _DMRG_
        if(ProgName(1:5).eq.'rassc' .and. command.eq.'ENDR')then
          qcmaquis_input = .false.
        else if(ProgName(1:5).eq.'dmrgs' .and. command.eq.'ENDO')then
          goto 9990
        end if
#endif

        If ( Command.eq.Cmd(iCmd) ) Then

#ifdef _DMRG_
          !> check for QCMaquis input section
          if(ProgName(1:5).eq.'rassc' .and. command.eq.'RGIN')then
            qcmaquis_input = .true.
            KeyFlags(iCmd) = .true.
          end if
          if(qcmaquis_input) goto 220
#endif

          write(6,*)' Understood keyword '''//Cmd(iCmd)//''''
          KeyFlags(iCmd)=.TRUE.
* Special case: Skip title line.
          If ( Command.eq.'TITL') Then
           write(6,*)' Dummy read title line.'
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
      write(6,*)' Tried to read a new line. Hit End of record.'
      write(6,*)' Last word was ',Command
      irc=_RC_INPUT_ERROR_
      GOTO 9990
*----------------------------------------------------
9920  CONTINUE
      write(6,*)' Tried, and failed, to read a new line.'
      write(6,*)' Last word was ',Command
      irc=_RC_INPUT_ERROR_
      GOTO 9990
*----------------------------------------------------
9990  CONTINUE
      Call qExit('Scan_Inp')
      Return
      End
