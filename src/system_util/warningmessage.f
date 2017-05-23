************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Valera Veryazov                                        *
************************************************************************
      Subroutine WarningMessage(Level,STR)
************************************************************
*
*   <DOC>
*     <Name>WarningMessage</Name>
*     <Syntax>Call WarningMessage(Level,STR)</Syntax>
*     <Arguments>
*       \Argument{Level}{Warning level}{Integer}{in}
*       \Argument{STR}{Message}{Character*(*)}{in}
*     </Arguments>
*     <Purpose>Print warning message</Purpose>
*     <Dependencies></Dependencies>
*     <Author>V. Veryazov</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects>Routine updates MaxWarnMess</Side_Effects>
*     <Description>
*       Print message in uniform format.
*       Recommended levels:
*       0 - for message without recording that it was an error
*       1 - for warnings
*       2 - for errors
*     </Description>
*    </DOC>
*
************************************************************
      character*(*) STR
      Integer Level
      common /WarnMess/ MaxWarnMess
      if(Level.gt.MaxWarnMess) MaxWarnMess=Level
      call SysPutsStart()
      call SysPuts(STR,' ',' ')
      call SysPutsEnd()

c      write(6,'(A)') '*** '
c      jj=1
c10    i=index(STR(jj:),';')
c      if(i.eq.0) then
c      goto 20
c      else
c      write(6,'(A,A)') '*** ',STR(jj:jj+i-2)
c      jj=i+jj
c      goto 10
c      endif
c20    write(6,'(A,A)') '*** ',STR(jj:)
c      write(6,'(A)') '*** '
      Return
      End
      Subroutine WarningInit
      common /WarnMess/ MaxWarnMess
c this is unfinished routine. it must be possible to keep more info
      MaxWarnMess=-1
      return
      end
      Subroutine WarningCheckOut(iWarn)
      common /WarnMess/ MaxWarnMess
      iWarn=MaxWarnMess
      return
      end
