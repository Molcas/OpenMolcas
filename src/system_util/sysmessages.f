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
* Copyright (C) 2001, Valera Veryazov                                  *
************************************************************************
************************************************************************
*                                                                      *
*     purpose:                                                         *
*       replace SYSDB                                                  *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     V.Veryazov University of Lund, 2001                              *
*                                                                      *
************************************************************************
*  SysHalt
*
*> @brief
*>   Quit calculation
*> @author V. Veryazov
*>
*> @details
*> A routine to stop calculation because of internal error or
*> inconsistency in the code.
*>
*> @param[in] Location routine name
************************************************************************
      Subroutine SysHalt(Location)
      Character *(*) Location
       Call SysAbendMsg(Location,'Internal error',' ')
      return
      end
************************************************************************
*  SysWarnMsg
*
*> @brief
*>   Print nice formatted warning message
*> @author V. Veryazov
*>
*> @details
*> Print formatted message.
*> For a set of standard messages (started from ``MSG:``) the aliases
*> (defined in ::SysExpand) will be used.
*>
*> @param[in] Location routine name
*> @param[in] Text1    message text
*> @param[in] Text2    message text
************************************************************************
      Subroutine SysWarnMsg(Location,Text1,Text2)
      Character*(*) Location,Text1, Text2
      Character Str*256
      common /WarnMess/ MaxWarnMess
      Level=1
c for these messages assume that Level is 1
      if(Level.gt.MaxWarnMess) MaxWarnMess=Level
      call SysPutsStart()
      call SysPuts('Location: ',Location,'\n\n\n')
      call SysExpand(Text1,Str,i)
      if(i.eq.0) then
        call SysPuts(Text1,' ',Text2)
      else
        call SysPuts(str(:i),' ',Text2)
      endif
      call SysPutsEnd()

      Return
      End
************************************************************************
      Subroutine SysWarnFileMsg(Location,TheFile,Text1,Text2)
      Character*(*) Location,Text1, Text2, TheFile
      Character Str*256
      call SysPutsStart()
      call SysPuts('Location: ',Location,'\n')
      call SysExpand(TheFile,Str,i)
      call SysPuts('File: ',TheFile,'\n\n\n')
      call SysExpand(Text1,Str,i)
      if(i.eq.0) then
        call SysPuts(Text1,' ',Text2)
      else
        call SysPuts(str(:i),' ',Text2)
      endif
      call SysPutsEnd()

      Return
      End
************************************************************************
*  SysAbendMsg
*
*> @brief
*>   Stop calculation
*> @author V. Veryazov
*>
*> @details
*> Print formatted message and stop the calculation.
*>
*> @param[in] Location routine name
*> @param[in] Text1    message text
*> @param[in] Text2    message text
************************************************************************
      Subroutine SysAbendMsg(Location,Text1,Text2)
      Character*(*) Location,Text1, Text2
      call SysWarnMsg(Location,Text1,Text2)
      Call Abend()
      Return
      End
************************************************************************
      Subroutine SysAbendFileMsg(Location,TheFile,Text1,Text2)
      Character*(*) Location,Text1, Text2, TheFile
      call SysWarnFileMsg(Location,TheFile,Text1,Text2)
      Call Abend()
      Return
      End
************************************************************************
      Subroutine SysQuitFileMsg(rc,Location,TheFile,Text1,Text2)
      Integer rc
      Character*(*) Location,Text1, Text2, TheFile
      call SysWarnFileMsg(Location,TheFile,Text1,Text2)
      Call Quit(rc)
      Return
      End
************************************************************************
      Subroutine SysQuitMsg(rc,Location,Text1,Text2)
      Integer rc
      Character*(*) Location,Text1, Text2
      call SysWarnMsg(Location,Text1,Text2)
      Call Quit(rc)
      Return
      End
************************************************************************
      Subroutine SysCondMsg(Text1,N1,Text2,N2)
      Character s*64
      Character*(*) Text1, Text2
      call SysPuts('Condition: ',Text1,' ')
      write(s,'(i16,a,i16)') N1,Text2,N2
      call SysPuts('Actual   : ',s,' ')
      call SysPutsEnd()
      call Abend()
      return
      end
************************************************************************
      Subroutine SysValueMsg(Text1,N1)
      Character*(*) Text1
      call SysValueWarnMsg(Text1,N1)
      call SysPutsEnd()
      call Abend()
      return
      end
************************************************************************
      Subroutine SysValueWarnMsg(Text1,N1)
      Character s*20
      Character*(*) Text1
      write(s,'(a,i16)') ' = ',N1
      call SysPuts('Value: ',Text1, s)
      return
      end
************************************************************************
      Subroutine SysPutsStart()
      character c
      c='#'
      write (6,'(a,79a1)') ' ',(c,i=1,79)
      write (6,'(a,79a1)') ' ',(c,i=1,79)
      write (6,'(a,73x,a)') ' ###','###'
      write (6,'(a,73x,a)') ' ###','###'
      return
      end
************************************************************************
      Subroutine SysPutsEnd()
      character c
      c='#'
      write (6,'(a,73x,a)') ' ###','###'
      write (6,'(a,73x,a)') ' ###','###'
      write (6,'(a,79a1)') ' ',(c,i=1,79)
      write (6,'(a,79a1)') ' ',(c,i=1,79)
      return
      end
************************************************************************
      Subroutine SysFileMsg(Location,Text1,Lu,Text2)
      Character*(*) Location,Text1, Text2
      character *256 str
      call SysPutsStart()
      call SysPuts('Location: ',Location,'\n')
      write(str,*) Lu
      call SysPuts('Unit    : ',str,' ')
      str=' '
      inquire (unit=lu, name=str)
       if(str.ne.' ') then
      call SysPuts('File    : ',str,'\n')
       endif

      call SysExpand(Text1,Str,i)
      if(i.eq.0) then
        call SysPuts(Text1,'\n',Text2)
      else
        call SysPuts(str(:i),'\n',Text2)
      endif
      call SysPutsEnd()
      call Abend()
      return
      end
