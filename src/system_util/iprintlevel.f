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
      Integer function iPrintLevel(Level)
************************************************************
*
*   <DOC>
*     <Name>iPrintLevel</Name>
*     <Syntax>iPrintLevel(Level)</Syntax>
*     <Arguments>
*       \Argument{Level}{Print Level}{Integer}{in}
*     </Arguments>
*     <Purpose>Check or set global print level</Purpose>
*     <Dependencies>getenvf</Dependencies>
*     <Author>V. Veryazov</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*       If Level is a number >= 0 - set print Level
*       else return current print level
*
*       In a first call of the routine environment variable
*          MOLCAS\_PRINT is used to set the initial print level.
*
*       Allowed values are: SILENT (0), TERSE (1), NORMAL (2),
*          VERBOSE (3), DEBUG (4) and INSANE (5)
*
*     </Description>
*    </DOC>
*
************************************************************
      Data isFirst/0/
      Save isFirst
      common /nPrintLevel/ nPrintLevel
      character *80 Name, Value
      if(Level.ge.0) Then
        nPrintLevel=Level
        isFirst=1
        iPrintLevel=Level
        return
      endif
      if (isFirst.eq.0) then
        Name='MOLCAS_PRINT'
        call getenvf(Name,Value)
            call UpCase(Value)
               if     (Value.eq.'SILENT'.or.Value.eq.'0') then
                       nPrintLevel=0
               elseif (Value.eq.'TERSE'.or.Value.eq.'1') then
                       nPrintLevel=1
               elseif (Value.eq.'NORMAL'.or.Value.eq.'2') then
                       nPrintLevel=2
               elseif (Value.eq.'VERBOSE'.or.Value.eq.'3') then
                       nPrintLevel=3
               elseif (Value.eq.'DEBUG'.or.Value.eq.'4') then
                       nPrintLevel=4
               elseif (Value.eq.'INSANE'.or.Value.eq.'5') then
                       nPrintLevel=5
               else
                       nPrintLevel=2
               endif
               iFirst=1
         endif
        iPrintLevel=nPrintLevel
      return
      end
