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
      Subroutine CollapseOutput(iSw, STR)
************************************************************
*
*   <DOC>
*     <Name>CollapseOutput</Name>
*     <Syntax>Call CollapseOutput(iSw,STR)</Syntax>
*     <Arguments>
*       \Argument{iSw}{Begin/End switch}{Integer}{in}
*       \Argument{STR}{Message}{Character*(*)}{in}
*     </Arguments>
*     <Purpose>Place markers for collapsed parts in output</Purpose>
*     <Dependencies></Dependencies>
*     <Author>V. Veryazov</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*         if iSw = 1, print Message preceded by a marker to
*          start collapsed part in the output.
*         if iSw = 0, place a marker to terminate collapsed part.
*          STR is not used, but recommended.
*     </Description>
*    </DOC>
*
************************************************************
      character*(*) STR
      Integer iSw
      Common /icolorize/icolorize
      if(icolorize.eq.1) then
       if(iSw.eq.1) then
         write(6,'(A,A)') '++ ',STR(:mylen(STR))
       else
        write(6,'(A)') '--'
       endif
      else
         if(iSw.eq.1) write(6,'(A)') STR(:mylen(STR))
      endif
      Return
      End
