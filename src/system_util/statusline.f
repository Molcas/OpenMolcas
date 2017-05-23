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
      Subroutine StatusLine(STR,STR1)
************************************************************
*
*   <DOC>
*     <Name>StatusLine</Name>
*     <Syntax>Call StatusLine(STR,STR1)</Syntax>
*     <Arguments>
*       \Argument{STR}{Status}{Character*(*)}{in}
*       \Argument{STR1}{Status}{Character*(*)}{in}
*     </Arguments>
*     <Purpose>Update status file</Purpose>
*     <Dependencies></Dependencies>
*     <Author>V. Veryazov</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*      Print parameters into status file. There is no
*      difference between params.
*      Call StatusLine(ModuleName,': just started')
*     </Description>
*    </DOC>
*
************************************************************
      character*(*) STR, STR1
      Integer Lu
#ifdef _MOLCAS_MPP_
      Logical King
      External King
      if(.Not.King()) return
#endif
      Lu=2
      call molcas_open(Lu,'status')
      write(Lu,'(A,A)') STR,STR1
      close(Lu)
      Return
      End
