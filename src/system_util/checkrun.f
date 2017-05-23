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
* Copyright (C) 2000,2016, Valera Veryazov                             *
************************************************************************
******************************************************************************
*                                                                            *
* Author:   Valera Veryazov 2000-2016                                        *
*           Theoretical Chemistry                                            *
*           Lund University                                                  *
*           Sweden                                                           *
*                                                                            *
******************************************************************************
      Subroutine CheckRun()
************************************************************
*
*   <DOC>
*     <Name>CheckRun</Name>
*     <Syntax>Call CheckRun()</Syntax>
*     <Arguments>
*     </Arguments>
*     <Purpose>Check that the program was started via molcas shell script</Purpose>
*     <Dependencies>Getenvf</Dependencies>
*     <Author>V. Veryazov</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*       The routine check MOLCAS\_STAMP and quit if this variable
*       has wrong value (to prevent the creation of dummy files)
*     </Description>
*    </DOC>
*
************************************************************
      Character Molcas*256
      Character Env*40
      Env='MOLCAS_STAMP'
*
      Molcas=' '
      Call getenvf(Env,Molcas)
      if (Molcas(1:1) .eq. 'A') then
      Return
      endif
      if (Molcas(1:1) .ne. '5') then
        write (6,*) 'Usage: molcas module_name input'
        Call Abend()
      endif
      Return
      End
