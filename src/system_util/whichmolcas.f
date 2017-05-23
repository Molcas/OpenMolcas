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
* Copyright (C) 2000-2016, Valera Veryazov                             *
************************************************************************
******************************************************************************
*                                                                            *
* Author:   Valera Veryazov 2000-2016                                        *
*           Theoretical Chemistry                                            *
*           Lund University                                                  *
*           Sweden                                                           *
*                                                                            *
******************************************************************************

      Subroutine WhichMolcas(Molcas)
c
c  if MOLCAS_STAMP='Alone'
c     get value of MOLCAS.
c
c
      Character*(*) Molcas
      Character Env*40
      Env='MOLCAS_STAMP'
*
      Molcas=' '
      Call getenvf(Env,Molcas)
      if (Molcas(1:1) .ne. 'A') then
      Molcas=' '
      Return
      endif
      Env='MOLCAS'
*
      Molcas=' '
      Call getenvf(Env,Molcas)
      Return
      End
