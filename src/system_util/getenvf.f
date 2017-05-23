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

*
*  This is a simple wrapper for getenv
*   Note:
*     Value must be defined as character*256
*     Name must be terminated by space
*
*  Sample:
*    Character Env*40
*    Character Value*256
*      Env='MOLCAS '
*      call getenvf(Env,Value)
*      print *, 'MOLCAS=',Value
*

      Subroutine getenvf(Name,Value)
      Character*(*) Name, Value
      value=' '
      ilen=len(Name)
      maxlen=len(Value)
      Call getenvf2c(Name,ilen,Value,maxlen,irl)
       if(irl.eq.0) then
        Value=' '
       else
        Value=Value(1:irl)
       endif
      Return
      End
