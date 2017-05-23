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
*
*  This function return 0 if module is running outside structure
*                       1 if module is running inside structure
*
      Function IsStructure ()
      character*256 Value
      Character*100 Get_SuperName
      External Get_SuperName
      Value=' '
      call getenvf('MOLCAS_STRUCTURE',Value)
      IsStructure=0
      if (Value.eq.'1')  IsStructure=1
      if (Get_SuperName().eq.'last_energy') IsStructure=1
      Return
      end
