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
      program main
      implicit real*8 (a-h,o-z)
      Character*20 Module_Name
      Parameter (Module_Name = 'scf')

      Call Start(Module_Name)
c      Call xml_Open('module',' ',' ','scf')
      Call scf(ireturn)
c      Call xml_Close('module')
      Call Finish(ireturn)

      end
