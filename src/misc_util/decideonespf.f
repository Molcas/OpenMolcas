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
      Subroutine DecideonESPF(Do_ESPF)

      Implicit Real*8 (a-h,o-z)
      Logical Do_ESPF


      Call Get_iScalar('System BitSwitch',iOption)
*
*  this is a logical variable (logical function)
      Do_ESPF=Iand(iOption,2048).Eq.2048
*
      Return
      End
