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
      Logical Function RF_On()
      Implicit Real*8 (a-h,o-z)
*     Call Get_iOption(iOption)
      Call Get_iScalar('System BitSwitch',iOption)
*
      RF_On=Iand(iOption,4).Eq.4
      Return
      End
