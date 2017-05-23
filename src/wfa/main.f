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
* Copyright (C) 2017, Felix Plasser                                    *
************************************************************************

      Program main
      Implicit None
      Integer ireturn
      Character(20), Parameter :: module_name = 'wfa'

      Call start(module_name)
      Call wfa(ireturn)
      Call finish(ireturn)

      End Program main
