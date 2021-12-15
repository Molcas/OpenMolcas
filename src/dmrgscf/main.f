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
* Copyright (C) 2017, Stefan Knecht                                    *
************************************************************************

      program main
      implicit none
      integer :: ireturn
      character(len=20), parameter :: module_name = 'dmrgscf'

      call start(module_name)
      call dmrgscf(ireturn)
      call finish(ireturn)

      end program main
