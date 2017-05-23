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
      subroutine rlsmem_ints
      implicit real*8 (a-h,o-z)
#include "setup.fh"
c
      If(MemMax_int.gt.0) then
        Call GetMem('MaxMem','Free','Real',ipMem_int,MemMax_int)
        MemMax_int=0
        ipMem_int=0
      End if
      return
      end
