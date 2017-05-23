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
      SubRoutine Free_HerRW()
      use Her_RW
      Implicit Real*8 (A-H,O-Z)
#include "stdalloc.fh"
      If (Allocated(iHerR)) Call mma_deallocate(iHerR)
      If (Allocated(iHerW)) Call mma_deallocate(iHerW)
      If (Allocated( HerR)) Call mma_deallocate( HerR)
      If (Allocated( HerW)) Call mma_deallocate( HerW)
      Return
      End
