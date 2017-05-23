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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      SubRoutine CloseR
************************************************************************
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             September '90                                            *
************************************************************************
      use vRys_RW
      Implicit Real*8 (A-H,O-Z)
#include "stdalloc.fh"
*
      If (.Not.Allocated(iHerW2)) Return
      Call mma_deallocate(iHerW2)
      Call mma_deallocate(iHerR2)
      Call mma_deallocate(HerW2)
      Call mma_deallocate(HerR2)
      Call mma_deallocate(Cff)
      Call mma_deallocate(x0)
      Call mma_deallocate(Map)
      call mma_deallocate(ddx)
      Call mma_deallocate(TMax)
*
      Return
      End
