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
* Copyright (C) 1992, Roland Lindh                                     *
************************************************************************
      SubRoutine ClsSew(iFrom)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
************************************************************************
      use Her_RW
      use Real_Spherical
      use EFP_module
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "timtra.fh"
#include "rctfld.fh"
#include "nsd.fh"
#include "setup.fh"
#include "status.fh"
*
      if(iFrom.eq.0) then
      iu=16
      iu=isfreeunit(iu)
      open(iu,file="ORB.std")
      Call Koor2file(iu)
      Call Basi2file(iu)
      close(iu)
      endif
      If (Seward_Status.eq.InActive) Return
*
      Call Term_Ints(.False.,.True.)
      Call Free_RctFld(iXPolType)
      Call Free_HerRW()
*
      If (Allocated(RSph)) Call mma_deallocate(RSph)
      If (Allocated(ipSph)) Call mma_deallocate(ipSph)
      If (Info_Status.eq.Active) Then
         Call GetMem('Info','Free','Real',LctInf,nInfo)
         Info_Status=InActive
      End If
*
      Call Free_iSD()
      Call Freek2()
      Call CloseR()
*
      If (EFP) Then
         Deallocate(FRAG_TYPE)
         Deallocate(ABC)
         Deallocate(EFP_COORS)
         EFP=.FALSE.
      End If
*
      Seward_Status=InActive
      Return
      End
