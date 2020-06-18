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
* Copyright (C) 1998, Roland Lindh                                     *
************************************************************************
      Subroutine Term_Ints(Verbose,Free_K2)
************************************************************************
*                                                                      *
*     Object: to deallocate memory in association with two-electron    *
*             calculations.                                            *
*                                                                      *
*     Author: Roland Lindh, Chemical Physics, University of Lund,      *
*             Sweden. January '98.                                     *
************************************************************************
      use k2_arrays, only: FT, Mem_DBLE, Mem_INT, Aux, iSOSym
      use Index_arrays
      Implicit Real*8 (A-H,O-Z)
*
#include "itmax.fh"
#include "info.fh"
#include "setup.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "nsd.fh"
#include "shinf.fh"
#include "status.fh"
      Logical Verbose, Free_K2
*
      If (ERI_Status.eq.Inactive) Return
      ERI_Status=Inactive
*     Call QEnter('T_I')
*
*     In case of semi-direct mode the memory is released externally.
*
      Call RlsMem_Ints
*
      If (Allocated(FT)) Call mma_deallocate(FT)
*
      If (Allocated(Mem_INT)) Then
         Call mma_deallocate(Mem_INT)
         Call mma_deallocate(Mem_DBLE)
         Call mma_deallocate(Aux)
      End If
*
      Call mma_deallocate(iSOSym)
*
*     Generate statistic of partioning
*
      If (Verbose) Call StatP(1)
*
      If (Indexation_Status.eq.Active) Then
         Indexation_Status=Inactive
         Call GetMem('nShBF','Free','Inte',ipShBF,mSkal*nIrrep)
         Call GetMem('ShLwC','Free','Inte',ipShLC,mSkal*nIrrep)
         Call GetMem('ShPSh','Free','Inte',ipShSh,mSkal*nIrrep)
         Call mma_deallocate(iSO2Sh)
         Call mma_deallocate(iCntr)
      End If
*
*---- Free memory for K2 data
      If (Free_K2) Call FreeK2
*
*     Call QExit('T_I')
      Return
      End
