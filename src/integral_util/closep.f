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
      SubRoutine CloseP
************************************************************************
*                                                                      *
* Object: to close the handling of the 2nd order density matrix.       *
*                                                                      *
* Called from: Seward                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
************************************************************************
      use aces_stuff
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "real.fh"
#include "pso.fh"
#include "setup.fh"
#include "mp2alaska.fh"
      Logical DoCholesky
*
      iRout = 249
      iPrint = nPrint(iRout)
*     Call qEnter('CloseP')
*
      If(case_mp2) then
         Call DecideOnCholesky(DoCholesky)
         If(.not. DoCholesky) Then
            Call DaClos(LuGam)
         End If
      End If
      If (Gamma_On) Then
************ columbus interface ****************************************
         Call DaClos(LuGamma)
         Call mma_deallocate(Bin)
         Call mma_deallocate(G_Toc)
         Call mma_deallocate(SO2cI)
      End If
*
      If (lPSO) Then
         Call GetMem(' G2 ','Free','Real',ipG2,nG2)
         Call GetMem(' G1 ','Free','Real',ipG1,nG1)
      End If
      Call GetMem('CMO  ','Free','Real',ipCMO  ,mCMO)
      Call GetMem('DSVar','Free','Real',ipDSVar,nDens)
      Call GetMem('DS   ','Free','Real',ipDS   ,nDens)
      Call GetMem('DVar ','Free','Real',ipDVar ,nDens)
      Call GetMem('D0   ','Free','Real',ipD0   ,nDens)
*
*      If (ip_Z_p_k.ne.ip_Dummy) Call Free_Work(ip_Z_p_k)
*
*     Call qExit('CloseP')
      Return
      End
