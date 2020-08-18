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
      SubRoutine ClsSew
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
      use External_Centers
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
      If (Seward_Status.eq.InActive) Return
*
      Call Term_Ints(.False.,.True.)
      Call Free_RctFld(iXPolType)
      Call Free_HerRW()
*
      If (Allocated(RSph)) Call mma_deallocate(RSph)
      If (Allocated(ipSph)) Call mma_deallocate(ipSph)
*                                                                      *
************************************************************************
*                                                                      *
      If (Info_Status.eq.Active) Then
         Call GetMem('Info','Free','Real',LctInf,nInfo)
         Info_Status=InActive
*                                                                      *
*        Process "external" centers
*
         If (Allocated(EF_Centers)) Then
            Call mma_deallocate(EF_Centers)
            nEF=0
         End If
         If (Allocated(OAM_Center)) Call mma_deallocate(OAM_Center)
         If (Allocated(OMQ_Center)) Call mma_deallocate(OMQ_Center)
         If (Allocated(DMS_Centers)) Then
            Call mma_deallocate(DMS_Centers)
            nDMS=0
         End If
         If (Allocated(Wel_Info)) Then
            Call mma_deallocate(Wel_Info)
            nWel=0
         End If
         If (Allocated(AMP_Center)) Call mma_deallocate(AMP_Center)
         If (Allocated(RP_Centers)) Then
            Call mma_deallocate(RP_Centers)
            nRP=0
         End If
         If (Allocated(XF)) Then
            Call mma_deallocate(XF)
            Call mma_deallocate(XMolnr)
            Call mma_deallocate(XEle)
            nData_XF=0
            nXF=0
            nXMolnr=0
         End If
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call Free_iSD()
      Call Freek2()
      Call CloseR()
*
      If (lEFP) Then
         Deallocate(FRAG_TYPE)
         Deallocate(ABC)
         Deallocate(EFP_COORS)
         lEFP=.FALSE.
      End If
*
      Seward_Status=InActive
      Return
      End
c
c This code originally was included into ClsSew
c occasionally it should be separated
      Subroutine DumpSagit()
      Implicit Real*8 (A-H,O-Z)
c      Character*8 sagit
c      Call getenvf('MOLCAS_SAGIT',sagit)
c      If (sagit(1:1).eq.'y'.or.sagit(1:1).eq.'Y') Then
CVV: dump info from runfile into ORB.std
C    note that changes in info.fh
C    should be reflected in sagit
        iutemp=16
        iutemp=isfreeunit(iutemp)
        Call molcas_open(iutemp,'ORB.std')
        Call Koor2file(iutemp)
        Call Basi2file(iutemp)
        close(iutemp)
c      End If
      End
