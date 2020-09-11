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
      use Real_Spherical, only: Sphere_Free
      use EFP_module
      use External_Centers
      use Basis_Info
      use Center_Info
      Use SOAO_Info
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
      Call Sphere_Free()
      Call SOAO_Info_Free()
      Call Basis_Info_Free()
      Call Center_Info_Free()
      Call External_Centers_Free()
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
