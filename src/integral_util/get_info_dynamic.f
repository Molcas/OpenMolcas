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
      SubRoutine Get_Info_Dynamic(Info,nInfo,ioffr,icase)
************************************************************************
*                                                                      *
* Object: to read all input information on the file INFO.              *
*                                                                      *
* Called from: Alaska                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              OpnCom                                                  *
*              ClsCom                                                  *
*              RdCom                                                   *
*              SetUp_RW                                                *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January 1992                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "status.fh"
*
      Logical     Found
*
*     Prologue
*
*     Call qEnter('GetInf')
*     Load the dynamic input area.
*
      Call qpg_dArray('SewXInfo',Found,Len2)
      If (.not.Found .or. Len2.eq.0) Then
         Call WarningMessage(1,
     &                        'Get_info_dynamic: Did not find SewXInfo')
      End If
      Len=Len2
      nInfo=Len
*     mod by M.Schuetz: LenInf used to broadcast dynamic input
*     area to servers (parallel distributed SCF)
*     LenInf is a member of IInfo common block
      LenInf=nInfo
      If (Info_Status.eq.Active) Then
         Call WarningMessage(2,'Info_Status already active!')
         Call Abend()
      End If
      If (Info_Status.ne.InActive) Then
         Call WarningMessage(2,'Info_Status not properly set!')
         Call Abend()
      End If
      Info_Status=Active
      Call GetMem(' SewXInfo ','Allo','Real',Info,nInfo)
      Call FZero(Work(Info),nInfo)
      Call Get_dArray('SewXInfo',Work(Info),Len)
*
*     Update all pointers to accomodate the new position in memory
*
      iOld = LctInf
*     mod by M.Schuetz: LctInf used to update pointers on servers
*     (parallel distributed SCF) after broadcast
*     LctInf is a member of IInfo common block
      LctInf=Info
      Call Gen_RelPointers(LctInf-1)
*
*     Epilogue, end
*
*     Call qExit('GetInf')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(ioffr)
         Call Unused_integer(icase)
      End If
      End
