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
*               2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine PLF_LDF_Gmax_S(TInt,nTInt,
     &                          AOint,ijkl,
     &                          iCmp,jCmp,kCmp,lCmp,
     &                          iBas,jBas,kBas,lBas)
************************************************************************
*                                                                      *
*  object: to sift and index the petite list format integrals.         *
*                                                                      *
*          the indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*                                                                      *
*  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
*          May '90                                                     *
*          Modified for Local DF, Thomas Bondo Pedersen, Oct. 2010     *
*                                                                      *
************************************************************************
      Implicit None
      Integer nTInt
      Real*8  TInt(nTInt)
      Integer ijkl, iCmp, jCmp, kCmp, lCmp
      Real*8  AOInt(ijkl,iCmp,jCmp,kCmp,lCmp)
      Integer iBas, jBas, kBas, lBas
#if defined (_DEBUG_)
#include "localdf_bas.fh"
#include "localdf_int3.fh"
#include "WrkSpc.fh"

      Character*14 SecNam
      Parameter (SecNam='PLF_LDF_Gmax_S')
#endif

      Integer i4
      Integer lSOl
      Integer nijkl

#if defined (_DEBUG_)
      Integer i
      Integer nBasSh
      nBasSh(i)=iWork(ip_nBasSh-1+i)

      If (iCmp*iBas.ne.1) Then
         Call WarningMessage(2,SecNam//': Shell dimension problem [1]')
         Call LDF_Quit(1)
      End If
      If (jCmp*jBas.ne.nBasSh(SHB)) Then
         Call WarningMessage(2,SecNam//': Shell dimension problem [2]')
         Call LDF_Quit(1)
      End If
      If (kCmp*kBas.ne.1) Then
         Call WarningMessage(2,SecNam//': Shell dimension problem [3]')
         Call LDF_Quit(1)
      End If
      If (lCmp*lBas.ne.nBasSh(SHD)) Then
         Call WarningMessage(2,SecNam//': Shell dimension problem [4]')
         Call LDF_Quit(1)
      End If
      If (SHB.ne.SHD) Then
         Call WarningMessage(2,SecNam//': SHB != SHD')
         Call LDF_Quit(1)
      End If
      If (nTInt.lt.2) Then
         Call WarningMessage(2,SecNam//' nTInt < 2')
         Call LDF_Quit(1)
      End If
#endif

      TInt(1)=0.0d0
      TInt(2)=0.0d0
      Do i4=1,lCmp
         Do lSOl=1,lBas
            nijkl=lBas*(lSOl-1)+lSOl
            TInt(1)=max(TInt(1),AOInt(nijkl,1,i4,1,i4))
            TInt(2)=TInt(2)+AOInt(nijkl,1,i4,1,i4)
         End Do
      End Do

#ifndef _DEBUG_
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iBas)
         Call Unused_integer(jBas)
         Call Unused_integer(kBas)
      End If
#endif
      End
