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
      Subroutine PLF_LDF_3Indx_1(TInt,nTInt,
     &                           AOint,ijkl,
     &                           iCmp,jCmp,kCmp,lCmp,
     &                           iAO,iAOst,iBas,jBas,kBas,lBas,kOp)
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
      use SOAO_Info, only: iAOtSO
      Implicit None
      Integer nTInt
      Real*8  TInt(nTInt)
      Integer ijkl, iCmp, jCmp, kCmp, lCmp
      Real*8  AOInt(ijkl,iCmp,jCmp,kCmp,lCmp)
      Integer iAO(4), iAOst(4)
      Integer iBas, jBas, kBas, lBas
      Integer kOp(4)
#include "localdf_bas.fh"
#include "localdf_int3.fh"
#include "WrkSpc.fh"

#if defined (_DEBUG_)
      Character*15 SecNam
      Parameter (SecNam='PLF_LDF_3Indx_1')
#endif

      Integer i1, i2, i3, i4
      Integer iShlI, iShlJ, iShlL
      Integer iSO, jSO, lSO
      Integer iSOi, jSOj, lSOl
      Integer ii, jj, ll
      Integer iRow, iRow1, iCol
      Integer nijkl, ip, ip0

      Integer i
      Integer iShlSO, nBasSh
#if defined (_DEBUG_)
      Integer iSOShl
      iSOShl(i)=iWork(ip_iSOShl-1+i)
#endif
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      iShlSO(i)=iWork(ip_iShlSO-1+i)

#if defined (_DEBUG_)
      If (iSOShl(iAOtSO(iAO(1)+1,kOp(1))+iAOst(1)).ne.SHA) Then
         Call WarningMessage(2,SecNam//': Shell problem [1]')
         Call LDF_Quit(1)
      End If
      If (iSOShl(iAOtSO(iAO(2)+1,kOp(2))+iAOst(2)).ne.SHB) Then
         Call WarningMessage(2,SecNam//': Shell problem [2]')
         Call LDF_Quit(1)
      End If
      If (iSOShl(iAOtSO(iAO(3)+1,kOp(3))+iAOst(3)).ne.SHC) Then
         Call WarningMessage(2,SecNam//': Shell problem [3]')
         Call LDF_Quit(1)
      End If
      If (iSOShl(iAOtSO(iAO(4)+1,kOp(4))+iAOst(4)).ne.SHD) Then
         Call WarningMessage(2,SecNam//': Shell problem [4]')
         Call LDF_Quit(1)
      End If
      If (iCmp*iBas.ne.nBasSh(SHA)) Then
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
#endif

      iShlI=SHA
      iShlJ=SHB
      iShlL=SHC
      i3=1
      Do i4=1,lCmp
         lSO=iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
         Do i2=1,jCmp
            jSO=iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
            Do i1=1,iCmp
               iSO=iAOtSO(iAO(1)+i1,kOp(1))+iAOst(1)
               nijkl=0
               Do lSOl=lSO,lSO+lBas-1
                  ll=iShlSO(lSOl)
                  iCol=iCol0+ll
                  ip0=nRow*(iCol-1)
                  Do jSOj=jSO,jSO+jBas-1
                     jj=iShlSO(jSOj)
                     iRow1=iRow0+nBasSh(iShlI)*(jj-1)
                     Do iSOi=iSO,iSO+iBas-1
                        ii=iShlSO(iSOi)
                        iRow=iRow1+ii
                        ip=ip0+iRow
                        nijkl=nijkl+1
                        TInt(ip)=AOInt(nijkl,i1,i2,i3,i4)
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do

c Avoid unused argument warnings
      If (.False.) Call Unused_integer(kBas)
      End
