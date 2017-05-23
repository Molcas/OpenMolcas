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
      Subroutine PLF_LDF_JK_2P_1(TInt,nTInt,Map,
     &                           AOint,ijkl,iCmp,jCmp,kCmp,lCmp,
     &                           iAO,iAOst,iBas,jBas,kBas,lBas,kOp,
     &                           iAOtSO,nAOtSO)
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
*          Modified for Local DF, Thomas Bondo Pedersen, September 2010*
*                                                                      *
************************************************************************
      Implicit None
      Integer nTInt
      Real*8  TInt(nTInt)
      Integer Map(4)
      Integer ijkl, iCmp, jCmp, kCmp, lCmp
      Real*8  AOInt(ijkl,iCmp,jCmp,kCmp,lCmp)
      Integer iAO(4), iAOst(4)
      Integer iBas, jBas, kBas, lBas
      Integer kOp(4)
      Integer nAOtSO
      Integer iAOtSO(nAOtSO,0:7)
#include "localdf_bas.fh"
#include "localdf_int2.fh"
#include "WrkSpc.fh"

      Character*15 SecNam
      Parameter (SecNam='PLF_LDF_JK_2P_1')

      Integer i1, i2, i3, i4
      Integer iShlJ, iShlL
      Integer jSO, lSO
      Integer jSOj, lSOl
      Integer ll, jj
      Integer nijkl
      Integer JL

      Integer i, j
      Integer iShlSO
      Integer iRow
      Integer iCol
#if defined (_DEBUG_)
      Integer iSOShl
      iSOShl(i)=iWork(ip_iSOShl-1+i)
#endif
      iShlSO(i)=iWork(ip_iShlSO-1+i)
      iRow(i,j)=iWork(ip_AB_IndxG-1+l_AB_IndxG_1*(j-1)+i)
      iCol(i,j)=iWork(ip_CD_IndxG-1+l_CD_IndxG_1*(j-1)+i)

      If (Map(1).eq.1 .and. Map(2).eq.2 .and.
     &    Map(3).eq.3 .and. Map(4).eq.4) Then
#if defined (_DEBUG_)
         If (iSOShl(iAOtSO(iAO(1)+1,kOp(1))+iAOst(1)).ne.SHA) Then
            Call WarningMessage(2,SecNam//': Shell problem [1.1]')
            Call LDF_Quit(1)
         End If
         If (iSOShl(iAOtSO(iAO(2)+1,kOp(2))+iAOst(2)).ne.SHB) Then
            Call WarningMessage(2,SecNam//': Shell problem [2.1]')
            Call LDF_Quit(1)
         End If
         If (iSOShl(iAOtSO(iAO(3)+1,kOp(3))+iAOst(3)).ne.SHC) Then
            Call WarningMessage(2,SecNam//': Shell problem [3.1]')
            Call LDF_Quit(1)
         End If
         If (iSOShl(iAOtSO(iAO(4)+1,kOp(4))+iAOst(4)).ne.SHD) Then
            Call WarningMessage(2,SecNam//': Shell problem [4.1]')
            Call LDF_Quit(1)
         End If
         If (iCmp*iBas.ne.1) Then
            Call WarningMessage(2,
     &                        SecNam//': Shell dimension problem [1.1]')
            Call LDF_Quit(1)
         End If
         If (kCmp*kBas.ne.1) Then
            Call WarningMessage(2,
     &                        SecNam//': Shell dimension problem [3.1]')
            Call LDF_Quit(1)
         End If
#endif
         iShlJ=SHB
         iShlL=SHD
         i1=1
         i3=1
         Do i4=1,lCmp
            lSO=iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
            Do i2=1,jCmp
               jSO=iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
               nijkl=0
               Do lSOl=lSO,lSO+lBas-1
                  ll=iShlSO(lSOl)
                  If (iCol(ll,iShlL).gt.0) Then
                     Do jSOj=jSO,jSO+jBas-1
                        jj=iShlSO(jSOj)
                        nijkl=nijkl+1
                        If (iRow(jj,iShlJ).gt.0) Then
                           JL=nAB*(iCol(ll,iShlL)-1)+iRow(jj,iShlJ)
                           TInt(JL)=AOInt(nijkl,i1,i2,i3,i4)
                        End If
                     End Do
                  Else
                     nijkl=nijkl+jBas
                  End If
               End Do
            End Do
         End Do
      Else If (Map(1).eq.3 .and. Map(2).eq.4 .and.
     &         Map(3).eq.1 .and. Map(4).eq.2) Then
#if defined (_DEBUG_)
         If (iSOShl(iAOtSO(iAO(1)+1,kOp(1))+iAOst(1)).ne.SHC) Then
            Call WarningMessage(2,SecNam//': Shell problem [1.2]')
            Call LDF_Quit(1)
         End If
         If (iSOShl(iAOtSO(iAO(2)+1,kOp(2))+iAOst(2)).ne.SHD) Then
            Call WarningMessage(2,SecNam//': Shell problem [2.2]')
            Call LDF_Quit(1)
         End If
         If (iSOShl(iAOtSO(iAO(3)+1,kOp(3))+iAOst(3)).ne.SHA) Then
            Call WarningMessage(2,SecNam//': Shell problem [3.2]')
            Call LDF_Quit(1)
         End If
         If (iSOShl(iAOtSO(iAO(4)+1,kOp(4))+iAOst(4)).ne.SHB) Then
            Call WarningMessage(2,SecNam//': Shell problem [4.2]')
            Call LDF_Quit(1)
         End If
         If (iCmp*iBas.ne.1) Then
            Call WarningMessage(2,
     &                        SecNam//': Shell dimension problem [1.2]')
            Call LDF_Quit(1)
         End If
         If (kCmp*kBas.ne.1) Then
            Call WarningMessage(2,
     &                        SecNam//': Shell dimension problem [3.2]')
            Call LDF_Quit(1)
         End If
#endif
         iShlJ=SHD
         iShlL=SHB
         i1=1
         i3=1
         Do i4=1,lCmp
            lSO=iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
            Do i2=1,jCmp
               jSO=iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
               nijkl=0
               Do lSOl=lSO,lSO+lBas-1
                  ll=iShlSO(lSOl)
                  If (iRow(ll,iShlL).gt.0) Then
                     Do jSOj=jSO,jSO+jBas-1
                        jj=iShlSO(jSOj)
                        nijkl=nijkl+1
                        If (iCol(jj,iShlJ).gt.0) Then
                           JL=nAB*(iCol(jj,iShlJ)-1)+iRow(ll,iShlL)
                           TInt(JL)=AOInt(nijkl,i1,i2,i3,i4)
                        End If
                     End Do
                  Else
                     nijkl=nijkl+jBas
                  End If
               End Do
            End Do
         End Do
      Else
         Call WarningMessage(2,
     &                    SecNam//': Shell combination not implemented')
         Call LDF_Quit(1)
      End If

#ifndef _DEBUG_
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iBas)
         Call Unused_integer(kBas)
      End If
#endif
      End
