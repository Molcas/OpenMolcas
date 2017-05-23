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
      Subroutine PLF_LDF_JK_2P_3(TInt,nTInt,Map,
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
*          Modified for Local DF, Thomas Bondo Pedersen, June 2010     *
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
      Parameter (SecNam='PLF_LDF_JK_2P_3')

      Integer i1, i2, i3, i4
      Integer iShlI, iShlJ, iShlK, iShlL
      Integer iSO, jSO, kSO, lSO
      Integer iSOi, jSOj, lSOl, kSOk
      Integer ii, jj, kk, ll
      Integer iShlKL, kl, kl0
      Integer iShlIJ, ij, ij0
      Integer nijkl
      Integer jkL, ijL

      Integer i, j
      Integer nBasSh, iShlSO, iRow, iCol
#if defined (_DEBUG_)
      Integer iSOShl
      iSOShl(i)=iWork(ip_iSOShl-1+i)
#endif
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      iShlSO(i)=iWork(ip_iShlSO-1+i)
      iRow(i,j)=iWork(ip_AB_IndxG2-1+l_AB_IndxG2_1*(j-1)+i)
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
         If (kCmp*kBas.ne.1) Then
            Call WarningMessage(2,
     &                        SecNam//': Shell dimension problem [1.1]')
            Call LDF_Quit(1)
         End If
#endif
         iShlI=SHA
         iShlJ=SHB
         iShlIJ=SPAB
         iShlL=SHD
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
                     If (iCol(ll,iShlL).gt.0) Then
                        Do jSOj=jSO,jSO+jBas-1
                           jj=iShlSO(jSOj)
                           ij0=nBasSh(iShlI)*(jj-1)
                           Do iSOi=iSO,iSO+iBas-1
                              ii=iShlSO(iSOi)
                              ij=ij0+ii
                              nijkl=nijkl+1
                              If (iRow(ij,iShlIJ).gt.0) Then
                                 ijL=nAB*(iCol(ll,iShlL)-1)
     &                              +iRow(ij,iShlIJ)
                                 TInt(ijL)=AOInt(nijkl,i1,i2,i3,i4)
                              End If
                           End Do
                        End Do
                     Else
                        nijkl=nijkl+iBas*jBas
                     End If
                  End Do
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
#endif
         iShlJ=SHD
         iShlK=SHA
         iShlL=SHB
         iShlKL=SPAB
         i1=1
         Do i4=1,lCmp
            lSO=iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
            Do i3=1,kCmp
               kSO=iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
               Do i2=1,jCmp
                  jSO=iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
                  nijkl=0
                  Do lSOl=lSO,lSO+lBas-1
                     ll=iShlSO(lSOl)
                     kl0=nBasSh(iShlK)*(ll-1)
                     Do kSOk=kSO,kSO+kBas-1
                        kk=iShlSO(kSOk)
                        kl=kl0+kk
                        If (iRow(kl,iShlKL).gt.0) Then
                           Do jSOj=jSO,jSO+jBas-1
                              jj=iShlSO(jSOj)
                              nijkl=nijkl+1
                              If (iCol(jj,iShlJ).gt.0) Then
                                 Jkl=nAB*(iCol(jj,iShlJ)-1)
     &                              +iRow(kl,iShlKL)
                                 TInt(Jkl)=AOInt(nijkl,i1,i2,i3,i4)
                              End If
                           End Do
                        Else
                           nijkl=nijkl+jBas
                        End If
                     End Do
                  End Do
               End Do
            End Do
         End Do
      Else
         Call WarningMessage(2,
     &                   SecNam//': shell combination not implemented!')
         Call LDF_Quit(1)
      End If

      End
