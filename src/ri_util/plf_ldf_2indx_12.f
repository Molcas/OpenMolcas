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
      Subroutine PLF_LDF_2Indx_12(TInt,nTInt,
     &                            AOint,ijkl,
     &                            iCmp,jCmp,kCmp,lCmp,
     &                            iAO,iAOst,iBas,jBas,kBas,lBas,kOp)
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
#include "localdf_int.fh"
#include "WrkSpc.fh"

#if defined (_DEBUG_)
      Character*16 SecNam
      Parameter (SecNam='PLF_LDF_2Indx_12')
#endif

      Integer iShlK, iShlKL
      Integer i1, i2, i3, i4
      Integer jSO, kSO, lSO
      Integer jSOj, kSOk, lSOl
      Integer kk, ll, kl0, kl
      Integer nijkl, ip0

      Integer i, j
      Integer iShlSO, nBasSh, IndxG2
#if defined (_DEBUG_)
      Integer iSOShl
      iSOShl(i)=iWork(ip_iSOShl-1+i)
#endif
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      iShlSO(i)=iWork(ip_iShlSO-1+i)
      IndxG2(i,j)=iWork(ip_IndxG2-1+l_IndxG2_1*(j-1)+i)

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
      If (iCmp*iBas.ne.1) Then
         Call WarningMessage(2,SecNam//': Shell dimension problem [1]')
         Call LDF_Quit(1)
      End If
      If (jCmp*jBas.ne.nBasSh(SHB)) Then
         Call WarningMessage(2,SecNam//': Shell dimension problem [2]')
         Call LDF_Quit(1)
      End If
      If (kCmp*kBas.ne.nBasSh(SHC)) Then
         Call WarningMessage(2,SecNam//': Shell dimension problem [3]')
         Call LDF_Quit(1)
      End If
      If (lCmp*lBas.ne.nBasSh(SHD)) Then
         Call WarningMessage(2,SecNam//': Shell dimension problem [4]')
         Call LDF_Quit(1)
      End If
#endif

      iShlK=SHC
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
                     If (IndxG2(kl,iShlKL).gt.0) Then
                        ip0=nRow_G*(IndxG2(kl,iShlKL)-1)+iOffuv
                        Do jSOj=jSO,jSO+jBas-1
                           nijkl=nijkl+1
                           TInt(ip0+iShlSO(jSOj))=AOInt(nijkl,
     &                                                  i1,i2,i3,i4)
                        End Do
                     Else
                        nijkl=nijkl+jBas
                     End If
                  End Do
               End Do
            End Do
         End Do
      End Do

#ifndef _DEBUG_
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iBas)
#endif
      End
