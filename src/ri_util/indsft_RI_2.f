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
************************************************************************
      SubRoutine IndSft_RI_2(iCmp,iShell,iBas,jBas,kBas,lBas,
     &                       Shijij, iAO, iAOst, ijkl,SOint,nSOint,
     &                       iSOSym,nSOs,TInt,nTInt,iOff,
     &                       iSO2Ind,iOffA)
************************************************************************
*  object: to sift and index the SO integrals.                         *
*                                                                      *
*          the indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
*          april '90                                                   *
*                                                                      *
************************************************************************
      use Basis_Info, only: nBas
      use SOAO_Info, only: iAOtSO
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
#include "srt0.fh"
#include "srt1.fh"
#include "WrkSpc.fh"
*
      Real*8 SOint(ijkl,nSOint), TInt(nTInt)
      Integer iCmp(4), iShell(4), iAO(4), iOffA(4,0:7),
     &        iAOst(4), iSOSym(2,nSOs), iSO2Ind(nSOs)
      Integer iOff(0:7)
      Logical Shijij, qijij
*     local array
      Integer jSym(0:7), lSym(0:7)
#ifdef _DEBUGPRINT_
      Data tr1,tr2/0.0d0,0.0d0/
      Save tr1,tr2
#endif
*                                                                      *
************************************************************************
*                                                                      *
      iTri(i,j)=Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      irout = 39
      iprint = nprint(irout)
*                                                                      *
************************************************************************
*                                                                      *
      k12=0
      k34=0
#ifdef _DEBUGPRINT_
      If (iPrint.ge.49) Then
         r1=DDot_(ijkl*nSOInt,SOInt,1,[One],0)
         r2=DDot_(ijkl*nSOInt,SOInt,1,SOInt,1)
         tr1=tr1+r1
         tr2=tr2+r2
         Write (6,*) ' Sum=',r1,tr1
         Write (6,*) ' Dot=',r2,tr2
         Call RecPrt(' in indsft:SOint ',' ',SOint,ijkl,nSOint)
      End If
#endif
      memSO2 = 0
*
*
*     quadruple loop over elements of the basis functions angular
*     description. loops are reduced to just produce unique SO integrals
*     observe that we will walk through the memory in AOint in a
*     sequential way.
*
      i1 = 1
      i3 = 1
      j1 = 0
      j3 = 0
      Do i2 = 1, iCmp(2)
         Do 201 j = 0, nIrrep-1
            ix = 0
            If (iAOtSO(iAO(2)+i2,j)>0) ix = 2**j
            jSym(j) = ix
201      Continue
         If (iShell(2).gt.iShell(1)) then
            i12 = iCmp(2)*(i1-1) + i2
         else
            i12 = iCmp(1)*(i2-1) + i1
         End If
         Do 400 i4 = 1, iCmp(4)
            Do 401 j = 0, nIrrep-1
               ix = 0
               If (iAOtSO(iAO(4)+i4,j)>0) ix = 2**j
               lSym(j) = ix
401         Continue
            If (iShell(4).gt.iShell(3)) then
               i34 = iCmp(4)*(i3-1) + i4
            else
               i34 = iCmp(3)*(i4-1) + i3
            End If
            If (Shijij .and. i34.gt.i12) go to 400
            qijij = Shijij .and. i12.eq.i34
C           Write (6,*) 'i1,i2,i3,i4=',i1,i2,i3,i4
*
*      loop over Irreps which are spanned by the basis function.
*      again, the loop structure is restricted to ensure unique
*      integrals.
*
             Do 210 j2 = 0, nIrrep-1
                If (jSym(j2).eq.0) go to 210
                j12 = ieor(j1,j2)
                If (qijij) then
                   If (iShell(1).gt.iShell(2)) then
                       k12 = nIrrep*j1 + j2+1
                   else
                       k12 = nIrrep*j2 + j1+1
                   End If
                End If
*
                iOffA_= iOffA(1,j2)
                iOffB_= iOffA(3,j2)
                If (j2.ne.0) iOffB_=iOffB_+1
                mm_ = iOffA(4,j2)
                nn  = mm_ - iOffA(2,j2)
                mx  = nn*(nn+1)/2
*
                j4 = ieor(j12,j3)
                If (lSym(j4).eq.0) go to 210
                If (qijij) then
                   If (iShell(3).gt.iShell(4)) then
                      k34 = nIrrep*j3 + j4+1
                   else
                      k34 = nIrrep*j4 + j3+1
                   End If
                   If (k34.gt.k12) go to 210
                End If
*
                memSO2 = memSO2 + 1
                If ( (nSkip(j2+1)+nSkip(j4+1) ).ne.0 ) GoTo 210
*
*               Compute absolute starting SO index
                jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
                lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)
*
                nijkl = 0
                Do lSOl = lSO, lSO+lBas-1
                   Do jSOj = jSO, jSO+jBas-1
                      nijkl = nijkl + 1
                      AInt=SOint(nijkl,memSO2)
                      iSO = jSOj-nBas(j2)
                      kSO = lSOl-nBas(j4)
*
                      iSO = iSO2Ind(iSO+iOffB_) + nn
                      ij= iTri(iSO,kSO) - mx + iOffA_
                      TInt(ij)=AInt
*
                   End Do
                End Do
*
210       Continue
*
400      Continue
      End Do
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iBas)
         Call Unused_integer(kBas)
         Call Unused_integer_array(iSOSym)
         Call Unused_integer_array(iOff)
      End If
      End
