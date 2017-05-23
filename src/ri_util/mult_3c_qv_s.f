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
* Copyright (C) Francesco Aquilante                                    *
************************************************************************
      Subroutine Mult_3C_Qv_s(A_3C,nA_3C,Qv,nQv,Rv,n_Rv,nVec,iOff_3C,
     &                        nIrrep,Out_of_Core,Lu_Q,QMode)
************************************************************************
*     Author:   F. Aquilante                                           *
*                                                                      *
*     Qv: is a symmetry blocked square matrix                          *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Real*8 A_3C(nA_3C), Qv(nQv), Rv(n_Rv)
      Integer iOff_3C(3,0:nIrrep-1), nVec(0:7), Lu_Q(0:nIrrep-1)
      Logical Out_of_Core
      Character*1 QMode
*                                                                      *
************************************************************************
*                                                                      *
      lstepA=0
      lstepR=1
      If (Qmode.eq.'T') Then
         Call FZero(Rv,n_Rv)
         lstepA=1
         lstepR=0
      End If
*
      iOffA=1
      iOffR=1

      If (Out_of_Core) Then
*                                                                      *
************************************************************************
*                                                                      *
         Do iIrrep = 0, nIrrep-1
            nI     = iOff_3C(2,iIrrep)
            nMuNu  = iOff_3C(1,iIrrep)
            If (nMuNu.le.0 .or. nI.le.0) Go To 999
*
            iOffR2=iOffR
            iOffA2=iOffA
            mQv = nI*nVec(iIrrep)
            iAddr=0
            Do While (mQv.ge.nI)
*
               nK=Min(mQv,nQv)/nI
               lQv=nI*nK
               Call dDaFile(Lu_Q(iIrrep),2,Qv,lQv,iAddr)
*
               Call A_3C_Qv_s(A_3C(iOffA2),
     &                        Qv,
     &                        Rv(iOffR2),nMuNu,nI,nK,QMode)
               mQv = mQv - lQv
               iOffR2 = iOffR2 + lstepR*nMuNu*nK
               iOffA2 = iOffA2 + lstepA*nMuNu*nK
            End Do
*
            iOffA = iOffA + nMuNu*nI
            iOffR = iOffR + nMuNu*nVec(iIrrep)
 999        Continue
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      Else  ! In-Core
*                                                                      *
************************************************************************
*                                                                      *
         iOffQ=1
         Do iIrrep = 0, nIrrep-1
            nI     = iOff_3C(2,iIrrep)
            nMuNu  = iOff_3C(1,iIrrep)
            If (nMuNu.le.0 .or. nI.le.0) Go To 998
*
            Call A_3C_Qv_s(A_3C(iOffA),
     &                     Qv(iOffQ),
     &                     Rv(iOffR),nMuNu,nI,nVec(iIrrep),QMode)
 998        Continue
            iOffA = iOffA + nMuNu*nI
            iOffR = iOffR + nMuNu*nVec(iIrrep)
            iOffQ = iOffQ + nI*nVec(iIrrep)
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
