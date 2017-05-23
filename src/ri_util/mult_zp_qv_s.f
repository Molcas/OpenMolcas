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
      Subroutine Mult_Zp_Qv_s(Zp,nZp,Qv,nQv,Rv,n_Rv,nVec,nMuNu,nI,
     &                        nIrrep,QMode)
************************************************************************
*     Author:   F. Aquilante                                           *
*                                                                      *
*     Qv: is a symmetry blocked square matrix                          *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Real*8 Zp(nZp), Qv(nQv), Rv(n_Rv)
      Integer nVec(0:nIrrep-1), nMuNu(0:nIrrep-1), nI(0:nIrrep-1)
      Character QMode*1, Name_Q*6
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

      Do iIrrep = 0, nIrrep-1
         nI_     = nI(iIrrep)
         If (iIrrep.eq.0) nI_=nI_-1
         nMuNu_  = nMuNu(iIrrep)
         If (nMuNu_.le.0 .or. nI_.le.0) Go To 999
*
         iSeed=55+iIrrep
         Lu_Q=IsFreeUnit(iSeed)
         Write(Name_Q,'(A4,I2.2)') 'QVEC',iIrrep
         Call DaName_MF_WA(Lu_Q,Name_Q)
         iAddr=0
*
         iOffR2=iOffR
         iOffA2=iOffA
         mQv = nI_*nVec(iIrrep)
         Do While (mQv.ge.nI_)
*
            nK=Min(mQv,nQv)/nI_
            lQv=nI_*nK
            Call dDaFile(Lu_Q,2,Qv,lQv,iAddr)
*
            Call A_3C_Qv_s(Zp(iOffA2),
     &                     Qv,
     &                     Rv(iOffR2),nMuNu_,nI_,nK,QMode)
            mQv = mQv - lQv
            iOffR2 = iOffR2 + lstepR*nMuNu_*nK
            iOffA2 = iOffA2 + lstepA*nMuNu_*nK
         End Do
*
         iOffA = iOffA + nMuNu_*nI_
         iOffR = iOffR + nMuNu_*nVec(iIrrep)
         Call DaClos(Lu_Q)
 999     Continue
      End Do
      Return
      End
