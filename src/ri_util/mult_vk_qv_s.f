************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Mult_Vk_Qv_s(V_k,nV_k,Qv,nQv,V_kQ,nV_kQ,nBas_Aux,
     &                        nVec,nIrrep,QMode)
      Implicit Real*8 (a-h,o-z)
      Real*8 Qv(nQv), V_k(nV_k), V_kQ(nV_kQ)
      Integer nBas_Aux(0:nIrrep-1), nVec(0:nIrrep-1)
      Logical Out_of_Core
      Character QMode*1, Name_Q*6
*                                                                      *
************************************************************************
*                                                                      *
      nMuNu=1
      kp_V_k=1
      lstepA=0
      lstepB=1
      If (Qmode.eq.'T') Then
         Call FZero(V_kQ,nV_kQ)
         lstepA=1
         lstepB=0
      End If

      Do iIrrep = 0, 0  ! loop is wisely restricted to tot. symm. irrep
         iSeed=55+iIrrep
         Lu_Q=IsFreeUnit(iSeed)
         Write(Name_Q,'(A4,I2.2)') 'QVEC',iIrrep
         Call DaName_MF_WA(Lu_Q,Name_Q)
         iAddr=0
*
         nI=nBas_Aux(iIrrep)
         If (iIrrep.eq.0) nI=nI-1
         nJ=nVec(iIrrep)
         mQv=nI*nJ
         Out_of_Core = mQv .gt. nQv
*                                                                      *
************************************************************************
*                                                                      *
         If (.Not.Out_of_Core) Then
*
*           in-core case
*
            Call dDaFile(Lu_Q,2,Qv,mQv,iAddr)
*
            Call A_3C_Qv_s(V_k(kp_V_k),
     &                     Qv,
     &                     V_kQ,nMuNu,nI,nJ,QMode)
*
         Else
*
*           Out-of-core case
*
            iOffA=kp_V_k
            iOffB=iOffA
            Do While (mQv.ge.nI)
*
               nK=Min(mQv,nQv)/nI
               lQv=nI*nK
               Call dDaFile(Lu_Q,2,Qv,lQv,iAddr)
*
               Call A_3C_Qv_s(V_k(iOffA),
     &                        Qv,
     &                        V_kQ(iOffB),nMuNu,nI,nK,QMode)
               mQv = mQv - lQv
               iOffA = iOffA + lstepA*nK
               iOffB = iOffB + lstepB*nK
            End Do
*
         End If
*                                                                      *
************************************************************************
*                                                                      *
         Call DaClos(Lu_Q)
*
         kp_V_k = kp_V_k + nBas_Aux(iIrrep)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
