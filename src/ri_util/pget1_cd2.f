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
* Copyright (C) 1992,2007, Roland Lindh                                *
*               2009, Francesco Aquilante                              *
************************************************************************
      SubRoutine PGet1_CD2(PAO,ijkl,nPAO,iCmp,
     &                 iAO,iAOst,Shijij,iBas,jBas,kBas,lBas,kOp,
     &                 ExFac,CoulFac,PMax,V_k,U_k,mV_k,Z_p_K,nnP1)
************************************************************************
*  Object: to assemble the 2nd order density matrix of a SCF wave      *
*          function from the 1st order density.                        *
*                                                                      *
*          The indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
* Called from: PGet0                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             January '92.                                             *
*                                                                      *
*             Modified for Cholesky 1-center gradients May 2007 by     *
*             R. Lindh                                                 *
*                                                                      *
*             Modified for RI-HF/CAS, Dec 2009 (F. Aquilante)          *
*                                                                      *
************************************************************************
*     use pso_stuff
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "print.fh"
#include "WrkSpc.fh"
#include "chomp2g_alaska.fh"
#include "exterm.fh"
      Real*8 PAO(ijkl,nPAO), V_k(mV_k), U_K(mV_K), Z_p_K(nnP1,mV_K),
     &       Fac_ij,Fac_kl
      Integer iAO(4), kOp(4), iAOst(4), iCmp(4)
      Logical Shijij
      External mn2K
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 39
      iPrint = nPrint(iRout)
*define _DEBUG_
#ifdef _DEBUG_
      Call qEnter('PGet1_CD2')
      Call RecPrt('PGet1_CD2: V_k',' ',V_k,1,mV_k)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Quadruple loop over elements of the basis functions angular
*     description.
*     Observe that we will walk through the memory in PAO in a
*     sequential way.
*
C     Fac = One / Four

      Call CWTime(Cpu1,Wall1)

      Fac = One
      PMax=Zero
      iPAO=0

      iSym = 1
      jSym = 1
      kSym = 1
      lSym = 1
      iSO = 1

      If(ExFac .eq. Zero) Then

         Do i1 = 1, iCmp(1)
            Do i2 = 1, iCmp(2)
               Do i3 = 1, iCmp(3)
                  Do i4 = 1, iCmp(4)
*
*               Unfold the way the eight indicies have been reordered.
                     iSO = iAOtSO(iAO(1)+i1,kOp(1))+iAOst(1)
                     jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
                     kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
                     lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
*
                     iPAO = iPAO + 1
                     nijkl = 0
*
                     Do lAOl = 0, lBas-1
                        lSOl = lSO + lAOl
                        Do kAOk = 0, kBas-1
                           kSOk = kSO + kAOk
                           Do jAOj = 0, jBas-1
                              jSOj = jSO + jAOj
                              Do iAOi = 0, iBas-1
                                 iSOi = iSO + iAOi
                                 nijkl = nijkl + 1
*
*---------------------------V_k(ij)*V_k(kl)
*
                                 Indi=Max(iSOi,jSOj)
                                 Indj=iSOi+jSOj-Indi
                                 Indk=Max(kSOk,lSOl)
                                 Indl=kSOk+lSOl-Indk
                                 Indij=(Indi-1)*Indi/2+Indj
                                 Indkl=(Indk-1)*Indk/2+Indl
                                 temp=V_k(Indij)*V_k(Indkl)*CoulFac
*-----Active space contribution (any factor?)
                                 ijVec=mn2K(Indij,1)
                                 klVec=mn2K(Indkl,1)
                                 If (ijVec.eq.0 .or. klVec.eq.0) GoTo 11
                                 Do jp=1,nnP1
                                    temp = temp
     &                                 + Z_p_K(jp,ijVec)*Z_p_K(jp,klVec)
                                 End Do
*
 11                              Continue
                                 PMax=Max(PMax,Abs(temp))
                                 PAO(nijkl,iPAO) = Fac * temp
*
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
      Else If(iMP2prpt .ne. 2) Then
         NumIK = nIJ1(iSym,lSym,iSO)
         If(NumIK.eq.0) Return
         ip_CijK2 = ip_CijK + NumIK

         Do i1 = 1, iCmp(1)
            Do i2 = 1, iCmp(2)
               iSO = iAOtSO(iAO(1)+i1,kOp(1))+iAOst(1)
               jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)

               Do i3 = 1, iCmp(3)
                  Do i4 = 1, iCmp(4)
                     kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
                     lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

                     iPAO = iPAO + 1
                     nijkl = 0

                     Do lAOl = 0, lBas-1
                        lSOl = lSO + lAOl
                        Do kAOk = 0, kBas-1
                           kSOk = kSO + kAOk
                           Indk=Max(kSOk,lSOl)
                           Indl=kSOk+lSOl-Indk
                           Indkl=(Indk-1)*Indk/2+Indl
                           klVec=mn2K(Indkl,1)
                           If(klvec.ne.0) Then
                              iAdrL = NumIK*(klVec-1)
     &                              + iAdrCVec(jSym,iSym,iSO)
                              Call dDaFile(LuCVector(jSym,iSO),2,
     &                                     Work(ip_CijK),
     &                                     NumIK,iAdrL)
                           End If


                           Do jAOj = 0, jBas-1
                              jSOj = jSO + jAOj
                              Do iAOi = 0, iBas-1
                                 iSOi = iSO + iAOi
                                 nijkl = nijkl + 1
*
                                 Indi=Max(iSOi,jSOj)
                                 Indj=iSOi+jSOj-Indi
                                 Indij=(Indi-1)*Indi/2+Indj
                                 ijVec=mn2K(Indij,1)

                                 If(ijVec.ne.klVec.and.ijvec.ne.0) Then
                                    iAdrJ = NumIK*(ijVec-1) +
     &                                      iAdrCVec(jSym,iSym,iSO)
                                    Call dDaFile(LuCVector(jSym,iSO),2,
     &                                           Work(ip_CijK2),NumIK,
     &                                           iAdrJ)
                                    ip_V2 = ip_CijK2
                                 Else
                                    ip_V2 = ip_CijK
                                 End If

                                 temp=V_k(Indij)*V_k(Indkl)*CoulFac

                                 If (ijVec.eq.0 .or. klVec.eq.0) GoTo 22
                                 If(Indi .eq. Indj) Then
                                    Fac_ij = 1.0d0
                                 Else
                                    Fac_ij = 0.5d0
                                 End If
                                 If(Indk .eq. Indl) Then
                                    Fac_kl = 1.0d0
                                 Else
                                    Fac_kl = 0.5d0
                                 End If

*----- Exchange contribution

                                 temp = temp - ExFac*Fac_ij*Fac_kl*
     &                                  dDot_(NumIK,Work(ip_CijK),1,
     &                                             Work(ip_V2),1)
*-----Active space contribution (any factor?)
                                 Do jp=1,nnP1
                                    temp = temp
     &                                + Z_p_K(jp,ijVec)*Z_p_K(jp,klVec)
                                 End Do
*
 22                              Continue
                                 PMax=Max(PMax,Abs(temp))
                                 PAO(nijkl,iPAO) = Fac * temp
*
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
      Else
         NumIK = nIJ1(iSym,lSym,iSO)
         If(NumIK.eq.0) Return
         ip_CijK2 = ip_CijK + NumIK

         Do i1 = 1, iCmp(1)
            Do i2 = 1, iCmp(2)
               iSO = iAOtSO(iAO(1)+i1,kOp(1))+iAOst(1)
               jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
               Do i3 = 1, iCmp(3)
                  Do i4 = 1, iCmp(4)
                     kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
                     lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

                     iPAO = iPAO + 1
                     nijkl = 0
*
                     Do lAOl = 0, lBas-1
                        lSOl = lSO + lAOl
                        Do kAOk = 0, kBas-1
                           kSOk = kSO + kAOk
                           Indk=Max(kSOk,lSOl)
                           Indl=kSOk+lSOl-Indk
                           Indkl=(Indk-1)*Indk/2+Indl
                           klVec=mn2K(Indkl,1)
                           If(klvec.ne.0) Then
                              iAdrL = NumIK*(klVec-1)
     &                              + iAdrCVec(jSym,iSym,iSO)
                              Call dDaFile(LuCVector(jSym,iSO),2,
     &                                     Work(ip_CijK),
     &                                     NumIK,iAdrL)
                           End If


                           Do jAOj = 0, jBas-1
                              jSOj = jSO + jAOj
                              Do iAOi = 0, iBas-1
                                 iSOi = iSO + iAOi
                                 nijkl = nijkl + 1
*
                                 Indi=Max(iSOi,jSOj)
                                 Indj=iSOi+jSOj-Indi
                                 Indij=(Indi-1)*Indi/2+Indj
                                 ijVec=mn2K(Indij,1)

                                 If(ijVec.ne.klVec.and.ijvec.ne.0) Then
                                    iAdrJ = NumIK*(ijVec-1) +
     &                                      iAdrCVec(jSym,iSym,iSO)
                                    Call dDaFile(LuCVector(jSym,iSO),2,
     &                                           Work(ip_CijK2),NumIK,
     &                                           iAdrJ)
                                    ip_V2 = ip_CijK2
                                 Else
                                    ip_V2 = ip_CijK
                                 End If

                                 temp= V_k(Indij)*V_k(Indkl)*CoulFac
     &                               + V_K(Indij)*U_K(Indkl)*CoulFac
     &                               + U_K(Indij)*V_K(Indkl)*CoulFac

                                 If (ijVec.eq.0 .or. klVec.eq.0) GoTo 33
                                 If(Indi .eq. Indj) Then
                                    Fac_ij = 1.0d0
                                 Else
                                    Fac_ij = 0.5d0
                                 End If
                                 If(Indk .eq. Indl) Then
                                    Fac_kl = 1.0d0
                                 Else
                                    Fac_kl = 0.5d0
                                 End If

*
*----- MP2 contribution
                                 Call Compute_A_jk_Mp2(1,ijVec,klVec,
     &                                                 tempJ_mp2,
     &                                                 Fac_ij,Fac_kl,
     &                                                 nAuxVe,2)
                                 temp = temp + tempJ_mp2*CoulFac
*
*----- Exchange contribution


                                 tempK = 2.0d0*Fac_ij*Fac_kl*
     &                                   dDot_(NumIK,Work(ip_CijK),1,
     &                                              Work(ip_V2),1)
                                 Call compute_A_jk_Mp2(1,ijVec,klVec,
     &                                tempK_mp2,
     &                                fac_ij,fac_kl,nAuxVe,1)

                                 tempK = tempK + tempK_mp2

                                 temp = temp - ExFac*tempK*Half
*----- Active space contribution (any factor?)
                                 Do jp=1,nnP1
                                    temp = temp
     &                                 + Z_p_K(jp,ijVec)*Z_p_K(jp,klVec)
                                 End Do
*
 33                              Continue
                                 PMax=Max(PMax,Abs(temp))
                                 PAO(nijkl,iPAO) = Fac * temp
*
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End If
*

      If (iPAO.ne.nPAO) Then
         Write (6,*) ' Error in PGet1_CD2!'
         Call Abend
      End If
*
#ifdef _DEBUG_
      Call RecPrt(' In PGet1_CD2:PAO ',' ',PAO,ijkl,nPAO)
      Call GetMem(' Exit PGet1_CD2','CHECK','REAL',iDum,iDum)
      Call qExit('PGet1_CD2')
#endif

      Call CWTime(Cpu2,Wall2)
      Cpu = Cpu2 - Cpu1
      Wall = Wall2 - Wall1
      tavec(1) = tavec(1) + Cpu
      tavec(2) = tavec(2) + Wall

      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_logical(Shijij)
      End If
      End
