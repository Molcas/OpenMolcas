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
* Copyright (C) 2010, Francesco Aquilante                              *
************************************************************************
      SubRoutine DrvN1_EMB(Grad,Temp,nGrad)
************************************************************************
*                                                                      *
* Object: to compute the molecular gradient contribution due to the    *
*         inter-subsystem nuclear repulsion energy                     *
*                                                                      *
* Called from: Alaska                                                  *
*                                                                      *
* Author : F. Aquilante, Geneva, Nov 2010                              *
************************************************************************
      Use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "SysDef.fh"
#include "print.fh"
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "disp.fh"
#include "rctfld.fh"
#include "WrkSpc.fh"
#include "status.fh"
      Real*8 A(3), B(3), RB(3), Grad(nGrad), Temp(nGrad)
      Integer iDCRR(0:7)
      Logical EQ, TstFnc
      Character*16 NamRfil
      Character*80 Lab
*
************************************************************************
*     Statement function for Charges of subsystem B
      Charge_B(i) = Work(ip_ChargeB+i-1)
************************************************************************
*
      If (nIrrep .gt. 1) Then
         Call WarningMessage(2,'Error in DrvN1_Emb')
         Write(6,*)'Sorry, OFE gradient code does not understand'
         Write(6,*)'the use of subsystem symmetry!'
         Call Quit_OnUserError()
      EndIf
      iRout = 33
      iPrint = nPrint(iRout)
*
      iIrrep = 0
      call dcopy_(nGrad,[Zero],0,Temp,1)
*
      Call Get_NameRun(NamRfil) ! save the old RUNFILE name
      Call NameRun('AUXRFIL')   ! switch RUNFILE name
*
      Call GetMem('B-Charges','Allo','Real',ip_ChargeB,nCnttp)
      Call Get_dArray('Nuclear charge',Work(ip_ChargeB),nCnttp)
*
      Call NameRun(NamRfil)   ! switch back to old RUNFILE name
*
      ZA=1.0d0
      iCnttp=1
      Do while (iCnttp.le.nCnttp .and. ZA.gt.0.0d0)
         ZA = Charge(iCnttp)
         iCnttp = iCnttp + 1
      End Do
      iCnttp_B=iCnttp-1  ! start of atoms of subsystem B
*
      If (iCnttp_B.eq.1) Then ! subsystem B comes first
         ZB=0.0d0
         nCnttp_B=1
         Do while (nCnttp_B.le.nCnttp .and. ZB.eq.0.0d0)
            ZB = Charge(nCnttp_B)
            nCnttp_B = nCnttp_B + 1
         End Do
         nCnttp_B=nCnttp_B-1  ! end of atoms of subsystem B
      Else
         nCnttp_B=nCnttp
      EndIf
*
      mdc = 0
*-----Loop over centers with the same charge (A-subsystem)
      Do iCnttp = 1, nCnttp
         ZA = Charge(iCnttp)
         If (iCnttp.ge.iCnttp_B .and. iCnttp.le.nCnttp_B
     &                          .and. ZA.gt.0.0d0) Then
            Call WarningMessage(2,'Internal error in DrvN1_Emb')
            Write(6,*)' Subsystems must come one after the other'
            Call Abend()
         EndIf
         If (ZA.eq.Zero) Go To 101
         ixyz = dbsc(iCnttp)%ipCntr
*--------Loop over all unique centers of this group (A-subsystem)
         Do iCnt = 1, nCntr(iCnttp)
            A(1) = Work(ixyz+(iCnt-1)*3)
            A(2) = Work(ixyz+(iCnt-1)*3+1)
            A(3) = Work(ixyz+(iCnt-1)*3+2)
*
            ndc = 0
            Do jCnttp = iCnttp_B, nCnttp_B  ! (B-subsystem)

               ZB = Charge_B(jCnttp)

               If (ZB.eq.Zero) Go To 201
               ZAZB = ZA * ZB
               jxyz = dbsc(jCnttp)%ipCntr
               jCntMx = nCntr(jCnttp)
               Do jCnt = 1, jCntMx
                  B(1) = Work(jxyz+(jCnt-1)*3  )
                  B(2) = Work(jxyz+(jCnt-1)*3+1)
                  B(3) = Work(jxyz+(jCnt-1)*3+2)
*
                  Fact = One
*                 Factor due to resticted summation
                  If (EQ(A,B)) Fact = Half
*
*                 Find the DCR for the two centers
*
                  Call DCR(LmbdR,iOper,nIrrep,
     &                     jStab(0,mdc+iCnt),nStab(mdc+iCnt),
     &                     jStab(0,ndc+jCnt),nStab(ndc+jCnt),
     &                     iDCRR,nDCRR)
*
                  PreFct = Fact*ZAZB*DBLE(nIrrep)/DBLE(LmbdR)
                  Do iR = 0, nDCRR-1
                     RB(1) = DBLE(iPhase(1,iDCRR(iR)))*B(1)
                     RB(2) = DBLE(iPhase(2,iDCRR(iR)))*B(2)
                     RB(3) = DBLE(iPhase(3,iDCRR(iR)))*B(3)
                     nOp = NrOpr(iDCRR(iR),iOper,nIrrep)
                     If (EQ(A,RB)) Go To 301
                     r12 = Sqrt((A(1)-RB(1))**2 +
     &                          (A(2)-RB(2))**2 +
     &                          (A(3)-RB(3))**2 )
*
*                    The factor u/g will ensure that the value of the
*                    gradient in symmetry adapted and no symmetry basis
*                    will have the same value.
*
                     fab=One
                     dfab=Zero
                     If (ECP(iCnttp)) Then
*-----------------------Add contribution from M1 operator
                        Cnt0M1=Zero
                        Cnt1M1=Zero
                        Do iM1xp=0, nM1(iCnttp)-1
                          Gamma =Work(ipM1xp(iCnttp)+iM1xp)
                          CffM1 =Work(ipM1cf(iCnttp)+iM1xp)
                          Cnt0M1=Cnt0M1+(CffM1*Exp(-Gamma*r12**2))
                          Cnt1M1=Cnt1M1+Gamma*(CffM1*Exp(-Gamma*r12**2))
                        End Do
                        fab=fab+Cnt0M1
                        dfab=dfab-Two*r12*Cnt1M1
*-----------------------Add contribution from M2 operator
                        Cnt0M2=Zero
                        Cnt1M2=Zero
                        Do iM2xp=0, nM2(iCnttp)-1
                          Gamma =Work(ipM2xp(iCnttp)+iM2xp)
                          CffM2 =Work(ipM2cf(iCnttp)+iM2xp)
                          Cnt0M2=Cnt0M2+(CffM2*Exp(-Gamma*r12**2))
                          Cnt1M2=Cnt1M2+Gamma*(CffM2*Exp(-Gamma*r12**2))
                        End Do
                        fab=fab+r12*Cnt0M2
                        dfab=dfab+(Cnt0M2-Two*r12**2*Cnt1M2)
                     End If
                     If (ECP(jCnttp)) Then
*-----------------------Add contribution from M1 operator
                        Cnt0M1=Zero
                        Cnt1M1=Zero
                        Do iM1xp=0, nM1(jCnttp)-1
                          Gamma =Work(ipM1xp(jCnttp)+iM1xp)
                          CffM1 =Work(ipM1cf(jCnttp)+iM1xp)
                          Cnt0M1=Cnt0M1+(CffM1*Exp(-Gamma*r12**2))
                          Cnt1M1=Cnt1M1+Gamma*(CffM1*Exp(-Gamma*r12**2))
                        End Do
                        fab=fab+Cnt0M1
                        dfab=dfab-Two*r12*Cnt1M1
*-----------------------Add contribution from M2 operator
                        Cnt0M2=Zero
                        Cnt1M2=Zero
                        Do iM2xp=0, nM2(jCnttp)-1
                          Gamma =Work(ipM2xp(jCnttp)+iM2xp)
                          CffM2 =Work(ipM2cf(jCnttp)+iM2xp)
                          Cnt0M2=Cnt0M2+(CffM2*Exp(-Gamma*r12**2))
                          Cnt1M2=Cnt1M2+Gamma*(CffM2*Exp(-Gamma*r12**2))
                        End Do
                        fab=fab+r12*Cnt0M2
                        dfab=dfab+(Cnt0M2-Two*r12**2*Cnt1M2)
                     End If
                     df_dr=(dfab*r12-fab)/r12**2
*
                     If (.Not.pChrg(iCnttp)) Then
                     nDisp = IndDsp(mdc+iCnt,iIrrep)
                     igu=nIrrep/nStab(mdc+iCnt)
                     Do iCar = 0, 2
                        dr_dA=(A(iCar+1)-RB(iCar+1))/r12
                        iComp = 2**iCar
                        If ( TstFnc(iOper,nIrrep,
     &                     iCoSet(0,0,mdc+iCnt),
     &                     nIrrep/nStab(mdc+iCnt),iChTbl,iIrrep,
     &                     iComp,nStab(mdc+iCnt)) ) Then
                           nDisp = nDisp + 1
                           If (Direct(nDisp)) Then
                              Temp(nDisp) = Temp(nDisp) +
     &                           One/DBLE(igu) * PreFct *
     &                           dr_dA * df_dr
                           End If
                        End If
                     End Do
                     End If
 301                 Continue
                  End Do
               End Do
 201           Continue
               ndc = ndc + nCntr(jCnttp)
            End Do
         End Do
 101     Continue
         mdc = mdc + nCntr(iCnttp)
      End Do
      If (iPrint.ge.15) Then
         Lab=' OFE Nuclear Repulsion Contribution'
         Call PrGrad(Lab,Temp,nGrad,lIrrep,ChDisp,5)
      End If
*
      Call GetMem('B-Charges','Free','Real',ip_ChargeB,nCnttp)
*
      Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
      Return
      End
