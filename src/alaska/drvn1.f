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
* Copyright (C) 1991,1995, Roland Lindh                                *
************************************************************************
      SubRoutine DrvN1(Grad,Temp,nGrad)
************************************************************************
*                                                                      *
* Object: to compute the molecular gradient contribution due to the    *
*         nuclear repulsion energy.                                    *
*                                                                      *
* Called from: Alaska                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             October '91                                              *
*                                                                      *
*             Modified for ECP's and external electric fields, May '95 *
************************************************************************
      use PCM_arrays, only: PCM_SQ, PCMTess, MM
      use External_Centers
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
      Real*8 A(3), B(3), RB(3), Grad(nGrad), Temp(nGrad), DA(3),
     &       Tempd(3)
      Integer iDCRR(0:7), jCoSet(8,8), iStb(0:7)
      Logical EQ, TstFnc, NoLoop
      Character Lab*80
*
*     Statement function for Cartesian index
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

      iRout = 33
      iPrint = nPrint(iRout)
*     Call qEnter('DrvN1')
*
      iIrrep = 0
*
************************************************************************
*                                                                      *
*            Compute the nuclear repulsion contribution                *
*                                                                      *
************************************************************************
*
      call dcopy_(nGrad,[Zero],0,Temp,1)
*
      mdc = 0
*-----Loop over centers with the same charge
      Do iCnttp = 1, nCnttp
         If (FragCnttp(iCnttp)) Then
           ZA = FragCharge(iCnttp)
         Else
           ZA = Charge(iCnttp)
         End If
         If (ZA.eq.Zero) Go To 101
         ixyz = ipCntr(iCnttp)
*--------Loop over all unique centers of this group
         Do iCnt = 1, nCntr(iCnttp)
            A(1) = Work(ixyz+(iCnt-1)*3)
            A(2) = Work(ixyz+(iCnt-1)*3+1)
            A(3) = Work(ixyz+(iCnt-1)*3+2)
*
            ndc = 0
            Do jCnttp = 1, iCnttp
               If (FragCnttp(jCnttp)) Then
                 ZB = FragCharge(jCnttp)
               Else
                 ZB = Charge(jCnttp)
               End If
               If (ZB.eq.Zero) Go To 201
               If (pChrg(iCnttp).and.pChrg(jCnttp)) Go To 201
               If (FragCnttp(iCnttp).and.FragCnttp(jCnttp)) Go To 201
               ZAZB = ZA * ZB
               jxyz = ipCntr(jCnttp)
               jCntMx = nCntr(jCnttp)
               If (iCnttp.eq.jCnttp) jCntMx = iCnt
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
*
                     If (.Not.pChrg(jCnttp)) Then
                     nDisp = IndDsp(ndc+jCnt,iIrrep)
                     igv=nIrrep/nStab(ndc+jCnt)
                     Do iCar = 0, 2
                        dr_dB=-(A(iCar+1)-RB(iCar+1))/r12
                        iComp = 2**iCar
                        If ( TstFnc(iOper,nIrrep,
     &                     iCoSet(0,0,ndc+jCnt),
     &                     nIrrep/nStab(ndc+jCnt),iChTbl,iIrrep,
     &                     iComp,nStab(ndc+jCnt)) ) Then
                           nDisp = nDisp + 1
                           If (Direct(nDisp)) Then
                              ps = DBLE(iPrmt(nOp,iChBas(2+iCar)))
                              Temp(nDisp) = Temp(nDisp) +
     &                           ps * One/DBLE(igv) * PreFct *
     &                           dr_dB * df_dr
                           End If
                        End If
                     End Do
                     End If
 301                 Continue
                  End Do
*
               End Do
 201           Continue
               ndc = ndc + nCntr(jCnttp)
            End Do
         End Do
 101     Continue
         mdc = mdc + nCntr(iCnttp)
      End Do
      If (iPrint.ge.15) Then
         Lab=' The Nuclear Repulsion Contribution'
         Call PrGrad(Lab,Temp,nGrad,lIrrep,ChDisp,5)
      End If
*
      Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
************************************************************************
*                                                                      *
*           Compute contribution due to the external field             *
*                                                                      *
************************************************************************
*

      If (.Not.lXF) Go To 666
*
      If (nIrrep.eq.8) Then
         nOper=3
      Else If (nIrrep.eq.4) Then
         nOper=2
      Else If (nIrrep.eq.2) Then
         nOper=1
      Else
         nOper=0
      End If
*
      If((nOrd_XF.gt.1).or.(iXPolType.gt.0)) Then
         Call WarningMessage(2,'Error in DrvN1')
         Write(6,*)'Sorry, gradients are not implemented for'
         Write(6,*)'higher XF than dipoles or for polarisabilities'
         Call Quit_OnUserError()
      EndIf

      call dcopy_(nGrad,[Zero],0,Temp,1)
*
      iDum=0
      Do iFd = 1, nXF
         ZA   = XF(4,iFd)
         If(nOrd_XF.eq.0) Then
            DA(1:3)=Zero
         Else
            DA(1:3)= XF(5:7,iFd)
         EndIf
         NoLoop = ZA.eq.Zero .and. DA(1).eq.Zero .and. DA(2).eq.Zero
     &            .and. DA(3).eq.Zero
         If (NoLoop) Go To 102
         A(1:3)=XF(1:3,iFd)
         iChxyz=iChAtm(A,iOper,nOper,iChBas(2))
         Call Stblz(iChxyz,iOper,nIrrep,nStb,iStb,iDum,jCoSet)
*
         ndc = 0
         Do jCnttp = 1, nCnttp
            ZB = Charge(jCnttp)
            If (ZB.eq.Zero) Go To 202
            If (pChrg(jCnttp)) Go To 202
            If (FragCnttp(jCnttp)) Go To 202
            ZAZB = ZA * ZB
            jxyz = ipCntr(jCnttp)
            Do jCnt = 1, nCntr(jCnttp)
               B(1) = Work(jxyz+(jCnt-1)*3  )
               B(2) = Work(jxyz+(jCnt-1)*3+1)
               B(3) = Work(jxyz+(jCnt-1)*3+2)
*
*              Find the DCR for the two centers
*
               Call DCR(LmbdR,iOper,nIrrep,
     &                  iStb,nStb,
     &                  jStab(0,ndc+jCnt),nStab(ndc+jCnt),
     &                  iDCRR,nDCRR)
*
               PreFct = DBLE(nIrrep)/DBLE(LmbdR)
               Do iR = 0, nDCRR-1
                  RB(1) = DBLE(iPhase(1,iDCRR(iR)))*B(1)
                  RB(2) = DBLE(iPhase(2,iDCRR(iR)))*B(2)
                  RB(3) = DBLE(iPhase(3,iDCRR(iR)))*B(3)
                  nOp = NrOpr(iDCRR(iR),iOper,nIrrep)
                  If (EQ(A,RB)) Go To 302
                  r12 = Sqrt((A(1)-RB(1))**2 +
     &                       (A(2)-RB(2))**2 +
     &                       (A(3)-RB(3))**2 )
                  DARB=DA(1)*(A(1)-RB(1))+
     &                 DA(2)*(A(2)-RB(2))+
     &                 DA(3)*(A(3)-RB(3))

*
*                 The factor u/g will ensure that the value of the
*                 gradient in symmetry adapted and no symmetry basis
*                 will have the same value.
*
                  fab0=One
                  fab1=One
                  fab2=Three
                  If (ECP(jCnttp)) Then
*--------------------Add contribution from M1 operator
                     Cnt0M1=Zero
                     Cnt1M1=Zero
                     Do iM1xp=0, nM1(jCnttp)-1
                       Gamma = Work(ipM1xp(jCnttp)+iM1xp)
                       CffM1 = Work(ipM1cf(jCnttp)+iM1xp)
                       Cnt0M1= Cnt0M1+(CffM1*Exp(-Gamma*r12**2))
                       Cnt1M1= Cnt1M1+Gamma*(CffM1*Exp(-Gamma*r12**2))
                     End Do
                     fab0=fab0+Cnt0M1-Two*r12**2*Cnt1M1
                     fab1=fab1+Cnt0M1
                     fab2=fab2+Three*Cnt0M1+Two*r12**2*Cnt1M1
*--------------------Add contribution from M2 operator
                     Cnt0M2=Zero
                     Cnt1M2=Zero
                     Do iM2xp=0, nM2(jCnttp)-1
                       Gamma = Work(ipM2xp(jCnttp)+iM2xp)
                       CffM2 = Work(ipM2cf(jCnttp)+iM2xp)
                       Cnt0M2= Cnt0M2+(CffM2*Exp(-Gamma*r12**2))
                       Cnt1M2= Cnt1M2+Gamma*(CffM2*Exp(-Gamma*r12**2))
                     End Do
                     fab0=fab0+Two*(r12*Cnt0M2-r12**3*Cnt1M2)
                     fab1=fab1+r12*Cnt0M2
                     fab2=fab2-Two*(r12*Cnt0M2-r12**3*Cnt1M2)
                  End If
*
                  nDisp = IndDsp(ndc+jCnt,iIrrep)
                  igv=nIrrep/nStab(ndc+jCnt)
                  Do iCar = 0, 2
                     iComp = 2**iCar
                     If ( TstFnc(iOper,nIrrep,
     &                  iCoSet(0,0,ndc+jCnt),
     &                  nIrrep/nStab(ndc+jCnt),iChTbl,iIrrep,
     &                  iComp,nStab(ndc+jCnt)) ) Then
                        nDisp = nDisp + 1
                        If (Direct(nDisp)) Then
                           ps = DBLE(iPrmt(nOp,iChBas(2+iCar)))
                           Temp(nDisp) = Temp(nDisp) +
     &                        ps * One/DBLE(igv) * PreFct * (
     &                        ZAZB*fab0*(A(iCar+1)-RB(iCar+1))/(r12**3)
     &                       +ZB *(fab1*DA(iCar+1)            /(r12**3)
     &                       -DARB*fab2*(A(iCar+1)-RB(iCar+1))/(r12**5))
     &                        )
                        End If
                     End If
                  End Do   ! End loop over cartesian components, iCar
*
 302              Continue
               End Do      ! End loop over DCR operators, iR
*
            End Do         ! End over centers, jCnt
 202        Continue
            ndc = ndc + nCntr(jCnttp)
         End Do            ! End over basis set types, jCnttp
 102     Continue
      End Do               ! End of centers of the external field, iFD
      If (iPrint.ge.15) Then
         Lab=' The Nuclear External Electric Field Contribution'
         Call PrGrad(Lab,Temp,nGrad,lIrrep,ChDisp,5)
      End If
*
      Call DaXpY_(nGrad,One,Temp,1,Grad,1)
 666  Continue
*
************************************************************************
*                                                                      *
*          Compute contributions due to the reaction field             *
*                 KirkWood Model                                       *
*                                                                      *
************************************************************************
*
      If (lRF.and..Not.lLangevin.and..Not.PCM) Then
         nCav=(lMax+1)*(lMax+2)*(lMax+3)/6
*
*------- Get the multipole moments
*
         Call Get_dArray('RCTFLD',MM,nCav*2)
         If (iPrint.ge.99) Call RecPrt('Total Multipole Moments',' ',
     &                                 MM(1,1),1,nCav)
         If (iPrint.ge.99) Call RecPrt('Total Electric Field',
     &                                 ' ',MM(1,2),1,nCav)
*
      call dcopy_(nGrad,[Zero],0,Temp,1)
*
      ip = 0
      Do ir = 0, lMax
         Do ix = ir, 0, -1
            Do iy = ir-ix, 0, -1
               iz = ir-ix-iy
               ip = ip + 1
               If (iPrint.ge.99) Write (6,*) ' ix,iy,iz=',ix,iy,iz
*
               mdc = 0
               Do iCnttp = 1, nCnttp
                  If (Charge(iCnttp).eq.Zero) Go To 103
                  If (FragCnttp(iCnttp)) Go To 103
                  ZA = Charge(iCnttp)
                  ixyz = ipCntr(iCnttp)
                  If (iPrint.ge.99) Then
                     Write (6,*) ' Charge=',ZA
                     Write (6,*) ' ixyz=',ixyz
                     Call RecPrt(' Centers',' ',Work(ixyz),3,
     &                            nCntr(iCnttp))
                  End If
                  Do iCnt = 1, nCntr(iCnttp)
                     A(1) = Work(ipCntr(iCnttp)+(iCnt-1)*3)
                     A(2) = Work(ipCntr(iCnttp)+(iCnt-1)*3+1)
                     A(3) = Work(ipCntr(iCnttp)+(iCnt-1)*3+2)

                     If (ix.eq.0) Then
                        CCoMx =One
                        CCoMxd=Zero
                     Else If (ix.eq.1) Then
                        CCoMx =A(1)
                        CCoMxd=One
                     Else
                        CCoMx =A(1)**ix
                        CCoMxd=DBLE(ix)*A(1)**(ix-1)
                     End If
*
                     If (iy.eq.0) Then
                        CCoMy =One
                        CCoMyd=Zero
                     Else If (iy.eq.1) Then
                        CCoMy =A(2)
                        CCoMyd=One
                     Else
                        CCoMy =A(2)**iy
                        CCoMyd=DBLE(iy)*A(2)**(iy-1)
                     End If
*
                     If (iz.eq.0) Then
                        CCoMz =One
                        CCoMzd=Zero
                     Else If (iz.eq.1) Then
                        CCoMz =A(3)
                        CCoMzd=One
                     Else
                        CCoMz =A(3)**iz
                        CCoMzd=DBLE(iz)*A(3)**(iz-1)
                     End If
                     tempd(1)= MM(ip,2) *ZA * CCoMxd* CCoMy * CCoMz
                     tempd(2)= MM(ip,2) *ZA * CCoMx * CCoMyd* CCoMz
                     tempd(3)= MM(ip,2) *ZA * CCoMx * CCoMy * CCoMzd
                     If (iPrint.ge.99) Then
                        Write (6,*) CCoMx, CCoMy, CCoMz
                        Write (6,*) 'Work(ip)=',Work(ip)
                        Write (6,*) 'tempd=',tempd
                     End If
*
*                    Distribute gradient
*
                     nDisp=IndDsp(mdc+iCnt,iIrrep)
                     Do iCar = 0, 2
                        iComp = 2**iCar
                        If ( TstFnc(iOper,nIrrep,
     &                     iCoSet(0,0,mdc+iCnt),
     &                     nIrrep/nStab(mdc+iCnt),iChTbl,iIrrep,
     &                     iComp,nStab(mdc+iCnt)) ) Then
                           nDisp = nDisp + 1
                           If (Direct(nDisp)) Then
                              Temp(nDisp) = Temp(nDisp) - Tempd(iCar+1)
                           End If
                        End If
                     End Do

*
                  End Do
 103              Continue
                  mdc = mdc + nCntr(iCnttp)
               End Do
*
            End Do
         End Do
      End Do
      If (iPrint.ge.15) Then
         Lab=' The Nuclear Reaction Field (KirkWood) Contribution'
         Call PrGrad(Lab,Temp,nGrad,lIrrep,ChDisp,5)
      End If
*
      Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
      Else If(lRF.and.PCM) Then
*
************************************************************************
*                                                                      *
*          Compute contributions due to the reaction field             *
*                      PCM Model                                       *
*                                                                      *
************************************************************************
*
      call dcopy_(nGrad,[Zero],0,Temp,1)
*
*---- Loop over tiles
*
*
      Do iTs = 1, nTs
         ZA   = PCM_SQ(1,iTs)+PCM_SQ(2,iTS)
         NoLoop = ZA.eq.Zero
         ZA = ZA / DBLE(nIrrep)
         If (NoLoop) Go To 112
         A(1:3) = PCMTess(1:3,iTs)
*
*------- Tile only stabilized by the unit operator
*
         nStb=1
         iStb(0)=0
*
         ndc = 0
         Do jCnttp = 1, nCnttp
            ZB = Charge(jCnttp)
            If (ZB.eq.Zero) Go To 212
            If (pChrg(jCnttp)) Go To 212
            If (FragCnttp(jCnttp)) Go To 212
            ZAZB = ZA * ZB
            jxyz = ipCntr(jCnttp)
            Do jCnt = 1, nCntr(jCnttp)
               B(1) = Work(jxyz+(jCnt-1)*3  )
               B(2) = Work(jxyz+(jCnt-1)*3+1)
               B(3) = Work(jxyz+(jCnt-1)*3+2)
*
*              Find the DCR for the two centers
*
               Call DCR(LmbdR,iOper,nIrrep,
     &                  iStb,nStb,
     &                  jStab(0,ndc+jCnt),nStab(ndc+jCnt),
     &                  iDCRR,nDCRR)
*
               PreFct = ZAZB*DBLE(nIrrep)/DBLE(LmbdR)
               Do iR = 0, nDCRR-1
                  RB(1) = DBLE(iPhase(1,iDCRR(iR)))*B(1)
                  RB(2) = DBLE(iPhase(2,iDCRR(iR)))*B(2)
                  RB(3) = DBLE(iPhase(3,iDCRR(iR)))*B(3)
                  nOp = NrOpr(iDCRR(iR),iOper,nIrrep)
                  If (EQ(A,RB)) Go To 312
                  r12 = Sqrt((A(1)-RB(1))**2 +
     &                       (A(2)-RB(2))**2 +
     &                       (A(3)-RB(3))**2 )
*
*                 The factor u/g will ensure that the value of the
*                 gradient in symmetry adapted and no symmetry basis
*                 will have the same value.
*
                  fab=One
                  dfab=Zero
                  If (ECP(jCnttp)) Then
*--------------------Add contribution from M1 operator
                     Cnt0M1=Zero
                     Cnt1M1=Zero
                     Do iM1xp=0, nM1(jCnttp)-1
                       Gamma = Work(ipM1xp(jCnttp)+iM1xp)
                       CffM1 = Work(ipM1cf(jCnttp)+iM1xp)
                       Cnt0M1= Cnt0M1+(CffM1*Exp(-Gamma*r12**2))
                       Cnt1M1= Cnt1M1+Gamma*(CffM1*Exp(-Gamma*r12**2))
                     End Do
                     fab=fab+Cnt0M1
                     dfab=dfab-Two*r12*Cnt1M1
*--------------------Add contribution from M2 operator
                     Cnt0M2=Zero
                     Cnt1M2=Zero
                     Do iM2xp=0, nM2(jCnttp)-1
                       Gamma = Work(ipM2xp(jCnttp)+iM2xp)
                       CffM2 = Work(ipM2cf(jCnttp)+iM2xp)
                       Cnt0M2= Cnt0M2+(CffM2*Exp(-Gamma*r12**2))
                       Cnt1M2= Cnt1M2+Gamma*(CffM2*Exp(-Gamma*r12**2))
                     End Do
                     fab=fab+r12*Cnt0M2
                     dfab=dfab+(Cnt0M2-Two*r12**2*Cnt1M2)
                  End If
                  df_dr=(dfab*r12-fab)/r12**2
*
                  nDisp = IndDsp(ndc+jCnt,iIrrep)
                  igv=nIrrep/nStab(ndc+jCnt)
                  Do iCar = 0, 2
                     dr_dB=-(A(iCar+1)-RB(iCar+1))/r12
                     iComp = 2**iCar
                     If ( TstFnc(iOper,nIrrep,
     &                  iCoSet(0,0,ndc+jCnt),
     &                  nIrrep/nStab(ndc+jCnt),iChTbl,iIrrep,
     &                  iComp,nStab(ndc+jCnt)) ) Then
                        nDisp = nDisp + 1
                        If (Direct(nDisp)) Then
                           ps = DBLE(iPrmt(nOp,iChBas(2+iCar)))
                           Temp(nDisp) = Temp(nDisp) +
     &                        ps * One/DBLE(igv) * PreFct *
     &                        dr_dB * df_dr
                        End If
                     End If
                  End Do   ! End loop over cartesian components, iCar
*
 312              Continue
               End Do      ! End loop over DCR operators, iR
*
            End Do         ! End over centers, jCnt
 212        Continue
            ndc = ndc + nCntr(jCnttp)
         End Do            ! End over basis set types, jCnttp
 112     Continue
      End Do               ! End of tiles
*
      If (iPrint.ge.15) Then
         Lab=' The Nuclear Reaction Field (PCM) Contribution'
         Call PrGrad(Lab,Temp,nGrad,lIrrep,ChDisp,5)
      End If
*
      Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
*---- Add contribution due to the tiles.
*
      Call FZero(Temp,nGrad)
      Call PCM_Cav_grd(Temp,nGrad)
      If (iPrint.ge.15) Then
         Lab=' The Cavity PCM Contribution'
         Call PrGrad(Lab,Temp,nGrad,lIrrep,ChDisp,5)
      End If
      Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
*---- Add contribution due to the electric field on the tiles.
*
      If (Conductor) Then
         Call FZero(Temp,nGrad)
         Call PCM_EF_grd(Temp,nGrad)
         If (iPrint.ge.15) Then
            Lab=' The EF PCM Contribution'
            Call PrGrad(Lab,Temp,nGrad,lIrrep,ChDisp,5)
         End If
         Call DaXpY_(nGrad,-One,Temp,1,Grad,1)
      End If
*
      End If
*
*     Call qExit('DrvN1')
      Return
      End
