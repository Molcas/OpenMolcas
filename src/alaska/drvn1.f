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
      use Basis_Info
      use Center_Info
      use PCM_arrays, only: PCM_SQ, PCMTess, MM
      use External_Centers
      use Symmetry_Info, only: iChBas
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
      iRout = 33
      iPrint = nPrint(iRout)
*     iPrint=15
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
         If (dbsc(iCnttp)%Frag) Then
           ZA = dbsc(iCnttp)%FragCharge
         Else
           ZA = dbsc(iCnttp)%Charge
         End If
         If (ZA.eq.Zero) Go To 101
*--------Loop over all unique centers of this group
         Do iCnt = 1, dbsc(iCnttp)%nCntr
            A(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
*
            ndc = 0
            Do jCnttp = 1, iCnttp
               If (dbsc(jCnttp)%Frag) Then
                 ZB = dbsc(jCnttp)%FragCharge
               Else
                 ZB = dbsc(jCnttp)%Charge
               End If
               If (ZB.eq.Zero) Go To 201
               If (dbsc(iCnttp)%pChrg.and.dbsc(jCnttp)%pChrg) Go To 201
               If (dbsc(iCnttp)%Frag.and.dbsc(jCnttp)%Frag) Go To 201
               ZAZB = ZA * ZB
               jCntMx = dbsc(jCnttp)%nCntr
               If (iCnttp.eq.jCnttp) jCntMx = iCnt
               Do jCnt = 1, jCntMx
                  B(1:3)=dbsc(jCnttp)%Coor(1:3,jCnt)
*
                  Fact = One
*                 Factor due to resticted summation
                  If (EQ(A,B)) Fact = Half
*
*                 Find the DCR for the two centers
*
                  Call DCR(LmbdR,
     &                     dc(mdc+iCnt)%iStab,dc(mdc+iCnt)%nStab,
     &                     dc(ndc+jCnt)%iStab,dc(ndc+jCnt)%nStab,
     &                     iDCRR,nDCRR)
*
                  PreFct = Fact*ZAZB*DBLE(nIrrep)/DBLE(LmbdR)
                  Do iR = 0, nDCRR-1
                     Call OA(iDCRR(iR),B,RB)
                     nOp = NrOpr(iDCRR(iR))
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
                     If (dbsc(iCnttp)%ECP) Then
*-----------------------Add contribution from M1 operator
                        Cnt0M1=Zero
                        Cnt1M1=Zero
                        Do iM1xp=1, dbsc(iCnttp)%nM1
                          Gamma =dbsc(iCnttp)%M1xp(iM1xp)
                          CffM1 =dbsc(iCnttp)%M1cf(iM1xp)
                          Cnt0M1=Cnt0M1+(CffM1*Exp(-Gamma*r12**2))
                          Cnt1M1=Cnt1M1+Gamma*(CffM1*Exp(-Gamma*r12**2))
                        End Do
                        fab=fab+Cnt0M1
                        dfab=dfab-Two*r12*Cnt1M1
*-----------------------Add contribution from M2 operator
                        Cnt0M2=Zero
                        Cnt1M2=Zero
                        Do iM2xp=1, dbsc(iCnttp)%nM2
                          Gamma =dbsc(iCnttp)%M2xp(iM2xp)
                          CffM2 =dbsc(iCnttp)%M2cf(iM2xp)
                          Cnt0M2=Cnt0M2+(CffM2*Exp(-Gamma*r12**2))
                          Cnt1M2=Cnt1M2+Gamma*(CffM2*Exp(-Gamma*r12**2))
                        End Do
                        fab=fab+r12*Cnt0M2
                        dfab=dfab+(Cnt0M2-Two*r12**2*Cnt1M2)
                     End If
                     If (dbsc(jCnttp)%ECP) Then
*-----------------------Add contribution from M1 operator
                        Cnt0M1=Zero
                        Cnt1M1=Zero
                        Do iM1xp=1, dbsc(jCnttp)%nM1
                          Gamma =dbsc(jCnttp)%M1xp(iM1xp)
                          CffM1 =dbsc(jCnttp)%M1cf(iM1xp)
                          Cnt0M1=Cnt0M1+(CffM1*Exp(-Gamma*r12**2))
                          Cnt1M1=Cnt1M1+Gamma*(CffM1*Exp(-Gamma*r12**2))
                        End Do
                        fab=fab+Cnt0M1
                        dfab=dfab-Two*r12*Cnt1M1
*-----------------------Add contribution from M2 operator
                        Cnt0M2=Zero
                        Cnt1M2=Zero
                        Do iM2xp=1, dbsc(jCnttp)%nM2
                          Gamma =dbsc(jCnttp)%M2xp(iM2xp)
                          CffM2 =dbsc(jCnttp)%M2cf(iM2xp)
                          Cnt0M2=Cnt0M2+(CffM2*Exp(-Gamma*r12**2))
                          Cnt1M2=Cnt1M2+Gamma*(CffM2*Exp(-Gamma*r12**2))
                        End Do
                        fab=fab+r12*Cnt0M2
                        dfab=dfab+(Cnt0M2-Two*r12**2*Cnt1M2)
                     End If
                     df_dr=(dfab*r12-fab)/r12**2
*
                     If (.Not.dbsc(iCnttp)%pChrg) Then
                     nDisp = IndDsp(mdc+iCnt,iIrrep)
                     igu=nIrrep/dc(mdc+iCnt)%nStab
                     Do iCar = 0, 2
                        dr_dA=(A(iCar+1)-RB(iCar+1))/r12
                        iComp = 2**iCar
                        If ( TstFnc(dc(mdc+iCnt)%iCoSet,
     &                     iIrrep,iComp,dc(mdc+iCnt)%nStab) ) Then
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
                     If (.Not.dbsc(jCnttp)%pChrg) Then
                     nDisp = IndDsp(ndc+jCnt,iIrrep)
                     igv=nIrrep/dc(ndc+jCnt)%nStab
                     Do iCar = 0, 2
                        dr_dB=-(A(iCar+1)-RB(iCar+1))/r12
                        iComp = 2**iCar
                        If ( TstFnc(dc(ndc+jCnt)%iCoSet,
     &                     iIrrep,
     &                     iComp,dc(ndc+jCnt)%nStab) ) Then
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
               ndc = ndc + dbsc(jCnttp)%nCntr
            End Do
         End Do
 101     Continue
         mdc = mdc + dbsc(iCnttp)%nCntr
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
         iChxyz=iChAtm(A)
         Call Stblz(iChxyz,nStb,iStb,iDum,jCoSet)
*
         ndc = 0
         Do jCnttp = 1, nCnttp
            ZB = dbsc(jCnttp)%Charge
            If (ZB.eq.Zero) Go To 202
            If (dbsc(jCnttp)%pChrg) Go To 202
            If (dbsc(jCnttp)%Frag) Go To 202
            ZAZB = ZA * ZB
            Do jCnt = 1, dbsc(jCnttp)%nCntr
               B(1:3)=dbsc(jCnttp)%Coor(1:3,jCnt)
*
*              Find the DCR for the two centers
*
               Call DCR(LmbdR,
     &                  iStb,nStb,
     &                  dc(ndc+jCnt)%iStab,dc(ndc+jCnt)%nStab,
     &                  iDCRR,nDCRR)
*
               PreFct = DBLE(nIrrep)/DBLE(LmbdR)
               Do iR = 0, nDCRR-1
                  Call OA(iDCRR(iR),B,RB)
                  nOp = NrOpr(iDCRR(iR))
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
                  If (dbsc(jCnttp)%ECP) Then
*--------------------Add contribution from M1 operator
                     Cnt0M1=Zero
                     Cnt1M1=Zero
                     Do iM1xp=1, dbsc(jCnttp)%nM1
                       Gamma = dbsc(jCnttp)%M1xp(iM1xp)
                       CffM1 = dbsc(jCnttp)%M1cf(iM1xp)
                       Cnt0M1= Cnt0M1+(CffM1*Exp(-Gamma*r12**2))
                       Cnt1M1= Cnt1M1+Gamma*(CffM1*Exp(-Gamma*r12**2))
                     End Do
                     fab0=fab0+Cnt0M1-Two*r12**2*Cnt1M1
                     fab1=fab1+Cnt0M1
                     fab2=fab2+Three*Cnt0M1+Two*r12**2*Cnt1M1
*--------------------Add contribution from M2 operator
                     Cnt0M2=Zero
                     Cnt1M2=Zero
                     Do iM2xp=1, dbsc(jCnttp)%nM2
                       Gamma = dbsc(jCnttp)%M2xp(iM2xp)
                       CffM2 = dbsc(jCnttp)%M2cf(iM2xp)
                       Cnt0M2= Cnt0M2+(CffM2*Exp(-Gamma*r12**2))
                       Cnt1M2= Cnt1M2+Gamma*(CffM2*Exp(-Gamma*r12**2))
                     End Do
                     fab0=fab0+Two*(r12*Cnt0M2-r12**3*Cnt1M2)
                     fab1=fab1+r12*Cnt0M2
                     fab2=fab2-Two*(r12*Cnt0M2-r12**3*Cnt1M2)
                  End If
*
                  nDisp = IndDsp(ndc+jCnt,iIrrep)
                  igv=nIrrep/dc(ndc+jCnt)%nStab
                  Do iCar = 0, 2
                     iComp = 2**iCar
                     If ( TstFnc(dc(ndc+jCnt)%iCoSet,
     &                  iIrrep,iComp,dc(ndc+jCnt)%nStab) ) Then
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
            ndc = ndc + dbsc(jCnttp)%nCntr
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
                  If (dbsc(iCnttp)%Charge.eq.Zero) Go To 103
                  If (dbsc(iCnttp)%Frag) Go To 103
                  ZA = dbsc(iCnttp)%Charge
                  If (iPrint.ge.99) Then
                     Write (6,*) ' Charge=',ZA
                     Write (6,*) ' ixyz=',ixyz
                     Call RecPrt(' Centers',' ',
     &                           dbsc(iCnttp)%Coor,3,
     &                           dbsc(iCnttp)%nCntr)
                  End If
                  Do iCnt = 1, dbsc(iCnttp)%nCntr
                     A(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)

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
                        If ( TstFnc(dc(mdc+iCnt)%iCoSet,
     &                     iIrrep,iComp,dc(mdc+iCnt)%nStab) ) Then
                           nDisp = nDisp + 1
                           If (Direct(nDisp)) Then
                              Temp(nDisp) = Temp(nDisp) - Tempd(iCar+1)
                           End If
                        End If
                     End Do

*
                  End Do
 103              Continue
                  mdc = mdc + dbsc(iCnttp)%nCntr
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
            ZB = dbsc(jCnttp)%Charge
            If (ZB.eq.Zero) Go To 212
            If (dbsc(jCnttp)%pChrg) Go To 212
            If (dbsc(jCnttp)%Frag) Go To 212
            ZAZB = ZA * ZB
            Do jCnt = 1, dbsc(jCnttp)%nCntr
               B(1:3)=dbsc(jCnttp)%Coor(1:3,jCnt)
*
*              Find the DCR for the two centers
*
               Call DCR(LmbdR,
     &                  iStb,nStb,
     &                  dc(ndc+jCnt)%iStab,dc(ndc+jCnt)%nStab,
     &                  iDCRR,nDCRR)
*
               PreFct = ZAZB*DBLE(nIrrep)/DBLE(LmbdR)
               Do iR = 0, nDCRR-1
                  Call OA(iDCRR(iR),B,RB)
                  nOp = NrOpr(iDCRR(iR))
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
                  If (dbsc(jCnttp)%ECP) Then
*--------------------Add contribution from M1 operator
                     Cnt0M1=Zero
                     Cnt1M1=Zero
                     Do iM1xp=1, dbsc(jCnttp)%nM1
                       Gamma =dbsc(jCnttp)%M1xp(iM1xp)
                       CffM1 =dbsc(jCnttp)%M1cf(iM1xp)
                       Cnt0M1= Cnt0M1+(CffM1*Exp(-Gamma*r12**2))
                       Cnt1M1= Cnt1M1+Gamma*(CffM1*Exp(-Gamma*r12**2))
                     End Do
                     fab=fab+Cnt0M1
                     dfab=dfab-Two*r12*Cnt1M1
*--------------------Add contribution from M2 operator
                     Cnt0M2=Zero
                     Cnt1M2=Zero
                     Do iM2xp=1, dbsc(jCnttp)%nM2
                       Gamma =dbsc(jCnttp)%M2xp(iM2xp)
                       CffM2 =dbsc(jCnttp)%M2cf(iM2xp)
                       Cnt0M2= Cnt0M2+(CffM2*Exp(-Gamma*r12**2))
                       Cnt1M2= Cnt1M2+Gamma*(CffM2*Exp(-Gamma*r12**2))
                     End Do
                     fab=fab+r12*Cnt0M2
                     dfab=dfab+(Cnt0M2-Two*r12**2*Cnt1M2)
                  End If
                  df_dr=(dfab*r12-fab)/r12**2
*
                  nDisp = IndDsp(ndc+jCnt,iIrrep)
                  igv=nIrrep/dc(ndc+jCnt)%nStab
                  Do iCar = 0, 2
                     dr_dB=-(A(iCar+1)-RB(iCar+1))/r12
                     iComp = 2**iCar
                     If ( TstFnc(dc(ndc+jCnt)%iCoSet,
     &                  iIrrep,iComp,dc(ndc+jCnt)%nStab) ) Then
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
            ndc = ndc + dbsc(jCnttp)%nCntr
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
      Return
      End
