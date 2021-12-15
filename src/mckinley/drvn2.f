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
* Copyright (C) 1991, Roland Lindh                                     *
*               1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine DrvN2(Hess,nGrad)
************************************************************************
*                                                                      *
* Object: to compute the molecular gradient contribution due to the    *
*         nuclear repulsion energy.                                    *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             October 1991                                             *
*             Anders Bernhardsson, Dept. of Theoretical Chemistry,     *
*             University of Lund, SWEDEN                               *
*             September 1995                                           *
************************************************************************
      use Basis_Info
      use Center_Info
      use PCM_arrays
      use Symmetry_Info, only: nIrrep, iChTbl
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "disp.fh"
#include "disp2.fh"
#include "WrkSpc.fh"
#include "rctfld.fh"
      Real*8 A(3), B(3), RB(3), Hess(nGrad*(nGrad+1)/2),prmt(0:7),
     &       C(3), D(3), SD(3)
      Integer iDCRR(0:7),IndGrd(3,2,0:7),ii(2), iStb(0:7),jStb(0:7),
     &        iDCRS(0:7),IndHss(2,3,2,3,0:7),nop(2),kop(2)
      Logical EQ, TstFnc,TF, NoLoop
      Data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
*                                                                      *
************************************************************************
*                                                                      *
*     Statement Function
*
      xPrmt(i,j) = Prmt(iAnd(i,j))
      TF(mdc,iIrrep,iComp) = TstFnc(dc(mdc)%iCoSet,
     &                              iIrrep,iComp,dc(mdc)%nStab)
      iTri(i1,i2)=Max(i1,i2)*(Max(i1,i2)-1)/2+Min(i1,i2)
*                                                                      *
************************************************************************
*                                                                      *
*
c     iRout = 33
c     iPrint = nPrint(iRout)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the nuclear repulsion contributions                      *
*                                                                      *
************************************************************************
*                                                                      *
      nHess = nGrad*(nGrad+1)/2
      call dcopy_(nHess,[Zero],0,Hess,1)
*
      mdc = 0
*-----Loop over centers with the same change
      Do iCnttp = 1, nCnttp
       ZA = dbsc(iCnttp)%Charge
       If (ZA.eq.Zero) Go To 101
*--------Loop over all unique centers of this group
       Do 110 iCnt = 1, dbsc(iCnttp)%nCntr
        A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)
*
         ndc = 0
         Do jCnttp = 1, iCnttp
           ZB=dbsc(jCnttp)%Charge
           If (ZB.eq.Zero) Go To 201
           ZAZB = ZA * ZB
           jCntMx = dbsc(jCnttp)%nCntr
           If (iCnttp.eq.jCnttp) jCntMx = iCnt
           Do jCnt = 1, jCntMx
             B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)
*
             Fact = One
*                 Factor due to resticted summation
             If (EQ(A,B)) Fact = Half
*
*                 Find the DCR for the two centers
*
             Call DCR(LmbdR,dc(mdc+iCnt)%iStab,dc(mdc+iCnt)%nStab,
     &                      dc(ndc+jCnt)%iStab,dc(ndc+jCnt)%nStab,
     &                      iDCRR,nDCRR)
*
             PreFct = Fact*ZAZB*DBLE(nIrrep)/DBLE(LmbdR)
             Do iR = 0, nDCRR-1
               Call OA(iDCRR(iR),B,RB)
               nOp(1) = NrOpr(0)
               nOp(2) = NrOpr(iDCRR(iR))
               kop(1)=0
               kop(2)=iDCRR(iR)
               If (EQ(A,RB)) Go To 301
               r12 = Sqrt((A(1)-RB(1))**2 +
     &                    (A(2)-RB(2))**2 +
     &                    (A(3)-RB(3))**2 )
*
*              The factor u/g will ensure that the value of the
*              gradient in symmetry adapted and no symmetry basis
*              will have the same value.
*
               fab = One
               dfab=Zero
               ddfab=Zero
*
               If (dbsc(iCnttp)%ECP) Then
*-----------------Add contibution from M1 operator
                  Cnt0M1=Zero
                  Cnt1M1=Zero
                  Cnt2M1=Zero
                  Do iM1xp = 1, dbsc(iCnttp)%nM1
                     Gamma =dbsc(iCnttp)%M1xp(iM1xp)
                     CffM1 =dbsc(iCnttp)%M1cf(iM1xp)
                     Cnt0M1=Cnt0M1+(CffM1*Exp(-Gamma*r12**2))
                     Cnt1M1=Cnt1M1+Gamma*(CffM1*Exp(-Gamma*r12**2))
                     Cnt2M1=Cnt2M1+Gamma**2*(CffM1*Exp(-Gamma*r12**2))
                  End Do
                  fab=fab+Cnt0M1
                  dfab=dfab-Two*r12*Cnt1M1
                  ddfab=-Two*Cnt1M1+Four*r12**2*Cnt2M1
*-----------------Add contibution from M2 operator
                  Cnt0M2=Zero
                  Cnt1M2=Zero
                  Cnt2M2=Zero
                  Do iM2xp = 1, dbsc(iCnttp)%nM2
                     Gamma =dbsc(iCnttp)%M2xp(iM2xp)
                     CffM2 =dbsc(iCnttp)%M2cf(iM2xp)
                     Cnt0M2=Cnt0M2+(CffM2*Exp(-Gamma*r12**2))
                     Cnt1M2=Cnt1M2+Gamma*(CffM2*Exp(-Gamma*r12**2))
                     Cnt2M2=Cnt2M2+Gamma**2*(CffM2*Exp(-Gamma*r12**2))
                  End Do
                  fab=fab+r12*Cnt0M2
                  dfab=dfab+Cnt0M2-Two*r12**2*Cnt1M2
                  ddfab=ddfab-Six**r12*Cnt1M2+Four*r12*Three*Cnt2M2
               End If
               If (dbsc(jCnttp)%ECP) Then
*-----------------Add contibution from M1 operator
                  Cnt0M1=Zero
                  Cnt1M1=Zero
                  Cnt2M1=Zero
                  Do iM1xp = 1, dbsc(jCnttp)%nM1
                     Gamma =dbsc(jCnttp)%M1xp(iM1xp)
                     CffM1 =dbsc(jCnttp)%M1cf(iM1xp)
                     Cnt0M1=Cnt0M1+(CffM1*Exp(-Gamma*r12**2))
                     Cnt1M1=Cnt1M1+Gamma*(CffM1*Exp(-Gamma*r12**2))
                     Cnt2M1=Cnt2M1+Gamma**2*(CffM1*Exp(-Gamma*r12**2))
                  End Do
                  fab=fab+Cnt0M1
                  dfab=dfab-Two*r12*Cnt1M1
                  ddfab=-Two*Cnt1M1+Four*r12**2*Cnt2M1
*-----------------Add contibution from M2 operator
                  Cnt0M2=Zero
                  Cnt1M2=Zero
                  Cnt2M2=Zero
                  Do iM2xp = 1, dbsc(jCnttp)%nM2
                     Gamma =dbsc(jCnttp)%M2xp(iM2xp)
                     CffM2 =dbsc(jCnttp)%M2cf(iM2xp)
                     Cnt0M2=Cnt0M2+(CffM2*Exp(-Gamma*r12**2))
                     Cnt1M2=Cnt1M2+Gamma*(CffM2*Exp(-Gamma*r12**2))
                     Cnt2M2=Cnt2M2+Gamma**2*(CffM2*Exp(-Gamma*r12**2))
                  End Do
                  fab=fab+r12*Cnt0M2
                  dfab=dfab+Cnt0M2-Two*r12**2*Cnt1M2
                  ddfab=ddfab-Six**r12*Cnt1M2+Four*r12*Three*Cnt2M2
               End If
*
               df_dr=(dfab*r12-fab)/r12**2
               d2f_dr2= ( (ddfab*r12)   * r12**2
     &                  - (dfab*r12-fab)* Two*r12 ) / r12**4
*
               Call ICopy(nirrep*36,[0],0,Indhss,1)
               Call ICopy(nirrep*6,[0],0,indgrd,1)
*
*          Determine which displacement in all IR's, each center is *
*          associated with
*
               nnIrrep=nIrrep
               If (sIrrep) nnIrrep=1

               Do  iIrrep=0,nnIrrep-1
                 nDisp1 = IndDsp(mdc+iCnt,iIrrep)
                 nDisp2 = IndDsp(ndc+jCnt,iIrrep)
                 Do  iCar = 0,2
                   iComp = 2**iCar
                   If ( TF(mdc+iCnt,iIrrep,iComp)) Then
                     nDisp1 = nDisp1 + 1
                     IndGrd(iCar+1,1,iIrrep) = nDisp1
                   Else
                     IndGrd(iCar+1,1,iIrrep)=0
                   End If
                   iComp = 2**iCar
                   If ( TF(ndc+jCnt,iIrrep,iComp)) Then
                     nDisp2 = nDisp2 + 1
                     IndGrd(iCar+1,2,iIrrep) = nDisp2
                   Else
                     IndGrd(iCar+1,2,iIrrep)=0
                   End If
                 End Do ! iCar
               End Do   ! iIrrep
*
*          Determine index for each 2'nd derivative
*
*
           Do  iIrrep=0,nnIrrep-1
             Do iAtom=1,2
               Do iCar=1,3
                 Do jAtom=1,iAtom
                   jCar_Max=3
                   if (iAtom.eq.jAtom) jCar_Max=iCar
                   Do jCar=1,jCar_Max
                     If ((IndGrd(iCar,iAtom,iIrrep).gt.0) .and.
     &                   (IndGrd(jCar,jAtom,iIrrep).gt.0)) Then
*
                       IndHss(iAtom,iCar,jAtom,jCar,iIrrep)=
     &                  iTri(IndGrd(iCar,iAtom,iIrrep),
     &                       IndGrd(jCar,jAtom,iIrrep))
*
                     Else
*
                       IndHss(iAtom,iCar,jAtom,jCar,iIrrep)=0
*
                     End If
                   End Do ! jCar
                 End Do   ! jAtom
               End Do     ! iCar
             End Do       ! iAtom
           End Do         ! iIrrep
*
           ii(1)=dc(mdc+icnt)%nStab
           ii(2)=dc(ndc+jcnt)%nStab
*
           Do iIrrep=0,nnIrrep-1
             Do iCent=1,2
               Do jCent=1,iCent
                 Do  iCar = 1, 3
                   jCar_Max=3
                   If (iCent.eq.jCent) jCar_Max=iCar
                   Do  jCar=1,jCar_Max
                    iCh1=2**(iCar-1)
                    iCh2=2**(jCar-1)
                    g=DBLE(iChTbl(iIrrep,nOp(icent)))*
     &                    xPrmt(kOp(icent),iCh1)*
     &                    DBLE(ii(icent))/
     &                    DBLE(nIrrep)
                    g=g*DBLE(iChTbl(iIrrep,nOp(jcent)))*
     &                     xPrmt(kOp(jcent),iCh2)*
     &                     DBLE(ii(jcent))/
     &                     DBLE(nIrrep)
                    g=g*(-One)**(icent+jcent)
                    if ((iCent.ne.jCent).and.(iCar.eq.jCar).and.
     &               (Abs(indgrd(iCar,iCent,iIrrep)).eq.
     &                Abs(indgrd(jCar,jCent,iIrrep)))) Then
                      ps=Two
                    Else
                      ps=One
                    End if

                    Index=indHss(iCent,iCar,jCent,jCar,iIrrep)
                    If (index.ne.0) Then
                       dr_dAi=(A(iCar)-RB(iCar))/r12
                       dr_dAj=(A(jCar)-RB(jCar))/r12
                       d2r_dAidAj=-(A(iCar)-RB(iCar))*dr_dAj
                       If (iCar.eq.jCar) d2r_dAidAj=d2r_dAidAj+r12
                       d2r_dAidAj=d2r_dAidAj/r12**2
                       Hess(Index) = Hess(index)+ g*PreFct*ps
     &                             *(d2r_dAidAj*df_dr +
     &                               dr_dAi*dr_dAj*d2f_dr2)
                    End If
                   End Do ! jCar
                 End Do   ! iCar
               End Do     ! jCent
             End Do       ! iCent
           End Do         ! iIrrep
*
*          call triprt(' ',' ',Hess,ldisp(0))
 301       Continue
           End Do
*
               End Do
 201           Continue
               ndc = ndc + dbsc(jCnttp)%nCntr
            End Do
 110     Continue
 101     Continue
         mdc = mdc + dbsc(iCnttp)%nCntr
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     PCM contributions
*                                                                      *
************************************************************************
*                                                                      *
      If (PCM) Then
*
*     We will have three contributions here.
*                                                                      *
************************************************************************
*                                                                      *
*
*     1) Process the contribution
*
*        Sum(i) q_i V_in^xy
*
*---- Loop over tiles
*
*
      Do iTs = 1, nTs
         ZA   = PCM_SQ(1,iTs)+PCM_SQ(2,iTs)
         NoLoop = ZA.eq.Zero
         If (NoLoop) Go To 112
         ZA = ZA / DBLE(nIrrep)
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
               B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)
*
*              Find the DCR for the two centers
*
               Call DCR(LmbdR,iStb,nStb,
     &                        dc(ndc+jCnt)%iStab,dc(ndc+jCnt)%nStab,
     &                        iDCRR,nDCRR)
*
               PreFct = ZAZB * DBLE(nIrrep)/DBLE(LmbdR)
               Do iR = 0, nDCRR-1
                  Call OA(iDCRR(iR),B,RB)
                  nOp(1) = NrOpr(0)
                  nOp(2) = NrOpr(iDCRR(iR))
                  r12 = Sqrt((A(1)-RB(1))**2 +
     &                       (A(2)-RB(2))**2 +
     &                       (A(3)-RB(3))**2 )
*
*                 The factor u/g will ensure that the value of the
*                 gradient in symmetry adapted and no symmetry basis
*                 will have the same value.
*
                  fab  =One
                  dfab =Zero
                  ddfab=Zero
                  If (dbsc(jCnttp)%ECP) Then
*--------------------Add contibution from M1 operator
                     Cnt0M1=Zero
                     Cnt1M1=Zero
                     Cnt2M1=Zero
                     Do iM1xp = 1, dbsc(jCnttp)%nM1
                        Gamma =dbsc(jCnttp)%M1xp(iM1xp)
                        CffM1 =dbsc(jCnttp)%M1cf(iM1xp)
                        Cnt0M1=Cnt0M1+(CffM1*Exp(-Gamma*r12**2))
                        Cnt1M1=Cnt1M1+Gamma*(CffM1*Exp(-Gamma*r12**2))
                        Cnt2M1=Cnt2M1
     &                        +Gamma**2*(CffM1*Exp(-Gamma*r12**2))
                     End Do
                     fab=fab+Cnt0M1
                     dfab=dfab-Two*r12*Cnt1M1
                     ddfab=-Two*Cnt1M1+Four*r12**2*Cnt2M1
*--------------------Add contibution from M2 operator
                     Cnt0M2=Zero
                     Cnt1M2=Zero
                     Cnt2M2=Zero
                     Do iM2xp = 1, dbsc(jCnttp)%nM2
                        Gamma =dbsc(jCnttp)%M2xp(iM2xp)
                        CffM2 =dbsc(jCnttp)%M2cf(iM2xp)
                        Cnt0M2=Cnt0M2+(CffM2*Exp(-Gamma*r12**2))
                        Cnt1M2=Cnt1M2+Gamma*(CffM2*Exp(-Gamma*r12**2))
                        Cnt2M2=Cnt2M2
     &                        +Gamma**2*(CffM2*Exp(-Gamma*r12**2))
                     End Do
                     fab=fab+r12*Cnt0M2
                     dfab=dfab+Cnt0M2-Two*r12**2*Cnt1M2
                     ddfab=ddfab-Six**r12*Cnt1M2+Four*r12*Three*Cnt2M2
                  End If
*
                  df_dr=(dfab*r12-fab)/r12**2
                  d2f_dr2= ( (ddfab*r12)   * r12**2
     &                     - (dfab*r12-fab)* Two*r12) / r12**4
*
               Call ICopy(nirrep*36,[0],0,Indhss,1)
               Call ICopy(nirrep*6,[0],0,indgrd,1)
*
*          Determine which displacement in all IRs, each center is *
*          associated with
*
               nnIrrep=nIrrep
               If (sIrrep) nnIrrep=1

               Do  iIrrep=0,nnIrrep-1
                 nDisp1 = IndDsp(ndc+jCnt,iIrrep)
                 Do  iCar = 0,2
                   iComp = 2**iCar
                   If ( TF(ndc+jCnt,iIrrep,iComp)) Then
                     nDisp1 = nDisp1 + 1
                     IndGrd(iCar+1,1,iIrrep) = nDisp1
                   Else
                     IndGrd(iCar+1,1,iIrrep)=0
                   End If
                 End Do ! iCar
               End Do   ! iIrrep
*
*          Determine index for each 2nd derivative
*
*          Note that each term is only associated with one basis
*          set center.
*
*
           iAtom=1
           jAtom=1
           Do  iIrrep=0,nnIrrep-1
             Do iCar=1,3
               jCar_Max=iCar
               Do jCar=1,jCar_Max
                 If ((IndGrd(iCar,iAtom,iIrrep).gt.0) .and.
     &               (IndGrd(jCar,jAtom,iIrrep).gt.0)) Then
*
                   IndHss(iAtom,iCar,jAtom,jCar,iIrrep)=
     &              iTri(IndGrd(iCar,iAtom,iIrrep),
     &                   IndGrd(jCar,jAtom,iIrrep))
*
                 Else
*
                   IndHss(iAtom,iCar,jAtom,jCar,iIrrep)=0
*
                 End If
               End Do ! jCar
             End Do     ! iCar
           End Do         ! iIrrep
*
           ii(1)=dc(ndc+jcnt)%nStab
*
           iCent=1
           jCent=1
           Do iIrrep=0,nnIrrep-1
              Do  iCar = 1, 3
                jCar_Max=iCar
                Do  jCar=1,jCar_Max
                 iCh1=2**(iCar-1)
                 iCh2=2**(jCar-1)
                 g=DBLE(iChTbl(iIrrep,nOp(icent)))*
     &                 xPrmt(kOp(icent),iCh1)*
     &                 DBLE(ii(icent))/
     &                 DBLE(nIrrep)
                 g=g*DBLE(iChTbl(iIrrep,nOp(jcent)))*
     &                  xPrmt(kOp(jcent),iCh2)*
     &                  DBLE(ii(jcent))/
     &                  DBLE(nIrrep)
                 g=g*(-One)**(icent+jcent)
*
                 Index=indHss(iCent,iCar,jCent,jCar,iIrrep)
                 If (Index.ne.0) Then
                    dr_dAi=(A(iCar)-RB(iCar))/r12
                    dr_dAj=(A(jCar)-RB(jCar))/r12
                    d2r_dAidAj=-(A(iCar)-RB(iCar))*dr_dAj
                    If (iCar.eq.jCar) d2r_dAidAj=d2r_dAidAj+r12
                    d2r_dAidAj=d2r_dAidAj/r12**2
                    Hess(Index) = Hess(Index)+ g*PreFct
     &                          *(d2r_dAidAj*df_dr +
     &                            dr_dAi*dr_dAj*d2f_dr2)
                 End If
                End Do ! jCar
              End Do   ! iCar
           End Do         ! iIrrep
*
*          Call TriPrt(' ',' ',Hess,ldisp(0))
*
               End Do      ! End loop over DCR operators, iR
*
            End Do         ! End over centers, jCnt
 212        Continue
            ndc = ndc + dbsc(jCnttp)%nCntr
         End Do            ! End over basis set types, jCnttp
 112     Continue
      End Do               ! End of tiles
*                                                                      *
************************************************************************
*                                                                      *
*     2) Process the contribution
*
*        Sum(i,j) V_i,n^x Q_ij V_j,n^y
*
      Do iTs = 1, nTs
         A(1:3) = PCMTess(1:3,iTs)
*
*------- Tile only stabilized by the unit operator
*
         nStb=1
         iStb(0)=0
*
      Do jTs = 1, iTs
         Fact=Two
         If (jTs.eq.iTs) Fact=One
         Q_ij = DMElm(nTs,iTs,jTs,PCMDM)
         NoLoop = Q_ij.eq.Zero
         If (NoLoop) Go To 122
         C(1:3) = PCMTess(1:3,jTs)
*
         mStb=1
         jStb(0)=0
*
*        Loop over the basis functions
*
         mdc = 0
         Do iCnttp = 1, nCnttp
            ZA = dbsc(iCnttp)%Charge
            If (ZA.eq.Zero) Go To 222
            If (dbsc(iCnttp)%pChrg) Go To 222
            If (dbsc(iCnttp)%Frag) Go To 222
*
            Do iCnt = 1, dbsc(iCnttp)%nCntr
               B(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)
*
*              Find the DCR for the two centers (
*
               Call DCR(LmbdR,iStb,nStb,
     &                        dc(mdc+iCnt)%iStab,dc(mdc+iCnt)%nStab,
     &                        iDCRR,nDCRR)
*
               PreFct_AB = DBLE(nIrrep)/DBLE(LmbdR)
               Do iR = 0, nDCRR-1
                  Call OA(iDCRR(iR),B,RB)
                  nOp(1) = NrOpr(iDCRR(iR))
                  r12_AB = Sqrt((A(1)-RB(1))**2 +
     &                          (A(2)-RB(2))**2 +
     &                          (A(3)-RB(3))**2 )
                  fab  =One
                  dfab =Zero
                  If (dbsc(iCnttp)%ECP) Then
*--------------------Add contibution from M1 operator
                     Cnt0M1=Zero
                     Cnt1M1=Zero
                     Do iM1xp = 1, dbsc(iCnttp)%nM1
                        Gamma =dbsc(iCnttp)%M1xp(iM1xp)
                        CffM1 =dbsc(iCnttp)%M1cf(iM1xp)
                        Cnt0M1=Cnt0M1+(CffM1*Exp(-Gamma*r12_AB**2))
                        Cnt1M1=Cnt1M1
     &                        +Gamma*(CffM1*Exp(-Gamma*r12_AB**2))
                     End Do
                     fab=fab+Cnt0M1
                     dfab=dfab-Two*r12_AB*Cnt1M1
*--------------------Add contibution from M2 operator
                     Cnt0M2=Zero
                     Cnt1M2=Zero
                     Do iM2xp = 1, dbsc(iCnttp)%nM2
                        Gamma =dbsc(iCnttp)%M2xp(iM2xp)
                        CffM2 =dbsc(iCnttp)%M2cf(iM2xp)
                        Cnt0M2=Cnt0M2+(CffM2*Exp(-Gamma*r12_AB**2))
                        Cnt1M2=Cnt1M2
     &                        +Gamma*(CffM2*Exp(-Gamma*r12_AB**2))
                     End Do
                     fab=fab+r12_AB*Cnt0M2
                     dfab=dfab+Cnt0M2-Two*r12_AB**2*Cnt1M2
                  End If
                  df_dr_AB=(dfab*r12_AB-fab)/r12_AB**2
*
         ndc = 0
         Do jCnttp = 1, nCnttp
            ZB = dbsc(jCnttp)%Charge
            If (ZB.eq.Zero) Go To 232
            If (dbsc(jCnttp)%pChrg) Go To 232
            If (dbsc(jCnttp)%Frag) Go To 232
*
            Do jCnt = 1, dbsc(jCnttp)%nCntr
               D(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)
*
*              Find the DCR for the two centers (
*
               Call DCR(LmbdS,iStb,nStb,
     &                        dc(ndc+jCnt)%iStab,dc(ndc+jCnt)%nStab,
     &                        iDCRS,nDCRS)
*
               PreFct_CD = DBLE(nIrrep)/DBLE(LmbdS)
               Do iS = 0, nDCRS-1
                  Call OA(iDCRS(iS),D,SD)
                  nOp(2) = NrOpr(iDCRS(iS))
                  r12_CD = Sqrt((C(1)-SD(1))**2 +
     &                          (C(2)-SD(2))**2 +
     &                          (C(3)-SD(3))**2 )
*
                  fcd  =One
                  dfcd =Zero
                  If (dbsc(jCnttp)%ECP) Then
*--------------------Add contibution from M1 operator
                     Cnt0M1=Zero
                     Cnt1M1=Zero
                     Do iM1xp = 1, dbsc(jCnttp)%nM1
                        Gamma =dbsc(jCnttp)%M1xp(iM1xp)
                        CffM1 =dbsc(jCnttp)%M1cf(iM1xp)
                        Cnt0M1=Cnt0M1+(CffM1*Exp(-Gamma*r12_CD**2))
                        Cnt1M1=Cnt1M1
     &                        +Gamma*(CffM1*Exp(-Gamma*r12_CD**2))
                     End Do
                     fcd=fcd+Cnt0M1
                     dfcd=dfcd-Two*r12_CD*Cnt1M1
*--------------------Add contibution from M2 operator
                     Cnt0M2=Zero
                     Cnt1M2=Zero
                     Do iM2xp = 1, dbsc(jCnttp)%nM2
                        Gamma =dbsc(jCnttp)%M2xp(iM2xp)
                        CffM2 =dbsc(jCnttp)%M2cf(iM2xp)
                        Cnt0M2=Cnt0M2+(CffM2*Exp(-Gamma*r12_CD**2))
                        Cnt1M2=Cnt1M2
     &                        +Gamma*(CffM2*Exp(-Gamma*r12_CD**2))
                     End Do
                     fcd=fcd+r12_CD*Cnt0M2
                     dfcd=dfcd+Cnt0M2-Two*r12_CD**2*Cnt1M2
                  End If
                  df_dr_CD=(dfcd*r12_CD-fcd)/r12_CD**2
*
               Call ICopy(nirrep*36,[0],0,Indhss,1)
               Call ICopy(nirrep*6,[0],0,indgrd,1)
*
*          Determine which displacement in all IR's, each center is *
*          associated with
*
               nnIrrep=nIrrep
               If (sIrrep) nnIrrep=1

               Do  iIrrep=0,nnIrrep-1
                 nDisp1 = IndDsp(mdc+iCnt,iIrrep)
                 nDisp2 = IndDsp(ndc+jCnt,iIrrep)
                 Do  iCar = 0,2
                   iComp = 2**iCar
*
                   If ( TF(mdc+iCnt,iIrrep,iComp)) Then
                     nDisp1 = nDisp1 + 1
                     IndGrd(iCar+1,1,iIrrep) = nDisp1
                   Else
                     IndGrd(iCar+1,1,iIrrep)=0
                   End If
*
                   If ( TF(ndc+jCnt,iIrrep,iComp)) Then
                     nDisp2 = nDisp2 + 1
                     IndGrd(iCar+1,2,iIrrep) = nDisp2
                   Else
                     IndGrd(iCar+1,2,iIrrep)=0
                   End If
                 End Do ! iCar
               End Do   ! iIrrep
*
*          Determine index for each 2'nd derivative
*
*
           Do  iIrrep=0,nnIrrep-1
             Do iAtom=1,2
               Do iCar=1,3
                 Do jAtom=1,iAtom
                   jCar_Max=3
                   if (iAtom.eq.jAtom) jCar_Max=iCar
                   Do jCar=1,jCar_Max
                     If ((IndGrd(iCar,iAtom,iIrrep).gt.0) .and.
     &                   (IndGrd(jCar,jAtom,iIrrep).gt.0)) Then
*
                       IndHss(iAtom,iCar,jAtom,jCar,iIrrep)=
     &                  iTri(IndGrd(iCar,iAtom,iIrrep),
     &                       IndGrd(jCar,jAtom,iIrrep))
*
                     Else
*
                       IndHss(iAtom,iCar,jAtom,jCar,iIrrep)=0
*
                     End If
                   End Do ! jCar
                 End Do   ! jAtom
               End Do     ! iCar
             End Do       ! iAtom
           End Do         ! iIrrep
*
           ii(1)=dc(mdc+icnt)%nStab
           ii(2)=dc(ndc+jcnt)%nStab
*
*          Note that we have two different cases here, depending on if
*          iTs=jTs or not!
*          For iTs=jTs and iCent.eq.jCent we do
*          dV_i/dx*dV_i/dy only and exclude dV_i/dy*dV_i/dx since they
*          are the same. For iTs=/=jTs we need both dV_i/dx*dV_j/dy and
*          dV_i/dy*dV_j/dx.
*
           Do iIrrep=0,nnIrrep-1
             Do iCent=1,2
               Do jCent=1,iCent
                 Do  iCar = 1, 3
                   jCar_Max=3
                   If (iCent.eq.jCent.and.iTs.eq.jTs) jCar_Max=iCar
                   Do  jCar=1,jCar_Max
                    iCh1=2**(iCar-1)
                    iCh2=2**(jCar-1)
                    g=DBLE(iChTbl(iIrrep,nOp(icent)))*
     &                    xPrmt(kOp(icent),iCh1)*
     &                    DBLE(ii(icent))/
     &                    DBLE(nIrrep)
                    g=g*DBLE(iChTbl(iIrrep,nOp(jcent)))*
     &                     xPrmt(kOp(jcent),iCh2)*
     &                     DBLE(ii(jcent))/
     &                     DBLE(nIrrep)
                    g=g*(-One)**(icent+jcent)
*
                    If ((iCent.ne.jCent).and.(iCar.eq.jCar).and.
     &               (Abs(IndGrd(iCar,iCent,iIrrep)).eq.
     &                Abs(IndGrd(jCar,jCent,iIrrep)))) Then
                      ps=Two
                    Else
                      ps=One
                    End if

                    Index=IndHss(iCent,iCar,jCent,jCar,iIrrep)
                    If (Index.ne.0) Then
                       dr_dB=-(A(iCar)-RB(iCar))/r12_AB
                       dr_dD=-(C(jCar)-SD(jCar))/r12_CD
                       Hess(Index) = Hess(Index)
     &                             + Fact*g*ps
     &                             * ZA * ZA * Q_ij
     &                             * PreFct_AB * dr_dB * df_dr_AB
     &                             * PreFct_CD * dr_dD * df_dr_CD
                    End If
                   End Do ! jCar
                 End Do   ! iCar
               End Do     ! jCent
             End Do       ! iCent
           End Do         ! iIrrep
*
*          Call TriPrt(' ',' ',Hess,ldisp(0))
*
               End Do      ! iS
*
            End Do         ! jCnt
 232        Continue
            ndc = ndc + dbsc(jCnttp)%nCntr
         End Do            ! jCnttp
*
               End Do      ! iR
*
            End Do         ! iCnt
 222        Continue
            mdc = mdc + dbsc(iCnttp)%nCntr
         End Do            ! jCnttp
 122     Continue
      End Do               ! jTs
      End Do               ! iTs
*                                                                      *
************************************************************************
*                                                                      *
*       Add additional contributions
*
        nPCMHss = nGrad * nGrad
        Call Get_nAtoms_All(nAtoms)
        Call GetMem('PCM_Hss','Allo','Real',ip_pcmhss,nPCMHss)
        Call GetMem('Der1','Allo','Real',ip_Der1,nTs)
        Call GetMem('DerDM','Allo','Real',ip_DerDM,nTs*nTs)
        Call GetMem('Temp','Allo','Real',ip_Temp,nTs*nTs)
        Call Cav_Hss(nAtoms,nGrad,nTs,nS,Eps,PCMSph,
     &               PCMiSph,PCM_N,PCMTess,PCM_SQ,
     &               PCMDM,Work(ip_Der1),Work(ip_DerDM),Work(ip_Temp),
     &               dTes,DPnt,dRad,dCntr,Work(ip_pcmhss),nPCMHss)
        Call GetMem('PCM_Hss','Free','Real',ip_pcmhss,nPCMHss)
        Call GetMem('Der1','Free','Real',ip_Der1,nTs)
        Call GetMem('DerDM','Free','Real',ip_DerDM,nTs*nTs)
        Call GetMem('Temp','Free','Real',ip_Temp,nTs*nTs)
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
