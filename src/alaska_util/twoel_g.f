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
* Copyright (C) 1990,1992, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine TwoEl_g(Coor,iAnga,iCmp,iShell,iShll,iAO,
     &                   iStb,jStb,kStb,lStb,nRys,
     &                   Data1,nab,nHmab,nData1,Data2,ncd,nHmcd,nData2,
     &                   Pren,Prem,
     &                   nAlpha,iPrInc, nBeta,jPrInc,
     &                   nGamma,kPrInc,nDelta,lPrInc,
     &                   Coeff1,iBasi,Coeff2,jBasj,
     &                   Coeff3,kBask,Coeff4,lBasl,
     &                   Zeta,ZInv,P,nZeta,Eta,EInv,Q,nEta,
     &                   xA,xB,xG,xD,Grad,nGrad,IfGrad,IndGrd,
     &                   PSO,nPSO,Wrk2,nWrk2,Aux,nAux,Shijij)
************************************************************************
*                                                                      *
* Object: to generate the SO integrals for four fixed centers and      *
*         fixed basis set types.                                       *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*          Roland Lindh, Dept. of Theoretical Chemistry, University of *
*          Lund, SWEDEN. Modified to gradients, January '92.           *
************************************************************************
      use Real_Spherical
      use Basis_Info
      use Center_Info
      use Phase_Info
      use Real_Info, only: ChiI2
      use Temporary_Parameters, only: IsChi
      use Symmetry_Info, only: nIrrep
      Use Iso_C_Binding
      Implicit Real*8 (A-H,O-Z)
      External TERI1, ModU2, vCff2D
#include "Molcas.fh"
#include "ndarray.fh"
#include "real.fh"
#include "print.fh"
#include "disp.fh"
      Real*8, Target ::
     &       Data1(nZeta*(nDArray+2*nab)+nDScalar+nHmab,nData1),
     &       Data2(nEta *(nDArray+2*ncd)+nDScalar+nHmcd,nData2)
      Integer, Pointer :: iData1(:), iData2(:)
      Real*8 Coor(3,4), CoorM(3,4), CoorAC(3,2),
     &       xA(nZeta),xB(nZeta), xG(nEta), xD(nEta), Grad(nGrad),
     &       Zeta(nZeta), ZInv(nZeta), P(nZeta,3),
     &       Eta (nEta),  EInv(nEta ), Q(nEta, 3),
     &       Coeff1(nAlpha,iBasi), Coeff2(nBeta,jBasj),
     &       Coeff3(nGamma,kBask), Coeff4(nDelta,lBasl),
     &       PSO(iBasi*jBasj*kBask*lBasl,nPSO), Wrk2(nWrk2),
     &       Aux(nAux)
      Integer iDCRR(0:7), iDCRS(0:7), iDCRT(0:7), iStabN(0:7),
     &        iStabM(0:7), IndGrd(3,4), iAO(4),
     &        iAnga(4), iCmp(4), iShell(4), iShll(4),
     &        nOp(4), kOp(4), JndGrd(3,4), iuvwx(4)
      Logical Shijij, AeqB, CeqD, AeqC, ABeqCD,
     &        EQ, lEmpty, IfGrad(3,4),
     &        JfGrad(3,4), PreScr
#ifdef _DEBUG_
      Character ChOper(0:7)*3
      Data ChOper/' E ',' x ',' y ',' xy',' z ',' xz',' yz','xyz'/
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function to compute canonical index
*
      nElem(i) = (i+1)*(i+2)/2
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 12
      iPrint = nPrint(iRout)
*
      la = iAnga(1)
      lb = iAnga(2)
      lc = iAnga(3)
      ld = iAnga(4)
      iCmpa = iCmp(1)
      jCmpb = iCmp(2)
      kCmpc = iCmp(3)
      lCmpd = iCmp(4)
      iShlla = iShll(1)
      jShllb = iShll(2)
      kShllc = iShll(3)
      lShlld = iShll(4)
      IncZet=nAlpha*jPrInc
      IncEta=nGamma*lPrInc
      LmbdT=0
      nijkl = iBasi*jBasj*kBask*lBasl
      iuvwx(1) = dc(iStb)%nStab
      iuvwx(2) = dc(jStb)%nStab
      iuvwx(3) = dc(kStb)%nStab
      iuvwx(4) = dc(lStb)%nStab
      mab = nElem(la)*nElem(lb)
      mcd = nElem(lc)*nElem(ld)
      iW4 = 1
      If (jPrInc.ne.nBeta .or. lPrInc.ne.nDelta) Then
         iW2 = 1 + mab*mcd*nijkl
      Else
         iW2 = 1
      End If
*
      If (l2DI) Then
         iffab = ip_abG(nZeta,nHmab)-1
         iffabG= ip_abG(nZeta,nHmab)+nab*nZeta-1
         iffcd = ip_abG( nEta,nHmcd)-1
         iffcdG= ip_abG( nEta,nHmcd) +ncd*nEta-1
      Else
*--------Dummy pointer to assure that we will not be out
*        off bounds.
         iffab = ip_abG(nZeta,nHmab)-nZeta-1
         iffabG= ip_abG(nZeta,nHmab)-nZeta-1
         iffcd = ip_abG( nEta,nHmcd)-nZeta-1
         iffcdG= ip_abG( nEta,nHmcd)-nZeta-1
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Find the Double Coset Representatives for center A and B
*
      If (nIrrep.eq.1) Then
         nDCRR=1
         iDCRR(0)=0
         LmbdR=1
      Else
         Call DCR(LmbdR,dc(iStb)%iStab,dc(iStb)%nStab,
     &                  dc(jStb)%iStab,dc(jStb)%nStab,iDCRR,nDCRR)
      End If
#ifdef _DEBUG_
      If (iPrint.ge.99) Write (6,'(20A)') ' {R}=(',
     &      (ChOper(iDCRR(i)),',',i=0,nDCRR-1),')'
#endif
      u = DBLE(dc(iStb)%nStab)
      v = DBLE(dc(jStb)%nStab)
*
*-----Find stabilizer for center A and B
*
      If (nIrrep.eq.1) Then
         lStabM=1
         iStabM(0)=0
      Else
         Call Inter(dc(iStb)%iStab,dc(iStb)%nStab,
     &              dc(jStb)%iStab,dc(jStb)%nStab,iStabM,lStabM)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
*-----Find the Double Coset Representatives for center C and D.
*     Take care of redundancy if {f(aA)f(bB)}={f(cC)f(dD)}. Hence
*     we will only use unique combinations of operators from the
*     double coset representatives {R} and {S}.
*
      If (nIrrep.eq.1) Then
         nDCRS=1
         iDCRS(0)=0
         LmbdS=1
      Else
         Call DCR(LmbdS,dc(kStb)%iStab,dc(kStb)%nStab,
     &                  dc(lStb)%iStab,dc(lStb)%nStab,
     &                               iDCRS,nDCRS)
      End If
#ifdef _DEBUG_
      If (iPrint.ge.99) Write (6,'(20A)') ' {S}=(',
     &      (ChOper(iDCRS(i)),',',i=0,nDCRS-1),')'
#endif
      w = DBLE(dc(kStb)%nStab)
      x = DBLE(dc(lStb)%nStab)
*
*-----Find stabilizer for center C and D
*
      If (nIrrep.eq.1) Then
         lStabN=1
         iStabN(0)=0
      Else
         Call Inter(dc(kStb)%iStab,dc(kStb)%nStab,
     &              dc(lStb)%iStab,dc(lStb)%nStab,iStabN,lStabN)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
*-----Find the Double Coset Representatives for the two charge
*     distributions.
*
      If (nIrrep.eq.1) Then
         nDCRT=1
         iDCRT(0)=0
         LmbdT=1
      Else
         Call DCR(LmbdT,iStabM,lStabM,iStabN,lStabN,iDCRT,nDCRT)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----------Factor due to summation over DCR
*
            If (MolWgh.eq.1) Then
               Fact = DBLE(nIrrep) / DBLE(LmbdT)
            Else If (MolWgh.eq.0) Then
               Fact = u*v*w*x / DBLE(nIrrep**3 * LmbdT)
            Else
               Fact = Sqrt(u*v*w*x)/DBLE(nIrrep*LmbdT)
            End If
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
      nOp(1)=NrOpr(0)
      call dcopy_(3,Coor(1,1),1,CoorM(1,1),1)
      Do 100 lDCRR = 0, nDCRR-1
         nOp(2)=NrOpr(iDCRR(lDCRR))
         Call OA(iDCRR(lDCRR),Coor(1:3,2),CoorM(1:3,2))
         AeqB = EQ(CoorM(1,1),CoorM(1,2))
*
         MxDCRS = nDCRS-1
         Do 200 lDCRS = 0, MxDCRS
            call dcopy_(3,Coor(1,3),1,CoorM(1,3),1)
            Call OA(iDCRS(lDCRS),Coor(1:3,4),CoorM(1:3,4))
            CeqD = EQ(Coor(1,3),CoorM(1,4))
*
            Do 300 lDCRT = nDCRT-1, 0, -1
#ifdef _DEBUG_
               If (iPrint.ge.99) Write (6,'(6A)')
     &         ' R=',ChOper(iDCRR(lDCRR)),
     &         ', S=',ChOper(iDCRS(lDCRS)),
     &         ', T=',ChOper(iDCRT(lDCRT))
#endif
*
               nOp(3) = NrOpr(iDCRT(lDCRT))
               iDCRTS=iEor(iDCRT(lDCRT),iDCRS(lDCRS))
               nOp(4) = NrOpr(iDCRTS)
*
               Call OA(iDCRTS,Coor(1:3,4),CoorM(1:3,4))
               Call OA(iDCRT(lDCRT),Coor(1:3,3),CoorM(1:3,3))
*
#ifdef _DEBUG_
               If (iPrint.ge.59)
     &            Call RecPrt(' CoorM in TwoEl',' ',CoorM,3,4)
#endif
               AeqC = EQ(CoorM(1,1),CoorM(1,3))
               ABeqCD = AeqB .and. CeqD .and. AeqC
*--------------No contribution to gradient from one-center integrals
               If (ABeqCD) Go To 301
*
*--------------Modify which center we will differetiate with
*              respect to using the translational invariance.
*
               mGrad = 0
               Do 3333 iCar = 1, 3
*-----------------Copy to temporary arrays
                  JfGrad(iCar,1)=IfGrad(iCar,1)
                  JfGrad(iCar,2)=IfGrad(iCar,2)
                  JfGrad(iCar,3)=IfGrad(iCar,3)
                  JfGrad(iCar,4)=IfGrad(iCar,4)
                  JndGrd(iCar,1)=IndGrd(iCar,1)
                  JndGrd(iCar,2)=IndGrd(iCar,2)
                  JndGrd(iCar,3)=IndGrd(iCar,3)
                  JndGrd(iCar,4)=IndGrd(iCar,4)
*-----------------In case of four differentiations use
*                 the translational invariance to remove
*                 one or several of them.
                  If (IfGrad(iCar,1) .and.
     &                IfGrad(iCar,2) .and.
     &                IfGrad(iCar,3) .and.
     &                IfGrad(iCar,4) ) Then
*
                     nIdent=0
                     kCent=0
                     lCent=0
                     Do iCent = 1, 3
                        Do jCent = iCent+1, 4
                           If (EQ(CoorM(1,iCent),
     &                            CoorM(1,jCent))) Then
                              nIdent=nIdent+1
                              kCent=iCent
                              lCent=jCent
                           End If
                        End Do
                     End Do
*
*--------------------Remove a center which is unique and if possible
*                    if it occurs several times.
*
                     If (nIdent.eq.0) Then
*---------------------- Four center case, remove a center
                        jCent=1
                        Do iCent = 1, 4
                           If (iAnga(iCent).ne.0) jCent=iCent
                        End Do
                        JndGrd(iCar,jCent) = -JndGrd(iCar,jCent)
                        JfGrad(iCar,jCent) = .False.
                     Else If (nIdent.eq.1) Then
*---------------------- Three center case
                        iiCent=kCent
                        jjCent=lCent
                        JndGrd(iCar,iiCent) = 0
                        JfGrad(iCar,iiCent) = .False.
                        JndGrd(iCar,jjCent) = -JndGrd(iCar,jjCent)
                        JfGrad(iCar,jjCent) = .False.
                     Else If (nIdent.eq.2) Then
*---------------------- Two center case
                        iiCent=kCent
                        jjCent=lCent
                        ikl=2**(kCent-1)+2**(lCent-1)
                        Do iC = 1, 4
                           If (iAnd(ikl,2**(iC-1)).eq.0) Then
                              iCent=iC
                           End If
                        End Do
                        ikl=ikl+2**(iCent-1)
                        Do iC = 1, 4
                           If (iAnd(ikl,2**(iC-1)).eq.0) Then
                              jCent=iC
                           End If
                        End Do
                        ijMax=Max(iAnga(iCent),iAnga(jCent))
                        ijMin=Min(iAnga(iCent),iAnga(jCent))
                        klMax=Max(iAnga(kCent),iAnga(lCent))
                        klMin=Min(iAnga(kCent),iAnga(lCent))
                        If (klMin.gt.0) Then
                           iiCent=kCent
                           jjCent=lCent
                        Else If (ijMin.gt.0) Then
                           iiCent=iCent
                           jjCent=jCent
                        Else If (klMax.gt.0) Then
                           iiCent=kCent
                           jjCent=lCent
                        Else If (ijMax.gt.0) Then
                           iiCent=iCent
                           jjCent=jCent
                        End If
                        JndGrd(iCar,iiCent) = 0
                        JfGrad(iCar,iiCent) = .False.
                        JndGrd(iCar,jjCent) = -JndGrd(iCar,jjCent)
                        JfGrad(iCar,jjCent) = .False.
                     Else If (nIdent.eq.3) Then
*---------------------- Two center case
                        If (lCent.ne.4) Then
                           mCent=1
                        Else If (kCent.ne.3) Then
                           mCent=1
                        Else If (EQ(CoorM(1,1),CoorM(1,4))) Then
                           mCent=1
                        Else
                           mCent=2
                        End If
                        JndGrd(iCar,mCent) = 0
                        JfGrad(iCar,mCent) = .False.
                        JndGrd(iCar,kCent) = 0
                        JfGrad(iCar,kCent) = .False.
                        JndGrd(iCar,lCent) = -JndGrd(iCar,lCent)
                        JfGrad(iCar,lCent) = .False.
                     Else
                        Call WarningMessage(2,'Error in Twoel_g')
                        Write (6,*) ' Twoel: nIdent too large!'
                        Call Abend()
                     End If
                  End If
                  Do 4444 ixSh = 1, 4
                     If (JfGrad(iCar,ixSh)) mGrad=mGrad+1
 4444             Continue
 3333          Continue
               If (mGrad.eq.0) Go To 301
*
*--------------Find the proper centers to start of with the angular
*              momentum on. If la.eq.lb there will excist an
*              ambiguity to which center that angular momentum should
*              be accumulated on. In that case we will use A and C of
*              the order as defined by the basis functions types.
*
               If (iAnga(1).ge.iAnga(2)) Then
                  call dcopy_(3,CoorM(1,1),1,CoorAC(1,1),1)
               Else
                  call dcopy_(3,CoorM(1,2),1,CoorAC(1,1),1)
               End If
               If (iAnga(3).ge.iAnga(4)) Then
                   call dcopy_(3,CoorM(1,3),1,CoorAC(1,2),1)
               Else
                   call dcopy_(3,CoorM(1,4),1,CoorAC(1,2),1)
               End If
               kOp(1) = nOp(1)
               kOp(2) = nOp(2)
               kOp(3) = nOp(3)
               kOp(4) = nOp(4)

*
*--------------Desymmetrize the second order density matrix
*
*--------------(faA fbR(B) | fcT(C) fdTS(D))ijkl
*
               Call DesymP(iAnga,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &                     Shijij,iShll,iShell,iAO,kOp,nijkl,
     &                     Aux,nAux,Wrk2(iW2),PSO,nPSO)
*
               If (Fact.ne.One) Call DScal_(nijkl*
     &             iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4),
     &             Fact,Wrk2(iW2),1)
*
*--------------Backtransform 2nd order density matrix from spherical
*              harmonic gaussians to cartesian gaussians.
*
               ijklab = nijkl * iCmp(1)*iCmp(2)
               nW2 = ijklab*Max(kCmpc*lCmpd,mcd)
               iW3_ = iW2 + nW2
               nWrk3_ = nWrk2 - ((iW2-iW4) + nW2)
               Call SphCr1(Wrk2(iW2),ijklab,
     &                     Wrk2(iW3_),nWrk3_,
     &                     RSph(ipSph(lc)),nElem(lc),kCmpc,
     &                     Shells(kShllc)%Transf,
     &                     Shells(kShllc)%Prjct,
     &                     RSph(ipSph(ld)),nElem(ld),lCmpd,
     &                     Shells(lShlld)%Transf,
     &                     Shells(lShlld)%Prjct,
     &                     Wrk2(iW2),mcd)
               If (iW2.eq.iW4) Then
                  nW2 = nijkl*mcd*Max(iCmpa*jCmpb,mab)
                  nW4 = 0
               Else
                  nW2 = nijkl*mcd*iCmpa*jCmpb
                  nW4 = nijkl*mcd*mab
               End If
               iW3_=iW2 + nW2
               nWrk3_ = nWrk2 - (nW2 + nW4)
               Call SphCr2(Wrk2(iW2),nijkl,mcd,
     &                     Wrk2(iW3_),nWrk3_,
     &                     RSph(ipSph(la)),nElem(la),iCmpa,
     &                     Shells(iShlla)%Transf,
     &                     Shells(iShlla)%Prjct,
     &                     RSph(ipSph(lb)),nElem(lb),jCmpb,
     &                     Shells(jShllb)%Transf,
     &                     Shells(jShllb)%Prjct,
     &                     Wrk2(iW4),mab)
*
*--------------Transpose the 2nd order density matrix
*
               If (mab*mcd.ne.1) Then
                  iW3_ = iW4 + nijkl*mab*mcd
                  Call DGetMO(Wrk2(iW4),nijkl,nijkl,mab*mcd,Wrk2(iW3_),
     &                        mab*mcd)
                  call dcopy_(mab*mcd*nijkl,Wrk2(iW3_),1,Wrk2(iW4),1)
               End If
*
               lDCR1=NrOpr(iDCRR(lDCRR))+1
               lDCR2=NrOpr(iDCRS(lDCRS))+1
               ix1 = 1
               iy1 = 1
               iz1 = 1
               ix2 = iPhase(1,iDCRT(lDCRT))
               iy2 = iPhase(2,iDCRT(lDCRT))
               iz2 = iPhase(3,iDCRT(lDCRT))
*
               Call C_F_Pointer(C_Loc(Data1(ip_IndZ(1,nZeta),lDCR1+1)),
     &                         iData1,[nZeta+1])
               Call C_F_Pointer(C_Loc(Data2(ip_IndZ(1,nEta ),lDCR2+1)),
     &                         iData2,[nEta +1])
               nZeta_Tot=iData1(nZeta+1)
               nEta_Tot =iData2(nEta +1)
*
*--------------Loops to partion the primitives
*
               Do 400 iZeta = 1, nZeta_Tot, IncZet
                  mZeta=Min(IncZet,nZeta_Tot-iZeta+1)
                  If (lEmpty(Coeff2,nBeta, nBeta ,jBasj)) Go To 400
*
               Do 410 iEta  = 1, nEta_Tot,  IncEta
                  mEta=Min(IncEta,nEta_Tot-iEta+1)
                  If (lEmpty(Coeff4,nDelta,nDelta,lBasl)) Go To 410
*
                  Pren = Pren + DBLE(mab*mcd*mZeta*mEta)
*
*-----------------Preprescreen
*
                  Call PrePre_g(nZeta,nEta,mZeta,mEta,lZeta,lEta,
     &                          Data1(ip_Z   (iZeta,nZeta),lDCR1),
     &                          Data2(ip_Z   (iEta, nEta), lDCR2),
     &                          PreScr,CutGrd)
                  If (lZeta*lEta.eq.0) Go To 410
*
*-----------------Decontract the 2nd order density matrix
*
                  If (iW4.eq.iW2) Then
                     nW4=0
                     nW2=Max(nijkl,mZeta*mEta)*mab*mcd
                  Else
                     nW4=nijkl*mab*mcd
                     nW2=mZeta*mEta*mab*mcd
                  End If
                  iW3_ = iW2 + nW2
                  nWrk3_ = nWrk2 - (nW4+nW2)
                  Call Tcrtnc(Coeff1,nAlpha,iBasi,
     &                        Coeff2,nBeta,jBasj,
     &                        Coeff3,nGamma,kBask,
     &                        Coeff4,nDelta,lBasl,
     &                        Wrk2(iW4),mab*mcd,Wrk2(iW3_),nWrk3_,
     &                        Wrk2(iW2),
     &                        iData1(iZeta),mZeta,
     &                        iData2(iEta ),mEta)
*
*-----------------Transfer k2 data and prescreen
*
                  iW3_ = iW2 + mZeta*mEta*mab*mcd
                  nWrk3_ = nWrk2 - mZeta*mEta*mab*mcd
                  Call Screen_g(Wrk2(iW2),Wrk2(iW3_),mab*mcd,nZeta,nEta,
     &                        mZeta,mEta,lZeta,lEta,
     &                        Zeta,ZInv,P,xA,xB,
     &                        Data1(iZeta,lDCR1),
     &                        nAlpha,jPrim,
     &                        iData1(iZeta),
     &                        Eta, EInv,Q,xG,xD,
     &                        Data2(iEta ,lDCR2),
     &                        nGamma,lPrim,
     &                        iData2(iEta ),
     &                        ix1,iy1,iz1,ix2,iy2,iz2,
     &                        CutGrd,l2DI,
     &                        Data1(iZeta+iffab ,lDCR1),
     &                        Data1(iZeta+iffabG,lDCR1),nab,
     &                        Data2(iEta +iffcd ,lDCR2),
     &                        Data2(iEta +iffcdG,lDCR2),ncd,
     &                        PreScr,nWrk3_,IsChi,ChiI2)
                  Prem = Prem + DBLE(mab*mcd*lZeta*lEta)
c                 Write (*,*) 'Prem=',Prem
                  If (lZeta*lEta.eq.0) Go To 410
*
*-----------------Compute integral derivative and accumulate
*                 contribution to the molecular gradient. Note that
*                 the PSO matrix now is stored in Wrk2(iW2).
*
                  iW3_ = iW2 + lZeta*lEta*mab*mcd
                  Call Rysg1(iAnga,nRys,lZeta*lEta,
     &                       xA,xB,xG,xD,
     &                       Zeta,ZInv,lZeta,
     &                       Eta,EInv,lEta,
     &                       P,nZeta,Q,nEta,
     &                       CoorM,CoorM,CoorAC,Wrk2(iW3_),nWrk3_,
     &                       TERI1,ModU2,vCff2D,
     &                       Wrk2(iW2),mab*mcd,
     &                       Grad,nGrad,JfGrad,JndGrd,kOp,iuvwx)
                  Aha=Sqrt(DDot_(nGrad,Grad,1,Grad,1))
                  If (Aha.gt.1.0D+5) Then
                      Write (6,*)
     &                      'Norm of gradient contribution is huge!'
                      Write (6,*) 'Probably due to wrong coordinates.'
                  End If
*
 410           Continue
 400           Continue
               Nullify(iData1,iData2)
*
#ifdef _DEBUG_
               If (iPrint.ge.19) Call PrGrad(' In TwoEl',
     &                     Grad,nGrad,ChDisp,5)
#endif
*
 301          Continue
 300        Continue
 200     Continue
 100  Continue
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(iAO)
         Call Unused_integer(iPrInc)
         Call Unused_integer(kPrInc)
      End If
      End
