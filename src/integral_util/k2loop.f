!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1990,1991,1993,1999, Roland Lindh                      *
!               1990, IBM                                              *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine k2Loop(Coor,
     &                  iAnga,iCmpa,iShll,
     &                  iDCRR,nDCRR,
     &                  k2data,
     &                  Alpha,nAlpha,Beta, nBeta,
     &                  Alpha_,Beta_,
     &                  Coeff1,iBasn,Coeff2,jBasn,
     &                  Zeta,ZInv,Kappab,P,IndP,nZeta,IncZZ,Con,
     &                  Wrk,nWork2,
     &                  Cmpct,nScree,mScree,iStb,jStb,
     &                  Dij,nDij,nDCR,ijCmp,DoFock,
     &                  Scr,nScr,
     &                  Knew,Lnew,Pnew,Qnew,nNew,DoGrad,HMtrx,nHrrMtrx)
!***********************************************************************
!                                                                      *
! Object: to compute zeta, kappa, P, and the integrals [nm|nm] for     *
!         prescreening. This is done for all unique pairs of centers   *
!         generated from the symmetry unique centers A and B.          *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN.                              *
!             June '91, modified to compute zeta, P, kappa and inte-   *
!             grals for Schwartz inequality in a k2 loop.              *
!             Modified for direct SCF, January '93                     *
!***********************************************************************
      use Real_Spherical, only: ipSph, rSph
      use Basis_Info, only: Shells
      use Symmetry_Info, only: nIrrep, iOper
      use Constants, only: Zero, One, Four
      use Disp, only: Dirct, IndDsp
      use k2_structure, only: k2_type
      Implicit None
      External Cmpct
      Integer nZeta, ijCmp,  nDCRR,
     &        nAlpha, iBasn, nBeta, jBasn, nWork2, nScree, mScree,
     &        iStb, jStb, nDij, nDCR, nScr, nNew, nHRRMtrx, IncZZ
      type(k2_type), intent(inout) :: k2data(nDCRR)
      Real*8 Coor(3,4),
     &       Alpha(nAlpha), Beta(nBeta), Alpha_(nZeta), Beta_(nZeta),
     &       Coeff1(nAlpha,iBasn), Coeff2(nBeta,jBasn),
     &       Zeta(nZeta), ZInv(nZeta), Kappab(nZeta), P(nZeta,3),
     &       Con(nZeta), Wrk(nWork2), Dij(nDij,nDCR), Scr(nScr,3),
     &       Knew(nNew), Lnew(nNew), Pnew(nNew*3), Qnew(nNew*3),
     &       HMtrx(nHrrMtrx,2)
      Logical DoFock, DoGrad
      Integer iAnga(4), iCmpa(4), iShll(2), iDCRR(0:7), IndP(nZeta)

      External TERIS, ModU2, Cff2DS, Rys2D
      Real*8 CoorM(3,4), Coori(3,4), Coora(3,4), CoorAC(3,2),
     &       Q(3), TA(3), TB(3)
      Logical AeqB, NoSpecial
      Logical, External:: EQ, TF
      Integer  mStb(2), la, lb, iSmAng, mabMin, mabMax, ne,
     &         iCmpa_, jCmpb_, iShlla, jShllb,
     &         mcdMin, mcdMax, mabcd, mZeta, nT,
     &         iw3, i_Int, iw2, Jnd, iOffZ, lZeta, nDisp,
     &         iCmp, ixyz, nabSz, lDCRR, iIrrep, iZeta, iCnt, iComp
      Real*8 Dummy(1), Tst, ZtMax, abMax, ZtMaxD, abMaxD, Tmp, Delta,
     &       TEMP
      Real*8, External :: EstI
!                                                                      *
!***********************************************************************
!                                                                      *
!     Statement function to compute canonical index
!
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
!                                                                      *
!***********************************************************************
!                                                                      *
      call dcopy_(3,[One],0,Q,1)
      mStb(1) = iStb
      mStb(2) = jStb
      la = iAnga(1)
      lb = iAnga(2)
      iSmAng=la+lb+la+lb
      iCmpa_= iCmpa(1)
      jCmpb_= iCmpa(2)
      iShlla = iShll(1)
      jShllb = iShll(2)
!                                                                      *
!***********************************************************************
!                                                                      *
      call dcopy_(3,Coor(1,1),1,CoorM(1,1),1)
!                                                                      *
!***********************************************************************
!                                                                      *
      Do 100 lDCRR = 0, nDCRR-1
!
         Call OA(iDCRR(lDCRR),Coor(1:3,2),CoorM(1:3,2))
         AeqB = EQ(CoorM(1,1),CoorM(1,2))
!        Branch out if integrals are zero by symmetry.
         If (AeqB .and. Mod(iSmAng,2).eq.1) Go To 100
         call dcopy_(6,CoorM(1,1),1,CoorM(1,3),1)
#ifdef _DEBUGPRINT_
         Call RecPrt(' Actual centers',' ',CoorM,3,4)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!        Compute zeta, P and kappa.
!
!        No triangulatization applied at this level
         Call DoZeta(Alpha,nAlpha,Beta,nBeta,
     &               CoorM(1,1),CoorM(1,2),
     &               P,
     &               Zeta,
     &               Kappab,
     &               ZInv,
     &               Alpha_,
     &               Beta_,
     &               IndP)
!                                                                      *
!***********************************************************************
!                                                                      *
!        Generate transformation matrix from intermediate integrals
!        to final angular composition.
!
         mabMin=nabSz(Max(la,lb)-1)+1
         If (EQ(CoorM(1,1),CoorM(1,2))) mabMin = nabSz(la+lb-1)+1
         mabMax=nabSz(la+lb)
         ne=(mabMax-mabMin+1)
         Do iIrrep = 0, nIrrep-1
            Call OA(iOper(iIrrep),CoorM(1:3,1),TA)
            Call OA(iOper(iIrrep),CoorM(1:3,2),TB)
            Call HrrMtrx(k2Data(lDCRR+1)%HrrMtrx(:,iIrrep+1),
     &                   ne,la,lb,TA,TB,
     &                   Shells(iShlla)%Transf,RSph(ipSph(la)),iCmpa_,
     &                   Shells(jShllb)%Transf,RSph(ipSph(lb)),jCmpb_)
         End Do
!                                                                      *
!***********************************************************************
!                                                                      *
!        Compute primitive integrals to be used in the prescreening
!        by the Schwartz inequality.
!
         call dcopy_(12,CoorM(1,1),1,Coora(1,1),1)
         call dcopy_(12,CoorM(1,1),1,Coori(1,1),1)
!
!        Compute actual size of [a0|c0] block
!
         mcdMin=mabMin
         mcdMax=mabMax
         mabcd=(mabMax-mabMin+1)*(mcdMax-mcdMin+1)
!
!        Find the proper centers to start of with the angular
!        momentum on. If la.eq.lb there will exist an
!        ambiguity to which center that angular momentum should
!        be accumulated on. In that case we will use A and C of
!        the order as defined by the basis functions types.
!
         If (iAnga(1).ge.iAnga(2)) Then
            call dcopy_(3,Coora(1,1),1,CoorAC(1,1),1)
         Else
            call dcopy_(3,Coora(1,2),1,CoorAC(1,1),1)
         End If
         call dcopy_(3,CoorAC(1,1),1,CoorAC(1,2),1)
!
!        Compute [a0|c0], ijkl,a,c
!
         Jnd = 0
         nScree = nScree + nZeta
         Do iZeta = 1, nZeta, IncZZ
            mZeta = Min(nZeta-iZeta+1,IncZZ)
!
            nT = mZeta*1
            NoSpecial=.True.
            Call Rys(iAnga,nT,
     &               Zeta(iZeta),ZInv(iZeta),mZeta,[One],[One],1,
     &               P(iZeta,1),nZeta,Q,1,Kappab(iZeta),[One],
     &               Coori,Coora,CoorAC,
     &               mabMin,mabMax,mcdMin,mcdMax,
     &               Wrk,nWork2,TERIS,ModU2,Cff2DS,
     &               Rys2D,NoSpecial)
#ifdef _DEBUGPRINT_
            Call RecPrt(' In k2Loop: ijkl,[a0|c0]',' ',Wrk,
     &                            mZeta,mabcd)
#endif
!
!---------- Apply a transpose prior to Tnsctl to fake the action
!           of Cntrct.
!
            iW3=1+mZeta*mabcd
            Call DGeTMO(Wrk,mZeta,mZeta,mabcd,Wrk(iW3),mabcd)
            call dcopy_(mabcd*mZeta,Wrk(iW3),1,Wrk,1)
            Call TnsCtl(Wrk,nWork2,mZeta,mabMax,mabMin,mabMax,mabMin,
     &                  k2data(lDCRR+1)%HrrMtrx(:,1),
     &                  k2data(lDCRR+1)%HrrMtrx(:,1),
     &                  la,lb,la,lb,
     &                  iCmpa_,jCmpb_,iCmpa_,jCmpb_,
     &                  iShlla,jShllb,iShlla,jShllb,i_Int)
            If (i_Int.eq.1) Then
               iW2=1
               iW3=1+mZeta*(iCmpa_*jCmpb_)**2
            Else
               iW2=i_Int
               iW3=1
            End If
!                                                                      *
!***********************************************************************
!                                                                      *
!-----------Store data in core
!
            Call Cmpct(Wrk(iW2),iCmpa_,jCmpb_,nZeta,mZeta,
     &                 Zeta(iZeta),Kappab(iZeta),
     &                 P(iZeta,1),IndP(iZeta),Con,
     &                 k2Data(lDCRR+1)%Zeta(:),
     &                 k2Data(lDCRR+1)%Kappa(:),
     &                 k2Data(lDCRR+1)%PCoor(:,:),
     &                 k2Data(lDCRR+1)%IndZ(:),iZeta-1,Jnd,
     &                 k2Data(lDCRR+1)%ZInv(:),
     &                 AeqB,
     &                 k2Data(lDCRR+1)%ab(:),
     &                 k2Data(lDCRR+1)%abCon(:),
     &                 Alpha_(iZeta),
     &                 k2Data(lDCRR+1)%Alpha(:),
     &                 Beta_(iZeta),
     &                 k2Data(lDCRR+1)%Beta(:))
!
         End Do ! iZeta
         mScree = mScree + Jnd
!                                                                      *
!***********************************************************************
!                                                                      *
!        Estimate the largest contracted integral.
!
         k2Data(lDCRR+1)%EstI=
     &                      EstI(k2Data(lDCRR+1)%Zeta(:),
     &                           k2Data(lDCRR+1)%Kappa(:),
     &                           nAlpha,nBeta,
     &                           Coeff1,iBasn,Coeff2,jBasn,
     &                           k2Data(lDCRR+1)%ab(:),
     &                           iCmpa_*jCmpb_,
     &                           Wrk,nWork2,
     &                           k2Data(lDCRR+1)%IndZ(:))
!                                                                      *
!***********************************************************************
!                                                                      *
!------- Find the largest integral estimate (AO Basis).
!
         Tst  = -One
         Do  iZeta = 1, nZeta
             Tst=Max(k2Data(lDCRR+1)%Zeta(iZeta),Tst)
         End Do
         k2Data(lDCRR+1)%ZetaM = tst
!
         iOffZ = nDij-nZeta-1
         ZtMax=One
         abMax=Zero
         ZtMaxD=One
         abMaxD=Zero
         Do iZeta = 0, Jnd-1
            tmp = k2Data(lDCRR+1)%abCon(iZeta+1)
            If (abMax.lt.tmp) Then
               abMax = tmp
               ZtMax = k2Data(lDCRR+1)%Zeta(iZeta+1)
            End If
            If (DoFock) Then
               tmp = k2Data(lDCRR+1)%ab(iZeta+1)
     &             * Dij(iOffZ+iZeta,lDCRR+1)
               If (abMaxD.lt.tmp) Then
                 abMaxD = tmp
                 ZtMaxD = k2Data(lDCRR+1)%Zeta(iZeta+1)
               End If
            Else
                 ZtMaxD=-One
                 abMaxD=Zero
            End If
         End Do
         k2Data(lDCRR+1)%ZtMax = ZtMax
         k2Data(lDCRR+1)%abMax = abMax
         k2Data(lDCRR+1)%ZtMaxD= ZtMaxD
         k2Data(lDCRR+1)%abMaxD= abMaxD
!                                                                      *
!***********************************************************************
!                                                                      *
!------- Compute data for gradient evaluation
!
         mZeta=Jnd
         If (DoGrad.and.mZeta.gt.0) Then
!
!---------- Compute primitive integrals to be used in the prescreening
!           by the Cauchy-Schwarz inequality.
!
            Do iZeta = 1, mZeta, IncZZ
               lZeta = Min(mZeta-iZeta+1,IncZZ)
               Call SchInt(CoorM,iAnga,iCmpa,lZeta,
     &                     k2Data(lDCRR+1)%Zeta(iZeta:),
     &                     k2Data(lDCRR+1)%ZInv(iZeta:),
     &                     k2Data(lDCRR+1)%Kappa(iZeta:),
     &                     k2Data(lDCRR+1)%PCoor(iZeta:,:),
     &                     k2Data(lDCRR+1)%Kappa(iZeta:),
     &                     k2Data(lDCRR+1)%PCoor(iZeta:,:),
     &                     nZeta,Wrk,nWork2,HMtrx,
     &                     nHrrMtrx,iShlla,jShllb,i_Int)
               Call PckInt(Wrk(i_Int),lZeta,ijCmp,
     &                     k2Data(lDCRR+1)%abG(iZeta:,1),
     &                     k2Data(lDCRR+1)%Kappa(iZeta:),.True.,
     &                     k2Data(lDCRR+1)%Zeta(iZeta:),nZeta,
     &                     Dummy)
            End Do
!
!---------- Second order numerical differentiation. The gradients are
!           restricted to only those with respect to symmetrical
!           displacements. The symmetric three point formula is used in
!           the numerical procedure.
!
            iIrrep=0
            Delta = 1.0D-03
            k2Data(lDCRR+1)%abG(:,2)=Zero
            Scr(1:nZeta*ijCmp,:)=Zero
!
!---------- Loop over center A and B.
!
            Do iCnt = 1, 2
!
               nDisp=IndDsp(mStb(iCnt),iIrrep)
               Do iComp = 1, 3
                  iCmp=2**(iComp-1)
                  If (TF(mStb(iCnt),iIrrep,iCmp) .and.
     &                Dirct(nDisp+1)) Then
                     nDisp = nDisp + 1
                     temp = CoorM(iComp,iCnt)
!
                     CoorM(iComp,iCnt  ) = temp + Delta
                     CoorM(iComp,iCnt+2) = temp + Delta
                     Call NewPK(CoorM(1,1),CoorM(1,2),
     &                          Pnew,mZeta,nZeta,
     &                          Knew,
     &                          k2Data(lDCRR+1)%Alpha(:),
     &                          k2Data(lDCRR+1)%Beta(:))
                     Do iZeta = 1, mZeta, IncZZ
                        lZeta = Min(mZeta-iZeta+1,IncZZ)
                        Call SchInt(CoorM,
     &                              iAnga,iCmpa,lZeta,
     &                              k2Data(lDCRR+1)%Zeta(iZeta:),
     &                              k2Data(lDCRR+1)%ZInv(iZeta:),
     &                              Knew(iZeta),
     &                              Pnew(iZeta),
     &                              Knew(iZeta),
     &                              Pnew(iZeta),
     &                              nZeta,Wrk,nWork2,
     &                              HMtrx,nHrrMtrx,iShlla,jShllb,i_Int)
                        Call PckInt(Wrk(i_Int),lZeta,ijCmp,
     &                              Scr(iZeta,1),
     &                              Knew(iZeta),
     &                              .False.,
     &                              k2Data(lDCRR+1)%Zeta(iZeta:),
     &                              nZeta,
     &                              Knew(iZeta))
                     End Do
!
                     CoorM(iComp,iCnt  ) = temp - Delta
                     CoorM(iComp,iCnt+2) = temp - Delta
                     Call NewPK(CoorM(1,1),CoorM(1,2),
     &                          Qnew,mZeta,nZeta,
     &                          Lnew,
     &                          k2Data(lDCRR+1)%Alpha(:),
     &                          k2Data(lDCRR+1)%Beta(:))
                     Do iZeta = 1, mZeta, IncZZ
                        lZeta = Min(mZeta-iZeta+1,IncZZ)
                        Call SchInt(CoorM,
     &                              iAnga,iCmpa,lZeta,
     &                              k2Data(lDCRR+1)%Zeta(iZeta:),
     &                              k2Data(lDCRR+1)%ZInv(iZeta:),
     &                              Lnew(iZeta),
     &                              Qnew(iZeta),
     &                              Lnew(iZeta),
     &                              Qnew(iZeta),
     &                              nZeta,Wrk,nWork2,
     &                              HMtrx,nHrrMtrx,iShlla,jShllb,i_Int)
                        Call PckInt(Wrk(i_Int),lZeta,ijCmp,
     &                              Scr(iZeta,2),
     &                              Lnew(iZeta),
     &                              .False.,
     &                              k2Data(lDCRR+1)%Zeta(iZeta:),
     &                              nZeta,
     &                              Lnew(iZeta))
                     End Do
!
                     Call DaXpY_(nZeta*ijCmp, One,Scr(1,2),1,
     &                                           Scr(1,1),1)
!
                     CoorM(iComp,iCnt  ) = temp + Delta
                     CoorM(iComp,iCnt+2) = temp - Delta
                     Do iZeta = 1, mZeta, IncZZ
                        lZeta = Min(mZeta-iZeta+1,IncZZ)
                        Call SchInt(CoorM,
     &                              iAnga,iCmpa,lZeta,
     &                              k2Data(lDCRR+1)%Zeta(iZeta:),
     &                              k2Data(lDCRR+1)%ZInv(iZeta:),
     &                              Knew(iZeta),
     &                              Pnew(iZeta),
     &                              Lnew(iZeta),
     &                              Qnew(iZeta),
     &                              nZeta,Wrk,nWork2,
     &                              HMtrx,nHrrMtrx,iShlla,jShllb,i_Int)
                        Call PckInt(Wrk(i_Int),lZeta,ijCmp,
     &                              Scr(iZeta,3),
     &                              Knew(iZeta),
     &                              .False.,
     &                              k2Data(lDCRR+1)%Zeta(iZeta:),
     &                              nZeta,
     &                              Lnew(iZeta))
                     End Do
!
                     Call DaXpY_(nZeta*ijCmp,-One,Scr(1,3),1,
     &                                           Scr(1,1),1)
!
                     CoorM(iComp,iCnt  ) = temp - Delta
                     CoorM(iComp,iCnt+2) = temp + Delta
                     Do iZeta = 1, mZeta, IncZZ
                        lZeta = Min(mZeta-iZeta+1,IncZZ)
                        Call SchInt(CoorM,
     &                              iAnga,iCmpa,lZeta,
     &                              k2Data(lDCRR+1)%Zeta(iZeta:),
     &                              k2Data(lDCRR+1)%ZInv(iZeta:),
     &                              Lnew(iZeta),
     &                              Qnew(iZeta),
     &                              Knew(iZeta),
     &                              Pnew(iZeta),
     &                              nZeta,Wrk,nWork2,
     &                              HMtrx,nHrrMtrx,iShlla,jShllb,i_Int)
                        Call PckInt(Wrk(i_Int),lZeta,ijCmp,
     &                              Scr(iZeta,3),
     &                              Lnew(iZeta),
     &                              .False.,
     &                              k2Data(lDCRR+1)%Zeta(iZeta:),
     &                              nZeta,
     &                              Knew(iZeta))
                     End Do
!
                     Call DaXpY_(nZeta*ijCmp,-One,Scr(1,3),1,
     &                                           Scr(1,1),1)
!
                     Call DScal_(nZeta*ijCmp,One/(Four*Delta**2),
     &                                           Scr(1,1),1)
                     Call AbsAdd(nZeta*ijCmp,    Scr(:,1),1,
     &                                k2data(lDCRR+1)%abG(:,2),1)
!
                     CoorM(iComp,iCnt  ) = temp
                     CoorM(iComp,iCnt+2) = temp
                  End If
               End Do
!
            End Do
         End If       ! DoGrad
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
         Write (6,*)
         Write (6,*) 'lDCRR=',lDCRR
         Call WrCheck('Zeta ',k2Data(lDCRR+1)%Zeta(:),nZeta)
         Call WrCheck('Kappa',k2Data(lDCRR+1)%Kappa(:),nZeta)
         Call WrCheck('P    ',K2Data(lDCRR+1)%PCoor(:,:),nZeta*3)
         Call WrCheck('xA   ',k2Data(lDCRR+1)%Alpha(:),nZeta)
         Call WrCheck('xB   ',k2Data(lDCRR+1)%Beta(:),nZeta)
         Call WrCheck('ZInv ',k2Data(lDCRR+1)%ZInv(:),nZeta)
         If (DoGrad) Then
            Call WrCheck('ab   ',
     &        k2data(lDCRR+1)%abG(:,1),nZeta*ijCmp)
            Call WrCheck('abG  ',
     &        k2data(lDCRR+1)%abG(:,2),nZeta*ijCmp)
         End If
         Write (6,*)
         Write (6,*) ' ERI(Max)=',k2Data(lDCRR+1)%EstI
         Write (6,*) ' ZtMax   =',k2Data(lDCRR+1)%ZtMax
         Write (6,*) ' abMax   =',k2Data(lDCRR+1)%abMax
         Write (6,*) ' ZtMaxD  =',k2Data(lDCRR+1)%ZtMaxD
         Write (6,*) ' abMaxD  =',k2Data(lDCRR+1)%abMaxD
         Call WrCheck(' HrrMtrx',
     &        k2Data(lDCRR+1)%HrrMtrx(:,:),
     &        ne*iCmpa_*jCmpb_)
#endif
 100  Continue ! lDCRR
!
      Return
!
      End SubRoutine k2Loop
