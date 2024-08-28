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
subroutine k2Loop(Coor,iAnga,iCmpa,iShll,iDCRR,nDCRR,k2data,Alpha,nAlpha,Beta,nBeta,Alpha_,Beta_,Coeff1,iBasn,Coeff2,jBasn,Zeta, &
                  ZInv,Kappab,P,IndP,nZeta,IncZZ,Con,Wrk,nWork2,Cmpct,nScree,mScree,iStb,jStb,Dij,nDij,nDCR,ijCmp,DoFock,Scr,nScr, &
                  Knew,Lnew,Pnew,Qnew,nNew,DoGrad,HMtrx,nHrrMtrx)
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
use Definitions, only: wp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
external Cmpct
integer nZeta, ijCmp, nDCRR, nAlpha, iBasn, nBeta, jBasn, nWork2, nScree, mScree, iStb, jStb, nDij, nDCR, nScr, nNew, nHRRMtrx, &
        IncZZ
type(k2_type), intent(inout) :: k2data(nDCRR)
real*8 Coor(3,4), Alpha(nAlpha), Beta(nBeta), Alpha_(nZeta), Beta_(nZeta), Coeff1(nAlpha,iBasn), Coeff2(nBeta,jBasn), Zeta(nZeta), &
       ZInv(nZeta), Kappab(nZeta), P(nZeta,3), Con(nZeta), Wrk(nWork2), Dij(nDij,nDCR), Scr(nScr,3), Knew(nNew), Lnew(nNew), &
       Pnew(nNew*3), Qnew(nNew*3), HMtrx(nHrrMtrx,2)
logical DoFock, DoGrad
integer iAnga(4), iCmpa(4), iShll(2), iDCRR(0:7), IndP(nZeta)
external TERIS, ModU2, Cff2DS, Rys2D
real*8 CoorM(3,4), Coori(3,4), Coora(3,4), CoorAC(3,2), Q(3), TA(3), TB(3)
logical AeqB, NoSpecial
logical, external :: EQ, TF
integer mStb(2), la, lb, iSmAng, mabMin, mabMax, ne, iCmpa_, jCmpb_, iShlla, jShllb, mcdMin, mcdMax, mabcd, mZeta, nT, iw3, i_Int, &
        iw2, Jnd, iOffZ, lZeta, nDisp, iCmp, ixyz, nabSz, lDCRR, iIrrep, iZeta, iCnt, iComp
real*8 Dummy(1), Tst, ZtMax, abMax, ZtMaxD, abMaxD, Tmp, Delta, TEMP
real*8, external :: EstI
! Statement function to compute canonical index
nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6-1

!                                                                      *
!***********************************************************************
!                                                                      *
call dcopy_(3,[One],0,Q,1)
mStb(1) = iStb
mStb(2) = jStb
la = iAnga(1)
lb = iAnga(2)
iSmAng = la+lb+la+lb
iCmpa_ = iCmpa(1)
jCmpb_ = iCmpa(2)
iShlla = iShll(1)
jShllb = iShll(2)
!                                                                      *
!***********************************************************************
!                                                                      *
call dcopy_(3,Coor(1,1),1,CoorM(1,1),1)
!                                                                      *
!***********************************************************************
!                                                                      *
do lDCRR=0,nDCRR-1

  call OA(iDCRR(lDCRR),Coor(1:3,2),CoorM(1:3,2))
  AeqB = EQ(CoorM(1,1),CoorM(1,2))
  ! Branch out if integrals are zero by symmetry.
  if (AeqB .and. (mod(iSmAng,2) == 1)) Go To 100
  call dcopy_(6,CoorM(1,1),1,CoorM(1,3),1)
# ifdef _DEBUGPRINT_
  call RecPrt(' Actual centers',' ',CoorM,3,4)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute zeta, P and kappa.

  ! No triangulatization applied at this level
  call DoZeta(Alpha,nAlpha,Beta,nBeta,CoorM(1,1),CoorM(1,2),P,Zeta,Kappab,ZInv,Alpha_,Beta_,IndP)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Generate transformation matrix from intermediate integrals
  ! to final angular composition.

  mabMin = nabSz(max(la,lb)-1)+1
  if (EQ(CoorM(1,1),CoorM(1,2))) mabMin = nabSz(la+lb-1)+1
  mabMax = nabSz(la+lb)
  ne = (mabMax-mabMin+1)
  do iIrrep=0,nIrrep-1
    call OA(iOper(iIrrep),CoorM(1:3,1),TA)
    call OA(iOper(iIrrep),CoorM(1:3,2),TB)
    call HrrMtrx(k2Data(lDCRR+1)%HrrMtrx(:,iIrrep+1),ne,la,lb,TA,TB,Shells(iShlla)%Transf,RSph(ipSph(la)),iCmpa_, &
                 Shells(jShllb)%Transf,RSph(ipSph(lb)),jCmpb_)
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute primitive integrals to be used in the prescreening
  ! by the Schwartz inequality.

  call dcopy_(12,CoorM(1,1),1,Coora(1,1),1)
  call dcopy_(12,CoorM(1,1),1,Coori(1,1),1)

  ! Compute actual size of [a0|c0] block

  mcdMin = mabMin
  mcdMax = mabMax
  mabcd = (mabMax-mabMin+1)*(mcdMax-mcdMin+1)

  ! Find the proper centers to start of with the angular
  ! momentum on. If la == lb there will exist an
  ! ambiguity to which center that angular momentum should
  ! be accumulated on. In that case we will use A and C of
  ! the order as defined by the basis functions types.

  if (iAnga(1) >= iAnga(2)) then
    call dcopy_(3,Coora(1,1),1,CoorAC(1,1),1)
  else
    call dcopy_(3,Coora(1,2),1,CoorAC(1,1),1)
  end if
  call dcopy_(3,CoorAC(1,1),1,CoorAC(1,2),1)

  ! Compute [a0|c0], ijkl,a,c

  Jnd = 0
  nScree = nScree+nZeta
  do iZeta=1,nZeta,IncZZ
    mZeta = min(nZeta-iZeta+1,IncZZ)

    nT = mZeta*1
    NoSpecial = .true.
    call Rys(iAnga,nT,Zeta(iZeta),ZInv(iZeta),mZeta,[One],[One],1,P(iZeta,1),nZeta,Q,1,Kappab(iZeta),[One],Coori,Coora,CoorAC, &
             mabMin,mabMax,mcdMin,mcdMax,Wrk,nWork2,TERIS,ModU2,Cff2DS,Rys2D,NoSpecial)
#   ifdef _DEBUGPRINT_
    call RecPrt(' In k2Loop: ijkl,[a0|c0]',' ',Wrk,mZeta,mabcd)
#   endif

    ! Apply a transpose prior to Tnsctl to fake the action of Cntrct.

    iW3 = 1+mZeta*mabcd
    call DGeTMO(Wrk,mZeta,mZeta,mabcd,Wrk(iW3),mabcd)
    call dcopy_(mabcd*mZeta,Wrk(iW3),1,Wrk,1)
    call TnsCtl(Wrk,nWork2,mZeta,mabMax,mabMin,mabMax,mabMin,k2data(lDCRR+1)%HrrMtrx(:,1),k2data(lDCRR+1)%HrrMtrx(:,1),la,lb,la, &
                lb,iCmpa_,jCmpb_,iCmpa_,jCmpb_,iShlla,jShllb,iShlla,jShllb,i_Int)
    if (i_Int == 1) then
      iW2 = 1
      iW3 = 1+mZeta*(iCmpa_*jCmpb_)**2
    else
      iW2 = i_Int
      iW3 = 1
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Store data in core

    call Cmpct(Wrk(iW2),iCmpa_,jCmpb_,nZeta,mZeta,Zeta(iZeta),Kappab(iZeta),P(iZeta,1),IndP(iZeta),Con,k2Data(lDCRR+1)%Zeta(:), &
               k2Data(lDCRR+1)%Kappa(:),k2Data(lDCRR+1)%PCoor(:,:),k2Data(lDCRR+1)%IndZ(:),iZeta-1,Jnd,k2Data(lDCRR+1)%ZInv(:), &
               AeqB,k2Data(lDCRR+1)%ab(:),k2Data(lDCRR+1)%abCon(:),Alpha_(iZeta),k2Data(lDCRR+1)%Alpha(:),Beta_(iZeta), &
               k2Data(lDCRR+1)%Beta(:))

  end do ! iZeta
  mScree = mScree+Jnd
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Estimate the largest contracted integral.

  k2Data(lDCRR+1)%EstI = EstI(k2Data(lDCRR+1)%Zeta(:),k2Data(lDCRR+1)%Kappa(:),nAlpha,nBeta,Coeff1,iBasn,Coeff2,jBasn, &
                              k2Data(lDCRR+1)%ab(:),iCmpa_*jCmpb_,Wrk,nWork2,k2Data(lDCRR+1)%IndZ(:))
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Find the largest integral estimate (AO Basis).

  Tst = -One
  do iZeta=1,nZeta
    Tst = max(k2Data(lDCRR+1)%Zeta(iZeta),Tst)
  end do
  k2Data(lDCRR+1)%ZetaM = tst

  iOffZ = nDij-nZeta-1
  ZtMax = One
  abMax = Zero
  ZtMaxD = One
  abMaxD = Zero
  do iZeta=0,Jnd-1
    tmp = k2Data(lDCRR+1)%abCon(iZeta+1)
    if (abMax < tmp) then
      abMax = tmp
      ZtMax = k2Data(lDCRR+1)%Zeta(iZeta+1)
    end if
    if (DoFock) then
      tmp = k2Data(lDCRR+1)%ab(iZeta+1)*Dij(iOffZ+iZeta,lDCRR+1)
      if (abMaxD < tmp) then
        abMaxD = tmp
        ZtMaxD = k2Data(lDCRR+1)%Zeta(iZeta+1)
      end if
    else
      ZtMaxD = -One
      abMaxD = Zero
    end if
  end do
  k2Data(lDCRR+1)%ZtMax = ZtMax
  k2Data(lDCRR+1)%abMax = abMax
  k2Data(lDCRR+1)%ZtMaxD = ZtMaxD
  k2Data(lDCRR+1)%abMaxD = abMaxD
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute data for gradient evaluation

  mZeta = Jnd
  if (DoGrad .and. (mZeta > 0)) then

    ! Compute primitive integrals to be used in the prescreening
    ! by the Cauchy-Schwarz inequality.

    do iZeta=1,mZeta,IncZZ
      lZeta = min(mZeta-iZeta+1,IncZZ)
      call SchInt(CoorM,iAnga,iCmpa,lZeta,k2Data(lDCRR+1)%Zeta(iZeta:),k2Data(lDCRR+1)%ZInv(iZeta:),k2Data(lDCRR+1)%Kappa(iZeta:), &
                  k2Data(lDCRR+1)%PCoor(iZeta:,:),k2Data(lDCRR+1)%Kappa(iZeta:),k2Data(lDCRR+1)%PCoor(iZeta:,:),nZeta,Wrk,nWork2, &
                  HMtrx,nHrrMtrx,iShlla,jShllb,i_Int)
      call PckInt(Wrk(i_Int),lZeta,ijCmp,k2Data(lDCRR+1)%abG(iZeta:,1),k2Data(lDCRR+1)%Kappa(iZeta:),.true., &
                  k2Data(lDCRR+1)%Zeta(iZeta:),nZeta,Dummy)
    end do

    ! Second order numerical differentiation. The gradients are
    ! restricted to only those with respect to symmetrical
    ! displacements. The symmetric three point formula is used in
    ! the numerical procedure.

    iIrrep = 0
    Delta = 1.0e-3_wp
    k2Data(lDCRR+1)%abG(:,2) = Zero
    Scr(1:nZeta*ijCmp,:) = Zero

    ! Loop over center A and B.

    do iCnt=1,2

      nDisp = IndDsp(mStb(iCnt),iIrrep)
      do iComp=1,3
        iCmp = 2**(iComp-1)
        if (TF(mStb(iCnt),iIrrep,iCmp) .and. Dirct(nDisp+1)) then
          nDisp = nDisp+1
          temp = CoorM(iComp,iCnt)

          CoorM(iComp,iCnt) = temp+Delta
          CoorM(iComp,iCnt+2) = temp+Delta
          call NewPK(CoorM(1,1),CoorM(1,2),Pnew,mZeta,nZeta,Knew,k2Data(lDCRR+1)%Alpha(:),k2Data(lDCRR+1)%Beta(:))
          do iZeta=1,mZeta,IncZZ
            lZeta = min(mZeta-iZeta+1,IncZZ)
            call SchInt(CoorM,iAnga,iCmpa,lZeta,k2Data(lDCRR+1)%Zeta(iZeta:),k2Data(lDCRR+1)%ZInv(iZeta:),Knew(iZeta),Pnew(iZeta), &
                        Knew(iZeta),Pnew(iZeta),nZeta,Wrk,nWork2,HMtrx,nHrrMtrx,iShlla,jShllb,i_Int)
            call PckInt(Wrk(i_Int),lZeta,ijCmp,Scr(iZeta,1),Knew(iZeta),.false.,k2Data(lDCRR+1)%Zeta(iZeta:),nZeta,Knew(iZeta))
          end do

          CoorM(iComp,iCnt) = temp-Delta
          CoorM(iComp,iCnt+2) = temp-Delta
          call NewPK(CoorM(1,1),CoorM(1,2),Qnew,mZeta,nZeta,Lnew,k2Data(lDCRR+1)%Alpha(:),k2Data(lDCRR+1)%Beta(:))
          do iZeta=1,mZeta,IncZZ
            lZeta = min(mZeta-iZeta+1,IncZZ)
            call SchInt(CoorM,iAnga,iCmpa,lZeta,k2Data(lDCRR+1)%Zeta(iZeta:),k2Data(lDCRR+1)%ZInv(iZeta:),Lnew(iZeta),Qnew(iZeta), &
                        Lnew(iZeta),Qnew(iZeta),nZeta,Wrk,nWork2,HMtrx,nHrrMtrx,iShlla,jShllb,i_Int)
            call PckInt(Wrk(i_Int),lZeta,ijCmp,Scr(iZeta,2),Lnew(iZeta),.false.,k2Data(lDCRR+1)%Zeta(iZeta:),nZeta,Lnew(iZeta))
          end do

          call DaXpY_(nZeta*ijCmp,One,Scr(1,2),1,Scr(1,1),1)

          CoorM(iComp,iCnt) = temp+Delta
          CoorM(iComp,iCnt+2) = temp-Delta
          do iZeta=1,mZeta,IncZZ
            lZeta = min(mZeta-iZeta+1,IncZZ)
            call SchInt(CoorM,iAnga,iCmpa,lZeta,k2Data(lDCRR+1)%Zeta(iZeta:),k2Data(lDCRR+1)%ZInv(iZeta:),Knew(iZeta),Pnew(iZeta), &
                        Lnew(iZeta),Qnew(iZeta),nZeta,Wrk,nWork2,HMtrx,nHrrMtrx,iShlla,jShllb,i_Int)
            call PckInt(Wrk(i_Int),lZeta,ijCmp,Scr(iZeta,3),Knew(iZeta),.false.,k2Data(lDCRR+1)%Zeta(iZeta:),nZeta,Lnew(iZeta))
          end do

          call DaXpY_(nZeta*ijCmp,-One,Scr(1,3),1,Scr(1,1),1)

          CoorM(iComp,iCnt) = temp-Delta
          CoorM(iComp,iCnt+2) = temp+Delta
          do iZeta=1,mZeta,IncZZ
            lZeta = min(mZeta-iZeta+1,IncZZ)
            call SchInt(CoorM,iAnga,iCmpa,lZeta,k2Data(lDCRR+1)%Zeta(iZeta:),k2Data(lDCRR+1)%ZInv(iZeta:),Lnew(iZeta),Qnew(iZeta), &
                        Knew(iZeta),Pnew(iZeta),nZeta,Wrk,nWork2,HMtrx,nHrrMtrx,iShlla,jShllb,i_Int)
            call PckInt(Wrk(i_Int),lZeta,ijCmp,Scr(iZeta,3),Lnew(iZeta),.false.,k2Data(lDCRR+1)%Zeta(iZeta:),nZeta,Knew(iZeta))
          end do

          call DaXpY_(nZeta*ijCmp,-One,Scr(1,3),1,Scr(1,1),1)

          call DScal_(nZeta*ijCmp,One/(Four*Delta**2),Scr(1,1),1)
          call AbsAdd(nZeta*ijCmp,Scr(:,1),1,k2data(lDCRR+1)%abG(:,2),1)

          CoorM(iComp,iCnt) = temp
          CoorM(iComp,iCnt+2) = temp
        end if
      end do

    end do
  end if       ! DoGrad
  !                                                                    *
  !*********************************************************************
  !                                                                    *
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'lDCRR=',lDCRR
  call WrCheck('Zeta ',k2Data(lDCRR+1)%Zeta(:),nZeta)
  call WrCheck('Kappa',k2Data(lDCRR+1)%Kappa(:),nZeta)
  call WrCheck('P    ',K2Data(lDCRR+1)%PCoor(:,:),nZeta*3)
  call WrCheck('xA   ',k2Data(lDCRR+1)%Alpha(:),nZeta)
  call WrCheck('xB   ',k2Data(lDCRR+1)%Beta(:),nZeta)
  call WrCheck('ZInv ',k2Data(lDCRR+1)%ZInv(:),nZeta)
  if (DoGrad) then
    call WrCheck('ab   ',k2data(lDCRR+1)%abG(:,1),nZeta*ijCmp)
    call WrCheck('abG  ',k2data(lDCRR+1)%abG(:,2),nZeta*ijCmp)
  end if
  write(u6,*)
  write(u6,*) ' ERI(Max)=',k2Data(lDCRR+1)%EstI
  write(u6,*) ' ZtMax   =',k2Data(lDCRR+1)%ZtMax
  write(u6,*) ' abMax   =',k2Data(lDCRR+1)%abMax
  write(u6,*) ' ZtMaxD  =',k2Data(lDCRR+1)%ZtMaxD
  write(u6,*) ' abMaxD  =',k2Data(lDCRR+1)%abMaxD
  call WrCheck(' HrrMtrx',k2Data(lDCRR+1)%HrrMtrx(:,:),ne*iCmpa_*jCmpb_)
# endif
100 continue ! lDCRR
end do

return

end subroutine k2Loop
