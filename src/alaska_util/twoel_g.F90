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
! Copyright (C) 1990,1992, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine TwoEl_g(Coor,iAnga,iCmp,iShell,iShll,iAO,iStb,jStb,kStb,lStb,nRys, &
                   k2Data1,k2Data2, &
                   nData1,nData2,Pren,Prem,nAlpha,iPrInc,nBeta,jPrInc,nGamma,kPrInc,nDelta,lPrInc, &
                   Coeff1,iBasi,Coeff2,jBasj,Coeff3,kBask,Coeff4,lBasl, &
                   nZeta,nEta,Grad,nGrad,IfGrad,IndGrd,PSO,nPSO, &
                   Wrk2,nWrk2,Aux,nAux,Shijij)
!***********************************************************************
!                                                                      *
! Object: to generate the SO integrals for four fixed centers and      *
!         fixed basis set types.                                       *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!          Roland Lindh, Dept. of Theoretical Chemistry, University of *
!          Lund, SWEDEN. Modified to gradients, January '92.           *
!***********************************************************************

use Real_Spherical, only: ipSph, RSph
use Basis_Info, only: MolWgh, Shells
use Center_Info, only: dc
use Phase_Info, only: iPhase
use Gateway_Info, only: ChiI2
use Gateway_global, only: IsChi
use Symmetry_Info, only: nIrrep
use Index_Functions, only: nTri_Elem1
use k2_structure, only: k2_type
use k2_arrays, only: BraKet
use Disp, only: CutGrd, l2DI
#ifdef _DEBUGPRINT_
use Disp, only: ChDisp
#endif
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iAnga(4), iCmp(4), iShell(4), iShll(4), iAO(4), iStb, jStb, kStb, lStb, nRys, nData1, nData2, &
                                 nAlpha, iPrInc, nBeta, jPrInc, nGamma, kPrInc, nDelta, lPrInc, iBasi, jBasj, kBask, lBasl, nZeta, &
                                 nEta, nGrad, IndGrd(3,4), nPSO, nWrk2, nAux
real(kind=wp), intent(in) :: Coor(3,4), Coeff1(nAlpha,iBasi), Coeff2(nBeta,jBasj), Coeff3(nGamma,kBask), Coeff4(nDelta,lBasl), &
                             PSO(iBasi*jBasj*kBask*lBasl,nPSO)
type(k2_type), intent(in) :: k2Data1(nData1), k2Data2(nData2)
real(kind=wp), intent(inout) :: Pren, Prem, Grad(nGrad)
logical(kind=iwp), intent(in) :: IfGrad(3,4), Shijij
real(kind=wp), intent(out) :: Wrk2(nWrk2), Aux(nAux)
integer(kind=iwp) :: iC, iCar, iCent, iCmpa, iDCRR(0:7), iDCRS(0:7), iDCRT(0:7), iDCRTS, iEta, iiCent, ijklab, ijMax, ijMin, ikl, &
                     IncEta, IncZet, iShlla, iStabM(0:7), iStabN(0:7), iuvwx(4), iW2, iW3, iW4, ix1, ix2, ixSh, iy1, iy2, iz1, &
                     iz2, iZeta, jCent, jCmpb, jjCent, JndGrd(3,4), jShllb, kCent, kCmpc, klMax, klMin, kOp(4), kShllc, la, lb, &
                     lc, lCent, lCmpd, ld, lDCR1, lDCR2, lDCRR, lDCRS, lDCRT, lEta, LmbdR, LmbdS, LmbdT, lShlld, lStabM, lStabN, &
                     lZeta, mab, mcd, mCent, mEta, mGrad, MxDCRS, mZeta, nDCRR, nDCRS, nDCRT, nEta_Tot, nIdent, nijkl, nOp(4), &
                     nW2, nW4, nWrk3, nZeta_Tot
real(kind=wp) :: Aha, CoorAC(3,2), CoorM(3,4), Fact, u, v, w, x
logical(kind=iwp) :: ABeqCD, AeqB, AeqC, CeqD, JfGrad(3,4), PreScr
integer(kind=iwp), external :: NrOpr
real(kind=wp), external :: DDot_
logical(kind=iwp), external :: EQ, lEmpty
external :: ModU2, TERI1, vCff2D
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i, iPrint, iRout
character(len=3), parameter :: ChOper(0:7) = [' E ',' x ',' y ',' xy',' z ',' xz',' yz','xyz']
#include "print.fh"
#endif

#include "macros.fh"
unused_var(iPrInc)
unused_var(kPrInc)

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
iRout = 12
iPrint = nPrint(iRout)
#endif
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
IncZet = nAlpha*jPrInc
IncEta = nGamma*lPrInc
LmbdT = 0
nijkl = iBasi*jBasj*kBask*lBasl
iuvwx(1) = dc(iStb)%nStab
iuvwx(2) = dc(jStb)%nStab
iuvwx(3) = dc(kStb)%nStab
iuvwx(4) = dc(lStb)%nStab
mab = nTri_Elem1(la)*nTri_Elem1(lb)
mcd = nTri_Elem1(lc)*nTri_Elem1(ld)
iW4 = 1
if ((jPrInc /= nBeta) .or. (lPrInc /= nDelta)) then
  iW2 = 1+mab*mcd*nijkl
else
  iW2 = 1
end if

!                                                                      *
!***********************************************************************
!                                                                      *
! Find the Double Coset Representatives for center A and B

if (nIrrep == 1) then
  nDCRR = 1
  iDCRR(0) = 0
  LmbdR = 1
else
  call DCR(LmbdR,dc(iStb)%iStab,dc(iStb)%nStab,dc(jStb)%iStab,dc(jStb)%nStab,iDCRR,nDCRR)
end if
#ifdef _DEBUGPRINT_
if (iPrint >= 99) write(u6,'(20A)') ' {R}=(',(ChOper(iDCRR(i)),',',i=0,nDCRR-1),')'
#endif
u = real(dc(iStb)%nStab,kind=wp)
v = real(dc(jStb)%nStab,kind=wp)

! Find stabilizer for center A and B

if (nIrrep == 1) then
  lStabM = 1
  iStabM(0) = 0
else
  call Inter(dc(iStb)%iStab,dc(iStb)%nStab,dc(jStb)%iStab,dc(jStb)%nStab,iStabM,lStabM)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Find the Double Coset Representatives for center C and D.
! Take care of redundancy if {f(aA)f(bB)}={f(cC)f(dD)}. Hence
! we will only use unique combinations of operators from the
! double coset representatives {R} and {S}.

if (nIrrep == 1) then
  nDCRS = 1
  iDCRS(0) = 0
  LmbdS = 1
else
  call DCR(LmbdS,dc(kStb)%iStab,dc(kStb)%nStab,dc(lStb)%iStab,dc(lStb)%nStab,iDCRS,nDCRS)
end if
#ifdef _DEBUGPRINT_
if (iPrint >= 99) write(u6,'(20A)') ' {S}=(',(ChOper(iDCRS(i)),',',i=0,nDCRS-1),')'
#endif
w = real(dc(kStb)%nStab,kind=wp)
x = real(dc(lStb)%nStab,kind=wp)

! Find stabilizer for center C and D

if (nIrrep == 1) then
  lStabN = 1
  iStabN(0) = 0
else
  call Inter(dc(kStb)%iStab,dc(kStb)%nStab,dc(lStb)%iStab,dc(lStb)%nStab,iStabN,lStabN)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Find the Double Coset Representatives for the two charge
! distributions.

if (nIrrep == 1) then
  nDCRT = 1
  iDCRT(0) = 0
  LmbdT = 1
else
  call DCR(LmbdT,iStabM,lStabM,iStabN,lStabN,iDCRT,nDCRT)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Factor due to summation over DCR

if (MolWgh == 1) then
  Fact = real(nIrrep,kind=wp)/real(LmbdT,kind=wp)
else if (MolWgh == 0) then
  Fact = u*v*w*x/real(nIrrep**3*LmbdT,kind=wp)
else
  Fact = sqrt(u*v*w*x)/real(nIrrep*LmbdT,kind=wp)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
nOp(1) = NrOpr(0)
CoorM(:,1) = Coor(:,1)
do lDCRR=0,nDCRR-1
  nOp(2) = NrOpr(iDCRR(lDCRR))
  call OA(iDCRR(lDCRR),Coor(:,2),CoorM(:,2))
  AeqB = EQ(CoorM(:,1),CoorM(:,2))

  MxDCRS = nDCRS-1
  do lDCRS=0,MxDCRS
    CoorM(:,3) = Coor(:,3)
    call OA(iDCRS(lDCRS),Coor(1:3,4),CoorM(1:3,4))
    CeqD = EQ(Coor(:,3),CoorM(:,4))

    do lDCRT=nDCRT-1,0,-1
#     ifdef _DEBUGPRINT_
      if (iPrint >= 99) write(u6,'(6A)') ' R=',ChOper(iDCRR(lDCRR)),', S=',ChOper(iDCRS(lDCRS)),', T=',ChOper(iDCRT(lDCRT))
#     endif

      nOp(3) = NrOpr(iDCRT(lDCRT))
      iDCRTS = ieor(iDCRT(lDCRT),iDCRS(lDCRS))
      nOp(4) = NrOpr(iDCRTS)

      call OA(iDCRTS,Coor(:,4),CoorM(:,4))
      call OA(iDCRT(lDCRT),Coor(:,3),CoorM(:,3))

#     ifdef _DEBUGPRINT_
      if (iPrint >= 59) call RecPrt(' CoorM in TwoEl',' ',CoorM,3,4)
#     endif
      AeqC = EQ(CoorM(:,1),CoorM(:,3))
      ABeqCD = AeqB .and. CeqD .and. AeqC
      ! No contribution to gradient from one-center integrals
      if (ABeqCD) cycle

      ! Modify which center we will differetiate with
      ! respect to using the translational invariance.

      mGrad = 0
      do iCar=1,3
        ! Copy to temporary arrays
        JfGrad(iCar,:) = IfGrad(iCar,:)
        JndGrd(iCar,:) = IndGrd(iCar,:)
        ! In case of four differentiations use
        ! the translational invariance to remove
        ! one or several of them.
        if (IfGrad(iCar,1) .and. IfGrad(iCar,2) .and. IfGrad(iCar,3) .and. IfGrad(iCar,4)) then

          nIdent = 0
          kCent = 0
          lCent = 0
          do iCent=1,3
            do jCent=iCent+1,4
              if (EQ(CoorM(1,iCent),CoorM(1,jCent))) then
                nIdent = nIdent+1
                kCent = iCent
                lCent = jCent
              end if
            end do
          end do

          ! Remove a center which is unique and if possible
          ! if it occurs several times.

          if (nIdent == 0) then
            ! Four center case, remove a center
            jCent = 1
            do iCent=1,4
              if (iAnga(iCent) /= 0) jCent = iCent
            end do
            JndGrd(iCar,jCent) = -JndGrd(iCar,jCent)
            JfGrad(iCar,jCent) = .false.
          else if (nIdent == 1) then
            ! Three center case
            iiCent = kCent
            jjCent = lCent
            JndGrd(iCar,iiCent) = 0
            JfGrad(iCar,iiCent) = .false.
            JndGrd(iCar,jjCent) = -JndGrd(iCar,jjCent)
            JfGrad(iCar,jjCent) = .false.
          else if (nIdent == 2) then
            ! Two center case
            iiCent = kCent
            jjCent = lCent
            ikl = ibset(0,kCent-1)
            ikl = ibset(ikl,lCent-1)
            do iC=4,1,-1
              if (.not. btest(ikl,iC-1)) then
                iCent = iC
                exit
              end if
            end do
            ikl = ibset(ikl,iCent-1)
            do iC=4,1,-1
              if (.not. btest(ikl,iC-1)) then
                jCent = iC
                exit
              end if
            end do
            ijMax = max(iAnga(iCent),iAnga(jCent))
            ijMin = min(iAnga(iCent),iAnga(jCent))
            klMax = max(iAnga(kCent),iAnga(lCent))
            klMin = min(iAnga(kCent),iAnga(lCent))
            if (klMin > 0) then
              iiCent = kCent
              jjCent = lCent
            else if (ijMin > 0) then
              iiCent = iCent
              jjCent = jCent
            else if (klMax > 0) then
              iiCent = kCent
              jjCent = lCent
            else if (ijMax > 0) then
              iiCent = iCent
              jjCent = jCent
            end if
            JndGrd(iCar,iiCent) = 0
            JfGrad(iCar,iiCent) = .false.
            JndGrd(iCar,jjCent) = -JndGrd(iCar,jjCent)
            JfGrad(iCar,jjCent) = .false.
          else if (nIdent == 3) then
            ! Two center case
            if (lCent /= 4) then
              mCent = 1
            else if (kCent /= 3) then
              mCent = 1
            else if (EQ(CoorM(1,1),CoorM(1,4))) then
              mCent = 1
            else
              mCent = 2
            end if
            JndGrd(iCar,mCent) = 0
            JfGrad(iCar,mCent) = .false.
            JndGrd(iCar,kCent) = 0
            JfGrad(iCar,kCent) = .false.
            JndGrd(iCar,lCent) = -JndGrd(iCar,lCent)
            JfGrad(iCar,lCent) = .false.
          else
            call WarningMessage(2,'Error in Twoel_g')
            write(u6,*) ' Twoel: nIdent too large!'
            call Abend()
          end if
        end if
        do ixSh=1,4
          if (JfGrad(iCar,ixSh)) mGrad = mGrad+1
        end do
      end do
      if (mGrad == 0) cycle

      ! Find the proper centers to start of with the angular
      ! momentum on. If la == lb there will exist an
      ! ambiguity to which center that angular momentum should
      ! be accumulated on. In that case we will use A and C of
      ! the order as defined by the basis functions types.

      if (iAnga(1) >= iAnga(2)) then
        CoorAC(:,1) = CoorM(:,1)
      else
        CoorAC(:,1) = CoorM(:,2)
      end if
      if (iAnga(3) >= iAnga(4)) then
        CoorAC(:,2) = CoorM(:,3)
      else
        CoorAC(:,2) = CoorM(:,4)
      end if
      kOp(:) = nOp(:)

      ! Desymmetrize the second order density matrix

      ! (faA fbR(B) | fcT(C) fdTS(D))ijkl

      call DesymP(iAnga,iCmp(1),iCmp(2),iCmp(3),iCmp(4),Shijij,iShll,iShell,iAO,kOp,nijkl,Aux,nAux,Wrk2(iW2),PSO,nPSO)

      if (Fact /= One) call DScal_(nijkl*iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4),Fact,Wrk2(iW2),1)

      ! Backtransform 2nd order density matrix from spherical
      ! harmonic gaussians to cartesian gaussians.

      ijklab = nijkl*iCmp(1)*iCmp(2)
      nW2 = ijklab*max(kCmpc*lCmpd,mcd)
      iW3 = iW2+nW2
      nWrk3 = nWrk2-((iW2-iW4)+nW2)
      call SphCr1(Wrk2(iW2),ijklab,Wrk2(iW3),nWrk3,RSph(ipSph(lc)),nTri_Elem1(lc),kCmpc,Shells(kShllc)%Transf, &
                  Shells(kShllc)%Prjct,RSph(ipSph(ld)),nTri_Elem1(ld),lCmpd,Shells(lShlld)%Transf,Shells(lShlld)%Prjct, &
                  Wrk2(iW2),mcd)
      if (iW2 == iW4) then
        nW2 = nijkl*mcd*max(iCmpa*jCmpb,mab)
        nW4 = 0
      else
        nW2 = nijkl*mcd*iCmpa*jCmpb
        nW4 = nijkl*mcd*mab
      end if
      iW3 = iW2+nW2
      nWrk3 = nWrk2-(nW2+nW4)
      call SphCr2(Wrk2(iW2),nijkl,mcd,Wrk2(iW3),nWrk3,RSph(ipSph(la)),nTri_Elem1(la),iCmpa,Shells(iShlla)%Transf, &
                  Shells(iShlla)%Prjct,RSph(ipSph(lb)),nTri_Elem1(lb),jCmpb,Shells(jShllb)%Transf,Shells(jShllb)%Prjct, &
                  Wrk2(iW4),mab)

      ! Transpose the 2nd order density matrix

      if (mab*mcd /= 1) then
        iW3 = iW4+nijkl*mab*mcd
        call DGetMO(Wrk2(iW4),nijkl,nijkl,mab*mcd,Wrk2(iW3),mab*mcd)
        call dcopy_(mab*mcd*nijkl,Wrk2(iW3),1,Wrk2(iW4),1)
      end if

      lDCR1 = NrOpr(iDCRR(lDCRR))+1
      lDCR2 = NrOpr(iDCRS(lDCRS))+1
      ix1 = 1
      iy1 = 1
      iz1 = 1
      ix2 = iPhase(1,iDCRT(lDCRT))
      iy2 = iPhase(2,iDCRT(lDCRT))
      iz2 = iPhase(3,iDCRT(lDCRT))

      nZeta_Tot = k2Data1(lDCR1)%IndZ(nZeta+1)
      nEta_Tot = k2Data2(lDCR2)%IndZ(nEta+1)

      ! Loops to partion the primitives

      do iZeta=1,nZeta_Tot,IncZet
        mZeta = min(IncZet,nZeta_Tot-iZeta+1)
        if (lEmpty(Coeff2,nBeta,nBeta,jBasj)) cycle

        do iEta=1,nEta_Tot,IncEta
          mEta = min(IncEta,nEta_Tot-iEta+1)
          if (lEmpty(Coeff4,nDelta,nDelta,lBasl)) cycle

          Pren = Pren+real(mab*mcd*mZeta*mEta,kind=wp)

          ! Preprescreen

          call PrePre_g(mZeta,mEta,lZeta,lEta,PreScr,CutGrd,iZeta-1,iEta-1,k2Data1(lDCR1),k2Data2(lDCR2))
          if (lZeta*lEta == 0) cycle

          ! Decontract the 2nd order density matrix

          if (iW4 == iW2) then
            nW4 = 0
            nW2 = max(nijkl,mZeta*mEta)*mab*mcd
          else
            nW4 = nijkl*mab*mcd
            nW2 = mZeta*mEta*mab*mcd
          end if
          iW3 = iW2+nW2
          nWrk3 = nWrk2-(nW4+nW2)
          call Tcrtnc(Coeff1,nAlpha,iBasi,Coeff2,nBeta,jBasj,Coeff3,nGamma,kBask,Coeff4,nDelta,lBasl, &
                      Wrk2(iW4),mab*mcd,Wrk2(iW3),nWrk3,Wrk2(iW2), &
                      k2Data1(lDCR1)%IndZ(iZeta:iZeta+mZeta-1),mZeta, &
                      k2Data2(lDCR2)%IndZ(iEta:iEta+mEta-1),mEta)

          ! Transfer k2 data and prescreen

          iW3 = iW2+mZeta*mEta*mab*mcd
          nWrk3 = nWrk2-mZeta*mEta*mab*mcd
          call Screen_g(iZeta-1,iEta-1,Wrk2(iW2),Wrk2(iW3),mab*mcd,nZeta,nEta,mZeta,mEta,lZeta,lEta, &
                        k2Data1(lDCR1),k2Data2(lDCR2), &
                        Braket%Zeta,Braket%ZInv,Braket%P,Braket%xA,Braket%xB, &
                        Braket%Eta,Braket%EInv,Braket%Q,Braket%xG,Braket%xD, &
                        ix1,iy1,iz1,ix2,iy2,iz2,CutGrd,l2DI, &
                        PreScr,nWrk3,IsChi,ChiI2)
          Prem = Prem+real(mab*mcd*lZeta*lEta,kind=wp)
          !write(u6,*) 'Prem=',Prem
          if (lZeta*lEta == 0) cycle

          ! Compute integral derivative and accumulate
          ! contribution to the molecular gradient. Note that
          ! the PSO matrix now is stored in Wrk2(iW2).

          iW3 = iW2+lZeta*lEta*mab*mcd
          call Rysg1(iAnga,nRys,lZeta*lEta,BraKet%xA,BraKet%xB,BraKet%xG,BraKet%xD, &
                     BraKet%Zeta,BraKet%ZInv,lZeta,BraKet%Eta,BraKet%EInv,lEta, &
                     BraKet%P,nZeta,BraKet%Q,nEta,CoorM,CoorM,CoorAC, &
                     Wrk2(iW3),nWrk3,TERI1,ModU2,vCff2D,Wrk2(iW2),mab*mcd,Grad,nGrad,JfGrad,JndGrd,kOp,iuvwx)
          Aha = sqrt(DDot_(nGrad,Grad,1,Grad,1))
          if (Aha > 1.0e5_wp) then
            write(u6,*) 'Norm of gradient contribution is huge!'
            write(u6,*) 'Probably due to wrong coordinates.'
          end if

        end do
      end do

#     ifdef _DEBUGPRINT_
      if (iPrint >= 19) call PrGrad(' In TwoEl',Grad,nGrad,ChDisp)
#     endif

    end do
  end do
end do

return

end subroutine TwoEl_g
