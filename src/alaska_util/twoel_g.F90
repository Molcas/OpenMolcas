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

subroutine TwoEl_g(Coor,nRys,Grad,nGrad,IfGrad,IndGrd,PSO,nijkl,nPSO,Wrk2,nWrk2,iSD4)
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

use setup, only: nAux
use Real_Spherical, only: ipSph, RSph
use Basis_Info, only: Shells
use Phase_Info, only: iPhase
use Gateway_Info, only: ChiI2
use iSD_data, only: nSD
use Gateway_global, only: IsChi
use Index_Functions, only: iTri, nTri_Elem1
use k2_structure, only: Indk2, k2_type, k2Data
use k2_arrays, only: Aux, BraKet
use Disp, only: CutGrd, l2DI
use Rys_interfaces, only: cff2d_kernel, modu2_kernel, tval1_kernel
#ifdef _DEBUGPRINT_
use Symmetry_Info, only: ChOper
#endif
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nRys, nGrad, IndGrd(3,4), nijkl, nPSO, nWrk2, iSD4(0:nSD,4)
real(kind=wp), intent(in) :: Coor(3,4), PSO(nijkl,nPSO)
real(kind=wp), intent(inout) :: Grad(nGrad)
logical(kind=iwp), intent(in) :: IfGrad(3,4)
real(kind=wp), intent(out) :: Wrk2(nWrk2)
integer(kind=iwp) :: iAnga(4), iAO(4), iAOst(4), iBasi, iC, iCar, iCent, iCmp(4), iCmpa, iDCRR(0:7), iDCRS(0:7), iDCRT(0:7), &
                     iDCRTS, iEta, iiCent, ijklab, ijMax, ijMin, ijS, ik2, ikl, IncEta, IncZet, iS, iShell(4), iShll(4), iShlla, &
                     iuvwx(4), iW2, iW3, iW4, ix1, ix2, ixSh, iy1, iy2, iz1, iz2, iZeta, jBasj, jCent, jCmpb, jjCent, jk2, &
                     JndGrd(3,4), jPrInc, jS, jShllb, kBask, kCent, kCmpc, klMax, klMin, klS, kOp(4), kS, kShllc, la, lb, lBasl, &
                     lc, lCent, lCmpd, ld, lDCR1, lDCR2, lDCRR, lDCRS, lDCRT, lEta, lPrInc, lS, lShlld, lZeta, mab, mcd, mCent, &
                     mEta, mGrad, MxDCRS, mZeta, nAlpha, nBeta, nDCR1, nDCR2, nDCRR, nDCRS, nDCRT, nDelta, nEta, nEta_Tot, nGamma, &
                     nIdent, nOp(4), nW2, nW4, nWrk3, nZeta, nZeta_Tot
real(kind=wp) :: Aha, CoorAC(3,2), CoorM(3,4), Fact
logical(kind=iwp) :: ABeqCD, AeqB, AeqC, CeqD, JfGrad(3,4), PreScr, Shijij
procedure(cff2d_kernel) :: vCff2D
procedure(modu2_kernel) :: ModU2
procedure(tval1_kernel) :: TERI1
real(kind=wp), pointer :: Coeff1(:,:), Coeff2(:,:), Coeff3(:,:), Coeff4(:,:)
type(k2_type), pointer :: k2data1(:), k2data2(:)
integer(kind=iwp), external :: NrOpr
real(kind=wp), external :: DDot_
logical(kind=iwp), external :: EQ

!                                                                      *
!***********************************************************************
!                                                                      *

jPrInc = iSD4(6,2)
lPrInc = iSD4(6,4)

iAnga(:) = iSD4(1,:)
iCmp(:) = iSD4(2,:)
iShll(:) = iSD4(0,:)
iShell(:) = iSD4(11,:)
iAO(:) = iSD4(7,:)
iAOst(:) = iSD4(8,:)

Shijij = (iSD4(11,1) == iSD4(11,3)) .and. (iSD4(11,2) == iSD4(11,4))

nAlpha = iSD4(5,1)
nBeta = iSD4(5,2)
nGamma = iSD4(5,3)
nDelta = iSD4(5,4)

iBasi = iSD4(19,1)
jBasj = iSD4(19,2)
kBask = iSD4(19,3)
lBasl = iSD4(19,4)

nZeta = nAlpha*nBeta
nEta = nGamma*nDelta
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
! Pick up pointers to k2 entities.

iS = iShell(1)
jS = iShell(2)
kS = iShell(3)
lS = iShell(4)
ijS = iTri(iS,jS)
klS = iTri(kS,lS)
nDCR1 = IndK2(2,ijS)
ik2 = IndK2(3,ijS)
nDCR2 = IndK2(2,klS)
jk2 = IndK2(3,klS)
k2data1(1:nDCR1) => k2Data(1:nDCR1,ik2)
k2data2(1:nDCR2) => k2Data(1:nDCR2,jk2)

Coeff1(1:nAlpha,1:iBasi) => Shells(iShll(1))%pCff(:,iAOst(1)+1:)
Coeff2(1:nBeta,1:jBasj) => Shells(iShll(2))%pCff(:,iAOst(2)+1:)
Coeff3(1:nGamma,1:kBask) => Shells(iShll(3))%pCff(:,iAOst(3)+1:)
Coeff4(1:nDelta,1:lBasl) => Shells(iShll(4))%pCff(:,iAOst(4)+1:)

call mk_DCRs_and_Stabilizers(Fact,iuvwx,nDCRR,nDCRS,nDCRT,iDCRR,iDCRS,iDCRT,nSD,iSD4)
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
      write(u6,'(6A)') ' R=',ChOper(iDCRR(lDCRR)),', S=',ChOper(iDCRS(lDCRS)),', T=',ChOper(iDCRT(lDCRT))
#     endif

      nOp(3) = NrOpr(iDCRT(lDCRT))
      iDCRTS = ieor(iDCRT(lDCRT),iDCRS(lDCRS))
      nOp(4) = NrOpr(iDCRTS)

      call OA(iDCRTS,Coor(:,4),CoorM(:,4))
      call OA(iDCRT(lDCRT),Coor(:,3),CoorM(:,3))

#     ifdef _DEBUGPRINT_
      call RecPrt(' CoorM in TwoEl',' ',CoorM,3,4)
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

      ! Loops to partition the primitives

      do iZeta=1,nZeta_Tot,IncZet
        mZeta = min(IncZet,nZeta_Tot-iZeta+1)
        if (all(Coeff2(:,:) == Zero)) cycle

        do iEta=1,nEta_Tot,IncEta
          mEta = min(IncEta,nEta_Tot-iEta+1)
          if (all(Coeff4(:,:) == Zero)) cycle

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
      call PrGrad(' In TwoEl',Grad,nGrad)
#     endif

    end do
  end do
end do

end subroutine TwoEl_g
