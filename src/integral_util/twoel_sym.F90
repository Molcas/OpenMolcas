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
! Copyright (C) 1990,1991,1993, Roland Lindh                           *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine TwoEl_Sym(iS_,jS_,kS_,lS_,Coor,iAnga,iCmp,iShell,iShll,iAO,iAOst,NoInts,iStabs,nAlpha,iPrInc,nBeta,jPrInc,nGamma, &
                     kPrInc,nDelta,lPrInc,nData1,nData2,k2Data1,k2Data2,IJeqKL,kOp,Dij,mDij,mDCRij,Dkl,mDkl,mDCRkl,Dik,mDik, &
                     mDCRik,Dil,mDil,mDCRil,Djk,mDjk,mDCRjk,Djl,mDjl,mDCRjl,Coeff1,iBasi,Coeff2,jBasj,Coeff3,kBask,Coeff4,lBasl, &
                     FckTmp,nFT,nZeta,nEta,SOInt,nSOInt,Wrk,nWork2,Shijij,Aux,nAux)
!***********************************************************************
!                                                                      *
! Object: to generate the SO integrals for four fixed centers and      *
!         fixed basis set types.                                       *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!          Roland Lindh, Dept. of Theoretical Chemistry, University of *
!          Lund, SWEDEN. Modified to use Schwartz inequality for pre-  *
!          screening, July 1991.                                       *
!          Modified for direct SCF, January '93                        *
!***********************************************************************
use Real_Spherical
use Basis_Info
use Center_Info
use Phase_Info
use Gateway_Info, only: ThrInt, CutInt
use Symmetry_Info, only: nIrrep
use Int_Options, only: DoIntegrals, DoFock, FckNoClmb, FckNoExch
use Int_Options, only: ExFac, Thize, W2Disc, IntOnly => PreSch
use Int_Options, only: Disc_Mx, Disc, Quad_ijkl
use k2_arrays, only: TwoHam => pFq, Dens => pDq
use Breit, only: nOrdOp, nComp
use Constants
#ifdef _DEBUGPRINT_
use Symmetry_Info, only: ChOper
#endif
use k2_structure, only: k2_type

implicit none
#include "twoswi.fh"
integer iS_, jS_, kS_, lS_, nAlpha, iPrInc, nBeta, jPrInc, nGamma, kPrInc, nDelta, lPrInc, nData1, nData2, mDij, mDCRij, mDkl, &
        mDCRkl, mDik, mDCRik, mDil, mDCRil, mDjk, mDCRjk, mDjl, mDCRjl, iBasi, jBasj, kBask, lBasl, nFT, nZeta, nEta, nSOInt, &
        nWork2, nAux
real*8 Coor(3,4)
integer iAnga(4), iCmp(4), iShell(4), iShll(4), iAO(4), iAOst(4)
logical NoInts
integer iStabs(4)
type(k2_type) k2data1(nData1), k2Data2(nData2)
logical IJeqKL
integer kOp(4)
real*8 Dij(mDij,mDCRij), Dkl(mDkl,mDCRkl), Dik(mDik,mDCRik), Dil(mDil,mDCRil), Djk(mDjk,mDCRjk), Djl(mDjl,mDCRjl), &
       Coeff1(nAlpha,iBasi), Coeff2(nBeta,jBasj), Coeff3(nGamma,kBask), Coeff4(nDelta,lBasl), FckTmp(nFT)
real*8 SOInt(iBasi*jBasj*kBask*lBasl,nSOInt), Wrk(nWork2)
logical Shijij
real*8 Aux(nAux)
real*8 CoorM(3,4), CoorAC(3,2), QInd(2)
integer iDCRR(0:7), iDCRS(0:7), iDCRT(0:7), iStabN(0:7), iStabM(0:7), iWR(2)
logical NoPInts, AeqB, CeqD, AeqC, ABeqCD, Do_TnsCtl, IeqK, JeqL, Prescreen_On_Int_Only, DoCoul, DoExch, Scrij, Scrkl, Scrik, &
        Scril, Scrjk, Scrjl, Batch_On_Disk, DoAOBatch, All_Spherical
integer iStb, jStb, kStb, lStb
logical :: Copy = .true., NoCopy = .false.
integer :: jOp(6) = [0,0,0,0,0,0]
#include "SysDef.fh"
logical, external :: EQ, lEmpty
integer la, lb, lc, ld, ISMAng, nab, ncd, nijkl, nInts, ipAOInt, iW3, iW4, kInts, mInts, mabMin, mabMax, mcdMin, mcdMax, mabcd, &
        IncZet, IncEta, mWork2, nZeta_Tot, nEta_Tot, kabcd, ipAOInt_, mZeta, mEta, iOpt, i_Int, nByte, nabcd, iW4_
real*8 RST_Triplet, vijkl, vij, vkl, vik, vil, vjk, vjl
integer ix1, iy1, iz1, ix2, iy2, iz2, iR, iS, iT, iTS, iRT, iRTS, ij1, ij2, ij3, ij4, kl1, kl2, kl3, kl4, ik1, ik2, ik3, ik4, il1, &
        il2, il3, il4, jk1, jk2, jk3, jk4, jl1, jl2, jl3, jl4, nDCRR, lDCR1, MxDCRS, nDCRS, lDCR2, nDCRT, lDCRE_, lDCRT_, mWork3, &
        LmbdR, LmbdS, LmbdT, iDCRTS, lStabM, lStabN, iZeta, iEta, lDCRR, lDCRS, lDCRT
integer, external :: NrOpr
real*8 u, v, w, x, vijij, RS_doublet, FactNd
#ifdef _DEBUGPRINT_
integer i
#endif
interface
  subroutine FckAcc(iAng,iCmp_,Shijij,iShll,iShell,kOp,nijkl,AOInt,TwoHam,nDens,Scrt,nScrt,iAO,iAOst,iBas,jBas,kBas,lBas,Dij,ij1, &
                    ij2,ij3,ij4,Dkl,kl1,kl2,kl3,kl4,Dik,ik1,ik2,ik3,ik4,Dil,il1,il2,il3,il4,Djk,jk1,jk2,jk3,jk4,Djl,jl1,jl2,jl3, &
                    jl4,FT,nFT,DoCoul,DoExch,ExFac)
    integer nijkl, nDens, nScrt, nFT, ij1, ij2, ij3, ij4, kl1, kl2, kl3, kl4, ik1, ik2, ik3, ik4, il1, il2, il3, il4, jk1, jk2, &
            jk3, jk4, jl1, jl2, jl3, jl4
    integer iBas, jBas, kBas, lBas
    real*8, target :: Scrt(nScrt)
    integer iCmp_(4)
    real*8 AOInt(nijkl,iCmp_(1),iCmp_(2),iCmp_(3),iCmp_(4)), TwoHam(nDens)
    real*8, target :: Dij(ij1*ij2+1,ij3,ij4), Dkl(kl1*kl2+1,kl3,kl4), Dik(ik1*ik2+1,ik3,ik4), Dil(il1*il2+1,il3,il4), &
                      Djk(jk1*jk2+1,jk3,jk4), Djl(jl1*jl2+1,jl3,jl4)
    real*8, target :: FT(nFT)
    logical Shijij, DoCoul, DoExch
    integer iAng(4), iShell(4), iShll(4), kOp(4), iAO(4), iAOst(4)
    real*8 ExFac
  end subroutine FckAcc
end interface
! Statement function
integer :: ixyz, nabSz
nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6-1

iStb = iStabs(1)
jStb = iStabs(2)
kStb = iStabs(3)
lStb = iStabs(4)
if (nOrdOp /= 0) then
  write(6,*) 'Breit two-electron integrals not implemented yet'
  write(6,*) 'Symmetry adaptation different since the operator'
  write(6,*) 'is not symmetric.'
end if

All_Spherical = Shells(iShll(1))%Prjct .and. Shells(iShll(2))%Prjct .and. Shells(iShll(3))%Prjct .and. Shells(iShll(4))%Prjct
QInd(1) = Quad_ijkl
RST_triplet = One

Do_tnsctl = .false.
kabcd = 0

#ifdef _DEBUGPRINT_
call RecPrt('Coeff1',' ',Coeff1,nAlpha,iBasi)
call RecPrt('Coeff2',' ',Coeff2,nBeta,jBasj)
call RecPrt('Coeff3',' ',Coeff3,nGamma,kBask)
call RecPrt('Coeff4',' ',Coeff4,nDelta,lBasl)
call RecPrt('Coor',' ',Coor,3,4)
#endif

la = iAnga(1)
lb = iAnga(2)
lc = iAnga(3)
ld = iAnga(4)
iSmAng = la+lb+lc+ld
LmbdT = 0
nijkl = iBasi*jBasj*kBask*lBasl*nComp
nab = iCmp(1)*iCmp(2)
ncd = iCmp(3)*iCmp(4)
nabcd = nab*ncd
nInts = nijkl*nabcd
SOInt(:,:) = Zero
ipAOInt = 1
iW3 = 1+nInts
iW4 = 1
!                                                                      *
!***********************************************************************
!                                                                      *
! Find the Double Coset Representatives for center A and B

call DCR(LmbdR,dc(iStb)%iStab,dc(iStb)%nStab,dc(jStb)%iStab,dc(jStb)%nStab,iDCRR,nDCRR)
u = dble(dc(iStb)%nStab)
v = dble(dc(jStb)%nStab)
#ifdef _DEBUGPRINT_
write(6,'(20A)') ' {R}=(',(ChOper(iDCRR(i)),',',i=0,nDCRR-1),')'
#endif

! Find stabilizer for center A and B

call Inter(dc(iStb)%iStab,dc(iStb)%nStab,dc(jStb)%iStab,dc(jStb)%nStab,iStabM,lStabM)

!                                                                      *
!***********************************************************************
!                                                                      *
! Find the Double Coset Representatives for center C and D.

call DCR(LmbdS,dc(kStb)%iStab,dc(kStb)%nStab,dc(lStb)%iStab,dc(lStb)%nStab,iDCRS,nDCRS)
w = dble(dc(kStb)%nStab)
x = dble(dc(lStb)%nStab)
#ifdef _DEBUGPRINT_
write(6,'(20A)') ' {S}=(',(ChOper(iDCRS(i)),',',i=0,nDCRS-1),')'
#endif

! Find stabilizer for center C and D

call Inter(dc(kStb)%iStab,dc(kStb)%nStab,dc(lStb)%iStab,dc(lStb)%nStab,iStabN,lStabN)
!                                                                      *
!***********************************************************************
!                                                                      *
! Find the Double Coset Representatives for the two charge
! distributions.

call DCR(LmbdT,iStabM,lStabM,iStabN,lStabN,iDCRT,nDCRT)
!                                                                      *
!***********************************************************************
!                                                                      *
kOp(1) = NrOpr(0)
call dcopy_(3,Coor(1,1),1,CoorM(1,1),1)
do lDCRR=0,nDCRR-1
  kOp(2) = NrOpr(iDCRR(lDCRR))
  call OA(iDCRR(lDCRR),Coor(1:3,2),CoorM(1:3,2))
  AeqB = EQ(CoorM(1,1),CoorM(1,2))

  lDCR1 = NrOpr(iDCRR(lDCRR))+1

  vijij = k2Data1(lDCR1)%abMax

  ! switch (to generate better start orbitals...)
  if (NDDO .and. (.not. AeqB)) Go To 100
  ! switch
  MxDCRS = nDCRS-1
  do lDCRS=0,MxDCRS
    RS_doublet = dble(lDCRS*nDCRR+lDCRR+1)
    call dcopy_(3,Coor(1,3),1,CoorM(1,3),1)
    call OA(iDCRS(lDCRS),Coor(1:3,4),CoorM(1:3,4))
    CeqD = EQ(Coor(1,3),CoorM(1,4))

    ! switch (to generate better start orbitals...)
    if (NDDO .and. (.not. CeqD)) Go To 200

    lDCR2 = NrOpr(iDCRS(lDCRS))+1

    ! Pick up estimated largest integral value (AO)

    vijkl = vijij*k2Data2(lDCR2)%abMax
    do lDCRT=0,nDCRT-1
      ipAOInt = 1
      iW3 = 1+nInts
      RS_doublet = dble(lDCRS*nDCRR+lDCRR+1)
      RST_triplet = dble(lDCRT*nDCRR*nDCRS)+RS_doublet
      QInd(2) = RST_triplet
      !write(6,*) QInd(1), QInd(2)
      iDCRTS = ieor(iDCRT(lDCRT),iDCRS(lDCRS))
      call OA(iDCRTS,Coor(1:3,4),CoorM(1:3,4))
      call OA(iDCRT(lDCRT),Coor(1:3,3),CoorM(1:3,3))
      AeqC = EQ(CoorM(1,1),CoorM(1,3))
      ABeqCD = AeqB .and. CeqD .and. AeqC
      if (ABeqCD .and. (mod(iSmAng,2) == 1)) Go To 300

      !lwj For Spherical Gaussians, batches like
      !lwj (DS|SS), (FP|SS) and (FS|PS) vanish as well
      if (ABeqCD .and. All_Spherical .and. (2*max(la,lb,lc,ld) > iSmAng)) Go To 300
      !lwj

#     ifdef _DEBUGPRINT_
      write(6,'(6A)') ' R=',ChOper(iDCRR(lDCRR)),', S=',ChOper(iDCRS(lDCRS)),', T=',ChOper(iDCRT(lDCRT))
#     endif

      kOp(3) = NrOpr(iDCRT(lDCRT))
      kOp(4) = NrOpr(ieor(iDCRT(lDCRT),iDCRS(lDCRS)))

      ix1 = 1
      iy1 = 1
      iz1 = 1
      ix2 = iPhase(1,iDCRT(lDCRT))
      iy2 = iPhase(2,iDCRT(lDCRT))
      iz2 = iPhase(3,iDCRT(lDCRT))
      lDCRE_ = 0
      lDCRT_ = iDCRT(lDCRT)

      ! Find index to desymmetrized Dij, Dkl, Dik, Dil, Djk, and
      ! Djl. Some care has to be taken here. Assume that there
      ! are two operators, T and S which generates the center
      ! pairs A,T(B) and A,S(B). If these pairs are symmetry
      ! related we will only

      if (DoFock) then
        ! Dij
        iR = iDCRR(lDCRR)
        jOp(1) = NrOpr(iR)+1
        ! Dkl
        iS = iDCRS(lDCRS)
        jOp(2) = NrOpr(iS)+1
        ! Dik
        iT = iDCRT(lDCRT)
        jOp(3) = NrOpr(iT)+1
        ! Dil
        iTS = ieor(iT,iS)
        jOp(4) = NrOpr(iTS)+1
        ! Djk
        iRT = ieor(iR,iT)
        jOp(5) = NrOpr(iRT)+1
        ! Djl
        iRTS = ieor(iRT,iS)
        jOp(6) = NrOpr(iRTS)+1

      end if ! DoFock
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Prescreening at group level

      ! Select if batch should be written on the disc at the
      ! first iteration and then later on read.

      Batch_On_Disk = (vijkl > Thize) .and. (Disc+dble(nInts+2+2/RtoI) <= Disc_Mx)

      ! Set prescreening level
      !
      ! IntOnly = T  prescreening on integral value only
      ! IntOnly = F  prescreening on integral & density matrix
      !
      ! If integral batch is written to disc, no prescreening
      ! since prescreening on integrals only was done in k2 loop.

      Prescreen_On_Int_Only = IntOnly
      if (DoIntegrals) Prescreen_On_Int_Only = .true.
      if (Batch_On_Disk) Prescreen_On_Int_Only = .true.

      ! Prescreening based on the 1st order density in AO
      ! basis (contracted). Observe that this is a rough
      ! estimate of if there will be any contribution due to
      ! this integral batch.

      ! Special care here for RS, RT and ST degeneracy.

      ! Get maximum density elements in AO basis.

      if (DoFock) then
        vij = Dij(mDij,jOp(1))
        vkl = Dkl(mDkl,jOp(2))
        vik = Dik(mDik,jOp(3))
        vil = Dil(mDil,jOp(4))
        vjk = Djk(mDjk,jOp(5))
        vjl = Djl(mDjl,jOp(6))
        ! Coulomb contributions
        Scrkl = vij*vijkl >= ThrInt
        Scrij = vkl*vijkl >= ThrInt
        DoCoul = Scrij .or. Scrkl

        ! Exchange contributions
        vik = vik/Four
        vil = vil/Four
        vjk = vjk/Four
        vjl = vjl/Four
        Scrjl = vik*vijkl >= ThrInt
        Scrjk = vil*vijkl >= ThrInt
        Scril = vjk*vijkl >= ThrInt
        Scrik = vjl*vijkl >= ThrInt
        DoExch = Scrjl .or. Scrjk .or. Scril .or. Scrik
        if (FckNoClmb) DoCoul = .false.
        if (FckNoExch) DoExch = .false.
      else
        DoCoul = .false.
        DoExch = .false.
      end if
      DoAOBatch = (DoIntegrals .and. (vijkl > CutInt)) .or. (DoFock .and. (DoCoul .or. DoExch)) .or. (Batch_On_Disk .and. W2Disc)

      ! Branch out if crude estimate indicates no contributions!

      if (.not. DoAOBatch) then
        ! AO batch is not on the disk
        if (.not. Batch_On_Disk) Go To 300
        if (.not. W2Disc) then

          ! AO batch is on disk! Do a no copy read to
          ! position the next batch on the disc.

1111      continue
          call iRBuf(iWR,2,Copy)
          call dRBuf(QInd,2,Copy)
          call Store_QLast(QInd)
          kInts = iWR(1)
          mInts = iWR(2)
          if ((QInd(1) == Quad_ijkl) .and. (QInd(2) == RST_triplet)) then
            if (kInts /= nInts) then
              call WarningMessage(2,'Twoel: kInts /= nInts!')
              write(6,*) 'Twoel: kInts,mInts,nInts=',kInts,mInts,nInts
              write(6,*) 'Index,1:',QInd(1),QInd(2),Quad_ijkl,RST_triplet
              call Abend()
            end if
            if (mInts > 0) call dRBuf(Wrk(iW3),mInts,NoCopy)
            Disc = Disc+dble(2/RtoI+2+mInts)
            Go To 300
          else if (QInd(1) <= Quad_ijkl) then
            if (mInts > 0) call dRBuf(Wrk(iW3),mInts,NoCopy)
            Disc = Disc+dble(2/RtoI+2+mInts)
            Go To 1111
          else
            call WarningMessage(2,'Twoel: batch is lost!')
            write(6,*) 'Index,1:',QInd(1),QInd(2),Quad_ijkl,RST_triplet
            call Abend()
          end if
        end if
      end if

      if (DoFock) then
        if (iShell(1) >= iShell(2)) then
          ij1 = iBasi
          ij2 = jBasj
          ij3 = iCmp(1)
          ij4 = iCmp(2)
        else
          ij1 = jBasj
          ij2 = iBasi
          ij3 = iCmp(2)
          ij4 = iCmp(1)
        end if
        if (iShell(3) >= iShell(4)) then
          kl1 = kBask
          kl2 = lBasl
          kl3 = iCmp(3)
          kl4 = iCmp(4)
        else
          kl1 = lBasl
          kl2 = kBask
          kl3 = iCmp(4)
          kl4 = iCmp(3)
        end if
        if (iShell(1) >= iShell(3)) then
          ik1 = iBasi
          ik2 = kBask
          ik3 = iCmp(1)
          ik4 = iCmp(3)
        else
          ik1 = kBask
          ik2 = iBasi
          ik3 = iCmp(3)
          ik4 = iCmp(1)
        end if
        if (iShell(1) >= iShell(4)) then
          il1 = iBasi
          il2 = lBasl
          il3 = iCmp(1)
          il4 = iCmp(4)
        else
          il1 = lBasl
          il2 = iBasi
          il3 = iCmp(4)
          il4 = iCmp(1)
        end if
        if (iShell(2) >= iShell(3)) then
          jk1 = jBasj
          jk2 = kBask
          jk3 = iCmp(2)
          jk4 = iCmp(3)
        else
          jk1 = kBask
          jk2 = jBasj
          jk3 = iCmp(3)
          jk4 = iCmp(2)
        end if
        if (iShell(2) >= iShell(4)) then
          jl1 = jBasj
          jl2 = lBasl
          jl3 = iCmp(2)
          jl4 = iCmp(4)
        else
          jl1 = lBasl
          jl2 = jBasj
          jl3 = iCmp(4)
          jl4 = iCmp(2)
        end if
      end if ! DoFock

      ! Branch point for partial integral storage.

      if (Batch_On_Disk .and. (.not. W2Disc)) Go To 6767
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Here if the AO batch will be computed !

      ! Compute actual size of the {a0|c0} block

      mabMin = nabSz(max(la,lb)-1)+1
      if (EQ(CoorM(1,1),CoorM(1,2))) mabMin = nabSz(la+lb-1)+1
      mabMax = nabSz(la+lb)
      mcdMin = nabSz(max(lc,ld)-1)+1
      if (EQ(CoorM(1,3),CoorM(1,4))) mcdMin = nabSz(lc+ld-1)+1
      mcdMax = nabSz(lc+ld)
      mabcd = (mabMax-mabMin+1)*(mcdMax-mcdMin+1)

      ! Find the proper centers to start of with the angular
      ! momentum on. If la == lb there will exist an
      ! ambiguity to which center that angular momentum should
      ! be accumulated on. In that case we will use A and C of
      ! the order as defined by the basis functions types.

      if (iAnga(1) >= iAnga(2)) then
        call dcopy_(3,CoorM(1,1),1,CoorAC(1,1),1)
      else
        call dcopy_(3,CoorM(1,2),1,CoorAC(1,1),1)
      end if
      if (iAnga(3) >= iAnga(4)) then
        call dcopy_(3,CoorM(1,3),1,CoorAC(1,2),1)
      else
        call dcopy_(3,CoorM(1,4),1,CoorAC(1,2),1)
      end if

      ! Set flags if triangularization will be used

      IeqK = EQ(CoorM(1,1),CoorM(1,3))
      JeqL = EQ(CoorM(1,2),CoorM(1,4))
      IJeqKL = IeqK .and. JeqL

      ! Loops to partition the primitives

      ! Reset pointer ipAOInt if we need to reserve speacial
      ! space for the contracted integrals.

      IncZet = nAlpha*jPrInc
      IncEta = nGamma*lPrInc
      if ((nZeta /= IncZet) .or. (nEta /= IncEta)) then
        mWork2 = nWork2-nijkl*mabcd
        ipAOInt = 1+nijkl*mabcd
      else
        mWork2 = nWork2
        ipAOInt = 1
      end if

      nZeta_Tot = k2Data1(lDCR1)%IndZ(nZeta+1)
      nEta_Tot = k2Data2(lDCR2)%IndZ(nEta+1)

      kabcd = 0
      Do_TnsCtl = .false.
      NoInts = .true.
      NoPInts = .true.
      ipAOInt_ = ipAOInt
      iW4_ = iW4

      do iZeta=1,nZeta_Tot,IncZet
        mZeta = min(IncZet,nZeta_Tot-iZeta+1)
        if (lEmpty(Coeff2,nBeta,nBeta,jBasj)) cycle

        do iEta=1,nEta_Tot,IncEta
          mEta = min(IncEta,nEta_Tot-iEta+1)
          if (lEmpty(Coeff4,nDelta,nDelta,lBasl)) cycle

          call DrvRys(iZeta,iEta,nZeta,nEta,mZeta,mEta,nZeta_Tot,nEta_Tot,k2data1(lDCR1),k2data2(lDCR2),nAlpha,nBeta,nGamma, &
                      nDelta,ix1,iy1,iz1,ix2,iy2,iz2,ThrInt,CutInt,vij,vkl,vik,vil,vjk,vjl,Prescreen_On_Int_Only,NoInts,iAnga, &
                      CoorM,CoorAC,mabMin,mabMax,mcdMin,mcdMax,nijkl/nComp,nabcd,mabcd,Wrk,ipAOInt_,iW4_,nWork2,mWork2, &
                      k2data1(lDCR1)%HrrMtrx(:,NrOpr(lDCRE_)+1),k2data2(lDCR2)%HrrMtrx(:,NrOpr(lDCRT_)+1),la,lb,lc,ld,iCmp,iShll, &
                      NoPInts,Dij(1,jOp(1)),mDij,Dkl(1,jOp(2)),mDkl,Do_TnsCtl,kabcd,Coeff1,iBasi,Coeff2,jBasj,Coeff3,kBask,Coeff4, &
                      lBasl)

        end do
      end do

      ! Integrals are now returned in Wrk(ipAOInt) or Wrk(iW4)

      if (NoPInts) then
        if (W2Disc) then
          if (Batch_On_Disk) then

            ! If the batch was supposed to be on disk make a mark.

            mInts = 0

            iWR(1) = nInts
            iWR(2) = mInts
            call iWBuf(iWR,RtoI)
            call dWBuf(QInd,2)
            call Store_QLast(QInd)

            Disc = Disc+dble(2/RtoI+2+mInts)
          end if
        end if
        Go To 300
      end if

      ! Multiply with factors due to summation over DCR

      if (MolWgh == 1) then
        FactNd = dble(nIrrep)/dble(LmbdT)
      else if (MolWgh == 0) then
        FactNd = u*v*w*x/dble(nIrrep**3*LmbdT)
      else
        FactNd = sqrt(u*v*w*x)/dble(nirrep*lmbdt)
      end if
      if (FactNd /= One) call DScal_(kabcd*nijkl,FactNd,Wrk(iW4),1)

      ! Apply the transfer equation and transform the spherical
      ! harmonic gaussian.

      if (Do_TnsCtl) then
        call TnsCtl(Wrk(iW4),nWork2,nijkl,mabMax,mabMin,mcdMax,mcdMin,k2data1(lDCR1)%HrrMtrx(:,NrOpr(lDCRE_)+1), &
                    k2data2(lDCR2)%HrrMtrx(:,NrOpr(lDCRT_)+1),la,lb,lc,ld,iCmp(1),iCmp(2),iCmp(3),iCmp(4),iShll(1),iShll(2), &
                    iShll(3),iShll(4),i_Int)
        ipAOInt = i_Int
        if (i_Int == 1) then
          iW3 = 1+nijkl*nabcd
        else
          iW3 = 1
        end if
      else

        ! Undo the late Cntrct

        call dcopy_(nabcd*nijkl,Wrk(ipAOInt),1,Wrk(iW3),1)
        call DGeTMO(Wrk(iW3),nabcd,nabcd,nijkl,Wrk(ipAOInt),nijkl)

      end if

      ! Branch point for partial integral storage

      if (Batch_On_Disk .and. W2Disc) then

        ! Write integrals to current position on disc.

        iOpt = 0 ! Always Packing, not run length
        call PkR8(iOpt,nInts,nByte,Wrk(ipAOInt),Wrk(iW3))
        mInts = (nByte+RtoB-1)/RtoB

        iWR(1) = nInts
        iWR(2) = mInts
        call iWBuf(iWR,2)
        call dWBuf(QInd,2)
        call Store_QLast(Qind)
        call dWBuf(Wrk(iW3),mInts)
        Disc = Disc+dble(2/RtoI+2+mInts)

      end if

6767  continue
      if (Batch_On_Disk .and. (.not. W2Disc)) then
1112    continue
        call iRBuf(iWR,2,Copy)
        call dRBuf(QInd,2,Copy)
        call Store_QLast(QInd)
        kInts = iWR(1)
        mInts = iWR(2)
        if ((QInd(1) == Quad_ijkl) .and. (QInd(2) == RST_triplet)) then
          if (kInts /= nInts) then
            call WarningMessage(2,'Twoel: kInts /= nInts!')
            write(6,*) 'Twoel: kInts,mInts,nInts=',kInts,mInts,nInts
            write(6,*) 'Index,1:',QInd(1),QInd(2),Quad_ijkl,RST_triplet
            call Abend()
          end if
          if (mInts > 0) call dRBuf(Wrk(iW3),mInts,Copy)
          Disc = Disc+dble(2/RtoI+2+mInts)
          if (mInts == 0) Go To 300
        else if (QInd(1) < Quad_ijkl) then
          if (mInts > 0) call dRBuf(Wrk(iW3),mInts,NoCopy)
          Disc = Disc+dble(2/RtoI+2+mInts)
          Go To 1112
        else
          call WarningMessage(2,'Twoel: batch is lost!')
          write(6,*) 'Index,1:',QInd(1),QInd(2),Quad_ijkl,RST_triplet
          call Abend()
        end if

        iOpt = 0 ! Always packing, not run length
        call UpkR8(iOpt,nInts,nByte,Wrk(iW3),Wrk(ipAOInt))
      end if

      ! Accumulate contributions directly to the symmetry
      ! adapted Fock matrix.

      mWork3 = nWork2-iW3+1
      if (DoFock) call FckAcc(iAnga,iCmp,Shijij,iShll,iShell,kOp,nijkl/nComp,Wrk(ipAOInt),TwoHam,size(TwoHam),Wrk(iW3),mWork3,iAO, &
                              iAOst,iBasi,jBasj,kBask,lBasl,Dij(1,jOp(1)),ij1,ij2,ij3,ij4,Dkl(1,jOp(2)),kl1,kl2,kl3,kl4, &
                              Dik(1,jOp(3)),ik1,ik2,ik3,ik4,Dil(1,jOp(4)),il1,il2,il3,il4,Djk(1,jOp(5)),jk1,jk2,jk3,jk4, &
                              Djl(1,jOp(6)),jl1,jl2,jl3,jl4,FckTmp,nFT,DoCoul,DoExch,ExFac)

      ! Transform from AO basis to SO basis

      if (DoIntegrals) call SymAdp(iAnga,iCmp(1),iCmp(2),iCmp(3),iCmp(4),Shijij,iShll,iShell,iAO,kOp,nijkl,Aux,nAux,Wrk(ipAOInt), &
                                   SOInt,nSOInt,NoInts)

300   continue
    end do
200 continue
  end do
100 continue
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(iS_)
  call Unused_integer(jS_)
  call Unused_integer(kS_)
  call Unused_integer(lS_)
  call Unused_integer(iPrInc)
  call Unused_integer(kPrInc)
  call Unused_real_array(Dens)
end if

end subroutine TwoEl_Sym
