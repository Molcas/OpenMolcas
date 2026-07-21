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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine OLagNS_Hel2(iCase,NBSQT,lT2AO,iSym,iSymA,iSymB,iSymI,iSymJ,nMaxOrb,ERI1,Amp1,Scr,DPT2C,T2AO)
! DMNS_{ijkl}*d(ij|kl)/dx -> (pj|kl)*D_{qjkl} + (ip|kl)*D_{iqkl}
!                          + (ij|pl)*D_{ijql} + (ij|kp)*D_{ijkq}

! Integrals needed:
! (OO|OO), (VO|OO), (VV|OO), (VO|VO), (VV|VO)
! -> <OO|OO>, <VO|OO>, <VO|VO>, <VV|OO>, <VV|VO>
! Here, O is occupied (doubly and partially) orbitals
!       V is not filled (partially and virtual) orbitals

! <**|VV> and <**|OO> are fetched by Exch
! <V*|V*> and <O*|O*> are fetched by Exch

! However, (*O|*O) is split into (*C|*C), (*C|*A), (*A|*A)
!          (V*|V*) is split into (A*|A*), (A*|V*), (V*|V*)
!          (*O|V*) is split into (*A|A*), (*A|V*)

use Symmetry_Info, only: Mul
use SUPERINDEX, only: KAGEB, KAGTB, KIGEJ, KIGTJ, KTGEU, KTGTU, KTU, KTUV
use EQSOLV, only: IVECC2
use fake_GA, only: GA_Arrays
use caspt2_global, only: OLag
use general_data, only: NACTEL, NASH
use caspt2_module, only: HZERO, NAES, NAGEB, NAGEBES, NAGTB, NAGTBES, NASUP, NBAS, NBAST, NFRO, NIES, NIGEJES, &
                         NIGTJES, NISH, NISUP, NSES, NSSH, NSYM, NTGEUES, NTGTUES, NTU, NTUES, NTUVES
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iCase, NBSQT, lT2AO, iSym, iSymA, iSymB, iSymI, iSymJ, nMaxOrb
real(kind=wp), intent(out) :: ERI1(NBSQT), Amp1(nMaxOrb,nMaxOrb), Scr(nMaxOrb,nMaxOrb)
real(kind=wp), intent(inout) :: DPT2C(NBSQT), T2AO(lT2AO)
integer(kind=iwp) :: iAabs, iAgeB, iAgtB, iAS, iASM, iASP, iAtot, iBabs, iBtot, IgeJ, IgtJ, iIabs, iIS, iISM, iISP, iItot, iJabs, &
                     iJtot, IO1, IO2, IOFF1(8), IOFF2(8), ipTC, ipTCM, ipTCP, iSymAB, iSymK, iTabs, iUabs, iVabs, iVaHM, iVaHP, &
                     iVaM, iVaP, iVHM, iVHP, iViHM, iViHM0, iViHP, iViHP0, iViM, iViP, iVjM, iVjP, IW1, nAS, nAshA, nAshB, nAshI, &
                     nAshJ, nASM, nASP, nBasA, nBasB, nBasI, nBasJ, nCorA, nCorB, nCorI, nCorJ, nFroA, nFroB, nFroI, nFroJ, nIS, &
                     nIshA, nIshB, nIshI, nIshJ, nISM, nISP, nJ, nOccA, nOccA2, nOccB, nOccB2, nOrbA, nSshA, nSshB
real(kind=wp) :: Fac, ONEADD, SQ2, SQ3, SQI2, ValA, ValBM, ValBP, ValC1, ValC2, ValD1, ValD2, ValEM, ValEP, ValFM, ValFP, ValGM, &
                 ValGP, ValHM, ValHP
logical(kind=iwp) :: PM
real(kind=wp), allocatable :: WRK1(:), WRK2(:)

! EXCH(ISYP,ISYI,ISYQ,ISYJ,II,IJ,ERI,SCR)

! rhs_mp2_help1

!! The amplitude is in the IC (internally contracted) basis, so
!! the active orbital index (indices) must be transformed to the
!! (quasi-)canonical MO (or contravatiant) basis.
!! IC = SR (why?), contravariant = C
!write(u6,*) 'icase = ',icase

if ((iCase == 2) .or. (iCase == 6) .or. (iCase == 8) .or. (iCase == 10) .or. (iCase == 12)) then
  PM = .true.
else
  PM = .false.
end if
if ((iCase == 3) .or. (iCase == 7) .or. (iCase == 9) .or. (iCase == 11) .or. (iCase == 13)) return

SQ2 = sqrt(Two)
SQI2 = One/SQ2
SQ3 = sqrt(Three)
!iVec = iVecX
IO1 = 0
IO2 = 0
do iSymK=1,nSym
  IOFF1(iSymK) = IO1
  IOFF2(iSymK) = IO2
  iSymAB = Mul(iSymK,iSym)
  IO1 = IO1+nIsh(iSymK)*nAgeB(iSymAB)
  IO2 = IO2+nIsh(iSymK)*nAgtB(iSymAB)
end do

!! Some setup
!! Read T-amplitude, hopefully in contravariant form
!nINP = 0
!nINM = 0
!nIN = 0
if (PM) then
  !nINP = nINDEP(iSym,iCase)
  nASP = nASup(iSym,iCase)
  nISP = nISup(iSym,iCase)
  !if (nINP /= 0) nVec = nINP*nISP
  !nINM = nINDEP(iSym,iCase+1)
  nASM = nASup(iSym,iCase+1)
  nISM = nISup(iSym,iCase+1)
  !if (nINM /= 0) nVec = nINM*nISM
  if (nASP*nISP /= 0) then
    call RHS_ALLO(nASP,nISP,ipTCP)
    call RHS_READ(nASP,nISP,ipTCP,iCase,iSym,iVecC2)
  end if
  if (nASM*nISM /= 0) then
    call RHS_ALLO(nASM,nISM,ipTCM)
    call RHS_READ(nASM,nISM,ipTCM,iCase+1,iSym,iVecC2)
  end if
else
  !nIN = nINDEP(iSym,iCase)
  nAS = nASup(iSym,iCase)
  nIS = nISup(iSym,iCase)
  !if (nIN /= 0) nVec = nIN*nIS
  if (nAS*nIS /= 0) then
    call RHS_ALLO(nAS,nIS,ipTC)
    call RHS_READ(nAS,nIS,ipTC,iCase,iSym,iVecC2)
  end if
end if

! finish if no contributions
! should we use the number of independent vectors?
if (PM) then
  if ((nASP*nISP == 0) .and. (nASM*nISM == 0)) return
else
  if (nAS*nIS == 0) return
end if

call mma_allocate(WRK1,nBasT*nBasT,Label='WRK1')
call mma_allocate(WRK2,nBasT*nBasT,Label='WRK2')

nFroI = nFro(iSymI)
nFroJ = nFro(iSymJ)
nFroA = nFro(iSymA)
nFroB = nFro(iSymB)
nIshI = nIsh(iSymI)
nIshJ = nIsh(iSymJ)
nIshA = nIsh(iSymA)
nIshB = nIsh(iSymB)
nAshI = nAsh(iSymI)
nAshJ = nAsh(iSymJ)
nAshA = nAsh(iSymA)
nAshB = nAsh(iSymB)
nSshA = nSsh(iSymA)
nSshB = nSsh(iSymB)
nBasI = nBas(iSymI)
nBasJ = nBas(iSymJ)
nBasA = nBas(iSymA)
nBasB = nBas(iSymB)

nCorI = nFroI+nIshI
nCorJ = nFroJ+nIshJ
nCorA = nFroA+nIshA
nCorB = nFroB+nIshB
nOccA = nCorA+nAshA
nOccB = nCorB+nAshB
nOccA2 = nOccA-nFroA
nOccB2 = nOccB-nFroB
nOrbA = nOccA+nSshA

!nOcc = nFro(iSym)+nIsh(iSym)+nAsh(iSym)

!! active+virtual part for the right index
!nJ = nFro(iSymJ)+nIsh(iSymJ)+nAsh(iSymJ)

if (iCase == 1) then
  call OLagNS_A(Amp1)
else if ((iCase == 2) .or. (iCase == 3)) then
  call OLagNS_B(Amp1)
else if (iCase == 4) then
  call OLagNS_C(Amp1)
else if (iCase == 5) then
  call OLagNS_D(Amp1)
else if ((iCase == 6) .or. (iCase == 7)) then
  call OLagNS_E(Amp1)
else if ((iCase == 8) .or. (iCase == 9)) then
  call OLagNS_F(Amp1)
else if ((iCase == 10) .or. (iCase == 11)) then
  call OLagNS_G(Amp1)
else if ((iCase == 12) .or. (iCase == 13)) then
  call OLagNS_H(Amp1)
end if

if ((iCase == 1) .and. (HZERO == 'DYALL')) call OLagNS_AAAA(Amp1)
call mma_deallocate(WRK1)
call mma_deallocate(WRK2)

if (PM) then
  if (nASP*nISP /= 0) call RHS_FREE(ipTCP)
  if (nASM*nISM /= 0) call RHS_FREE(ipTCM)
else
  if (nAS*nIS /= 0) call RHS_FREE(ipTC)
end if

return

contains

subroutine OLagNS_A(AmpL1)

  real(kind=wp), intent(out) :: AmpL1(nAshA,nAshB)
  integer(kind=iwp) :: iA, iB, iI, iJ

  if (nAshI*nIshJ*nAshA*nAshB == 0) return

  !nJ = nIshJ
  do iI=1,nAshI
    iIabs = iI+nIshI+nAes(iSymI)
    !iItot = iI + nCorI
    !if (iSymI == iSymJ) nJ = iI
    do iJ=1,nIshJ
      iJabs = iJ+nIes(iSymJ)
      iJtot = iJ+nFroJ
      !if (iIabs < iJabs) Cycle
      Fac = One
      !if ((iI /= iJ) .and. (iSymI == iSymJ)) Fac = Two
      if (iSymI /= iSymJ) Fac = Two

      call Exch(iSymA,iSymI,iSymB,iSymJ,iI+nCorI,iJ+nFroJ,ERI1,Scr)

      AmpL1(:,:) = Zero

      do iA=1,nAshA
        iAabs = iA+nAes(iSymA)
        !iAtot = iA+nCorA
        do iB=1,nAshB
          iBabs = iB+nAes(iSymB)
          iBtot = iB+nCorB

          iTabs = iBabs
          iUabs = iI+nAes(iSymI)
          iVabs = iAabs
          IW1 = kTUV(iTabs,iUabs,iVabs)-nTUVes(iSym)
          !ValA = Zero
          !do iICB=1,nIN
          !  iVA = iICB+nIN*(iJabs-1)
          !  ValA = ValA+Work(ipT+iVA-1)*Work(LST+IW1-1+nAS*(iICB-1))
          !end do
          !ValA = ValA*Two
          iIS = iJabs
          iAS = IW1
          ValA = GA_Arrays(ipTC)%A(iAS+nAS*(iIS-1))*Two

          !! For FIMO derivative
          if (iUabs == iVabs) DPT2C(iBtot+nOrbA*(iJtot-1)) = DPT2C(iBtot+nOrbA*(iJtot-1))+ValA

          AmpL1(iA,iB) = AmpL1(iA,iB)+ValA
        end do
      end do

      AmpL1(:,:) = AmpL1(:,:)*Fac

      !! Calculate the actual contributions
      !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
      call OLagNS_post1(1,1,ERI1,AmpL1)

      !! Prepare for implicit (VV|VO) integrals
      !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
      !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu}
      call OLagNS_post2(nAshA,nAshB,nCorA,nCorB,AmpL1,WRK2)
      !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
      call OLagNS_post3(iIabs,iJabs,T2AO,WRK2)
    end do
  end do

end subroutine OLagNS_A

subroutine OLagNS_B(AmpL1)

  real(kind=wp), intent(out) :: AmpL1(nAshA,nAshB)
  integer(kind=iwp) :: iA, iB, iI, iJ

  if (nIshI*nIshJ*nAshA*nAshB == 0) return

  nJ = nIshJ
  do iI=1,nIshI
    iIabs = iI+nIes(iSymI)
    !iItot = iI+nFroI
    if (iSymI == iSymJ) nJ = iI
    do iJ=1,nJ
      iJabs = iJ+nIes(iSymJ)
      !iJtot = iJ+nFroJ
      if (iIabs < iJabs) cycle
      Fac = One
      !if ((iI /= iJ) .and. (iSymI == iSymJ)) Fac = Two
      if (iSymI /= iSymJ) Fac = Two

      call Exch(iSymA,iSymI,iSymB,iSymJ,iI+nFroI,iJ+nFroJ,ERI1,Scr)

      AmpL1(:,:) = Zero

      do iA=1,nAshA
        iAabs = iA+nAes(iSymA)
        !iAtot = iA+nCorA
        do iB=1,nAshB
          iBabs = iB+nAes(iSymB)
          !iBtot = iB+nCorB

          if (iaabs > ibabs) then
            iTabs = iAabs
            iUabs = iBabs
          else
            iTabs = iBabs
            iUabs = iAabs
          end if
          iViP = kIgeJ(iIabs,iJabs)-nIgeJes(iSym) ! inactive
          iVaP = kTgeU(iTabs,iUabs)-nTgeUes(iSym) !   active
          iViM = kIgtJ(iIabs,iJabs)-nIgtJes(iSym)
          iVaM = kTgtU(iTabs,iUabs)-nTgtUes(iSym)
          !! transform internally contracted (SR)
          !!        to contravariant (C)
          ValBP = Zero
          ValBM = Zero
          !do iICB=1,nINP
          !  iVP  = iICB+nINP*(iViP-1)
          !  ValBP = ValBP+Work(ipTP+iVP-1)*Work(LSTP+iVaP-1+nASP*(iICB-1))
          !end do
          iIS = iViP
          iAS = iVaP
          ValBP = GA_Arrays(ipTCP)%A(iAS+nASP*(iIS-1))
          if ((iAabs /= iBabs) .and. (iIabs /= iJabs)) then
            if (iIabs /= iJabs) then
              !do iICB=1,nINM
              !  iVM = iICB+nINM*(iViM-1)
              !  ValBM = ValBM+Work(ipTM+iVM-1)*Work(LSTM+iVaM-1+nASM*(iICB-1))
              !end do
              iIS = iViM
              iAS = iVaM
              ValBM = GA_Arrays(ipTCM)%A(iAS+nASM*(iIS-1))
            end if
            !! permutated
            if (iAabs < iBabs) ValBM = -ValBM
          end if
          if (iIabs == iJabs) ValBP = ValBP*SQI2

          AmpL1(iA,iB) = AmpL1(iA,iB)+ValBP+ValBM
        end do
      end do

      AmpL1(:,:) = AmpL1(:,:)*Fac

      !! Calculate the actual contributions
      !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
      call OLagNS_post1(1,1,ERI1,AmpL1)

      !! Prepare for implicit (VV|VO) integrals
      !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
      !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu}
      call OLagNS_post2(nAshA,nAshB,nCorA,nCorB,AmpL1,WRK2)
      !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
      call OLagNS_post3(iIabs,iJabs,T2AO,WRK2)
    end do
  end do

end subroutine OLagNS_B

subroutine OLagNS_C(AmpL1)

  real(kind=wp), intent(out) :: AmpL1(nAshA+nSshA,nAshB+nSshB)
  integer(kind=iwp) :: iA, iB, iI, iJ, iXabs

  if (nAshI*nAshJ*nSshA*nAshB == 0) return

  !nJ = nIshJ
  do iI=1,nAshI
    iIabs = iI+nIshI+nAes(iSymI)
    !iItot = iI+nCorI
    !if (iSymI == iSymJ) nJ = iI
    do iJ=1,nAshJ
      iJabs = iJ+nIshJ+nAes(iSymJ)
      !iJtot = iJ+nCorJ
      if (iIabs < iJabs) cycle
      Fac = One
      !if ((iI /= iJ) .and. (iSymI == iSymJ)) Fac = Two
      if (iSymI /= iSymJ) Fac = Two

      call Exch(iSymA,iSymI,iSymB,iSymJ,iI+nCorI,iJ+nCorJ,ERI1,Scr)

      AmpL1(:,:) = Zero

      do iA=1,nSshA
        iAabs = iA+nSes(iSymA)
        iAtot = iA+nFroA+nIshA+nAshA
        do iB=1,nAshB
          iBabs = iB+nAes(iSymB)
          iBtot = iB+nCorB

          iTabs = iI+nAes(iSymI)
          iUabs = iBabs
          iVabs = iJ+nAes(iSymJ)
          !write(u6,*) itabs,iuabs,ivabs
          !! (at|uv) -> (ai|bj) -> (at|uv)
          !IW1 = kTUV(iTabs,iUabs,iVabs)-nTUVes(iSym)
          ValC1 = Zero
          ValC2 = Zero
          !do iICB=1,nIN
          !  iV = iICB+nIN*(iAabs-1)
          !  ValC1 = ValC1+Work(ipT+iV-1)*Work(LST+IW1-1+nAS*(iICB-1))
          !end do
          !if (iIabs /= iJabs) then
          !  IW1 = kTUV(iVabs,iUabs,iTabs)-nTUVes(iSym)
          !  do iICB=1,nIN
          !    iV = iICB+nIN*(iAabs-1)
          !    ValC2 = ValC2+Work(ipT+iV-1)*Work(LST+IW1-1+nAS*(iICB-1))
          !  end do
          !end if

          !ValC1 = ValC1*Two
          !ValC2 = ValC2*Two

          iIS = iAabs
          iAS = kTUV(iTabs,iUabs,iVabs)-nTUVes(iSym)
          ValC1 = GA_Arrays(ipTC)%A(iAS+nAS*(iIS-1))*Two
          if (iIabs /= iJabs) then
            iAS = kTUV(iVabs,iUabs,iTabs)-nTUVes(iSym)
            ValC2 = GA_Arrays(ipTC)%A(iAS+nAS*(iIS-1))*Two
          end if

          iTabs = iBabs
          iUabs = iI+nAes(iSymI)
          iVabs = iJ+nAes(iSymJ)
          if (iUabs == iVabs) then
            !! For FIMO derivative
            ONEADD = Zero
            !IW1 = kTUV(iTabs,iUabs,iVabs)-nTUVes(iSym)
            !do iICB=1,nIN
            !  iV = iICB+nIN*(iAabs-1)
            !  ONEADD = ONEADD+Work(ipT+iV-1)*Work(LST+IW1-1+nAS*(iICB-1))
            !end do
            !ONEADD = ONEADD*Two
            iAS = kTUV(iTabs,iUabs,iVabs)-nTUVes(iSym)
            ONEADD = GA_Arrays(ipTC)%A(iAS+nAS*(iIS-1))*Two
            DPT2C(iAtot+nOrbA*(iBtot-1)) = DPT2C(iAtot+nOrbA*(iBtot-1))+ONEADD

            !! For -sum(y)(ay,yt) -> (ay,ty) derivative
            !! It is correct, but should be rewritten
            !ONEADD = Zero
            !do iXabs=1,nAshI !?
            !  IW1 = kTUV(iTabs,iXabs,iXabs)-nTUVes(iSym)
            !  do iICB=1,nIN
            !    iV = iICB+nIN*(iAabs-1)
            !    ONEADD = ONEADD+Work(ipT+iV-1)*Work(LST+IW1-1+nAS*(iICB-1))
            !  end do
            !end do
            !ONEADD = Two*ONEADD/real(MAX(1,NACTEL),kind=wp)
            ONEADD = Zero
            do iXabs=1,nAshI !?
              iAS = kTUV(iTabs,iXabs,iXabs)-nTUVes(iSym)
              ONEADD = ONEADD+GA_Arrays(ipTC)%A(iAS+nAS*(iIS-1))
            end do
            ONEADD = Two*ONEADD/real(max(1,NACTEL),kind=wp)
            AmpL1(iAtot-nCorA,iBtot-nCorA) = AmpL1(iAtot-nCorA,iBtot-nCorA)-ONEADD
          end if
          AmpL1(iAtot-nCorA,iBtot-nCorB) = AmpL1(iAtot-nCorA,iBtot-nCorB)+ValC1
          AmpL1(iBtot-nCorB,iAtot-nCorA) = AmpL1(iBtot-nCorB,iAtot-nCorA)+ValC2
        end do
      end do

      AmpL1(:,:) = AmpL1(:,:)*Fac

      !! Calculate the actual contributions
      !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
      call OLagNS_post1(3,3,ERI1,AmpL1)

      !! Prepare for implicit (VV|VO) integrals
      !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
      !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu}
      call OLagNS_post2(nAshA+nSshA,nAshB+nSshB,nCorA,nCorB,AmpL1,WRK2)
      !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
      call OLagNS_post3(iIabs,iJabs,T2AO,WRK2)
    end do
  end do

end subroutine OLagNS_C

subroutine OLagNS_D(AmpL1)

  real(kind=wp), intent(out) :: AmpL1(nAshA+nSshA,nAshB+nSshB)
  integer(kind=iwp) :: iA, iB, iI, iJ

  if (nAshI*nIshJ*nSshA*nAshB == 0) return

  !nJ = nIshJ
  do iI=1,nAshI
    iIabs = iI+nIshI+nAes(iSymI)
    iItot = iI+nCorI
    !if (iSymI == iSymJ) nJ = iI
    do iJ=1,nIshJ
      iJabs = iJ+nIes(iSymJ)
      iJtot = iJ+nFroJ
      !if (iIabs < iJabs) cycle
      Fac = One
      !if ((iI /= iJ) .and. (iSymI == iSymJ)) Fac = Two
      if (iSymI /= iSymJ) Fac = Two

      call Exch(iSymA,iSymI,iSymB,iSymJ,iI+nCorI,iJ+nFroJ,ERI1,Scr)

      AmpL1(:,:) = Zero

      do iA=1,nSshA
        iAabs = iA+nSes(iSymA)
        iAtot = iA+nFroA+nIshA+nAshA
        do iB=1,nAshB
          !iBabs = iB+nAes(iSymB)
          iBtot = iB+nCorB

          !iVi = iJabs+nIshA*(iAabs-1)+IOFF1(iSymA)
          !iVa1 = iB+nAshB*(iI+nAes(iSymI)-1)+IOFF1(iSymB)
          !iVa2 = iB+nAshB*(iI+nAes(iSymI)-1)+IOFF1(iSymB)+nAshT*nAshT
          !! transform internally contracted (SR)
          !!        to contravariant (C)
          !ValD1 = Zero
          !ValD2 = Zero
          !do iICB=1,nIN
          !  iVD = iICB+nIN*(iVi-1)
          !  ValD1 = ValD1+Work(ipT+iVD-1)*Work(LST+iVa1-1+nAS*(iICB-1))
          !  ValD2 = ValD2+Work(ipT+iVD-1)*Work(LST+iVa2-1+nAS*(iICB-1))
          !end do
          !ValD1 = ValD1*Two
          !ValD2 = ValD2*Two
          iIS = iJabs+nIshA*(iAabs-1)+iOFF1(iSymA)
          iAS = kTU(iB,iI)-nTUes(iSymA)
          ValD1 = GA_Arrays(ipTC)%A(iAS+nAS*(iIS-1))*Two
          iAS = iAS+nTU(iSymA)
          ValD2 = GA_Arrays(ipTC)%A(iAS+nAS*(iIS-1))*Two

          !! Fock contributions from the inactive density
          if (iItot == iBtot) DPT2C(iAtot+nOrbA*(iJtot-1)) = DPT2C(iAtot+nOrbA*(iJtot-1))+ValD1

          AmpL1(iA+nAshA,iB) = AmpL1(iA+nAshA,iB)+ValD2
          AmpL1(iB,iA+nAshB) = AmpL1(iB,iA+nAshB)+ValD1
        end do
      end do

      AmpL1(:,:) = AmpL1(:,:)*Fac

      !! Calculate the actual contributions
      !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
      call OLagNS_post1(3,3,ERI1,AmpL1)

      !! Prepare for implicit (VV|VO) integrals
      !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
      !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu}
      call OLagNS_post2(nAshA+nSshA,nAshB+nSshB,nCorA,nCorB,AmpL1,WRK2)
      !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
      call OLagNS_post3(iIabs,iJabs,T2AO,WRK2)
    end do
  end do

end subroutine OLagNS_D

subroutine OLagNS_E(AmpL1)

  real(kind=wp), intent(out) :: AmpL1(nAshA+nSshA,nAshB+nSshB)
  integer(kind=iwp) :: iA, iB, iI, iJ

  if (nIshI*nIshJ*nSshA*nAshB == 0) return

  nJ = nIshJ
  do iI=1,nIshI
    iIabs = iI+nIes(iSymI)
    !iItot = iI
    if (iSymI == iSymJ) nJ = iI
    do iJ=1,nJ
      iJabs = iJ+nIes(iSymJ)
      !iJtot = iJ
      if (iIabs < iJabs) cycle
      Fac = One
      !if ((iI /= iJ) .and. (iSymI == iSymJ)) Fac = Two
      if (iSymI /= iSymJ) Fac = Two

      call Exch(iSymA,iSymI,iSymB,iSymJ,iI+nFroI,iJ+nFroJ,ERI1,Scr)

      AmpL1(:,:) = Zero

      IgeJ = kIgeJ(iIabs,iJabs)-nIgeJes(iSym) ! iSymIJ
      IgtJ = kIgtJ(iIabs,iJabs)-nIgtJes(iSym) ! iSymIJ
      do iA=1,nSshA
        iAabs = iA+nSes(iSymA)
        !iAtot = iA+nIshA+nAshA
        do iB=1,nAshB
          iBabs = iB+nAes(iSymB)
          !iBtot = iB+nIshB

          iASP = iBabs
          iISP = iAabs+nSshA*(IgeJ-1)+iOFF1(iSymA)
          ValEP = GA_Arrays(ipTCP)%A(iASP+nASP*(iISP-1))
          ValEM = Zero
          if (iIabs > iJabs) then
            ValEP = ValEP*SQ2
            iASM = iBabs
            iISM = iAabs+nSshA*(IgtJ-1)+iOFF1(iSymA)
            ValEM = GA_Arrays(ipTCM)%A(iASM+nASM*(iISM-1))*SQ2*SQ3
          end if

          AmpL1(iA+nAshA,iB) = AmpL1(iA+nAshA,iB)+ValEP+ValEM
          AmpL1(iB,iA+nAshB) = AmpL1(iB,iA+nAshB)+ValEP-ValEM
        end do
      end do

      AmpL1(:,:) = AmpL1(:,:)*Fac

      !! Calculate the actual contributions
      !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
      call OLagNS_post1(3,3,ERI1,AmpL1)

      !! Prepare for implicit (VV|VO) integrals
      !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
      !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu}
      call OLagNS_post2(nAshA+nSshA,nAshB+nSshB,nCorA,nCorB,AmpL1,WRK2)
      !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
      call OLagNS_post3(iIabs,iJabs,T2AO,WRK2)
    end do
  end do

end subroutine OLagNS_E

subroutine OLagNS_F(AmpL1)

  real(kind=wp), intent(out) :: AmpL1(nSshA,nSshB)
  integer(kind=iwp) :: iA, iB, iI, iJ

  if (nAshI*nAshJ*nSshA*nSshB == 0) return

  !nJ = nIshJ
  do iI=1,nAshI
    iIabs = iI+nIshI+nAes(iSymI)
    !iItot = iI+nCorI
    !if (iSymI == iSymJ) nJ = iI
    do iJ=1,nAshI
      iJabs = iJ+nIshJ+nIes(iSymJ)
      !iJtot = iJ+nCorJ
      if (iIabs < iJabs) cycle
      Fac = One
      !if ((iI /= iJ) .and. (iSymI == iSymJ)) Fac = Two
      if (iSymI /= iSymJ) Fac = Two

      call Exch(iSymA,iSymI,iSymB,iSymJ,iI+nCorI,iJ+nCorJ,ERI1,Scr)

      AmpL1(:,:) = Zero

      iTabs = iI+nAes(iSymI)
      iUabs = iJ+nAes(iSymJ)
      do iA=1,nSshA
        iAabs = iA+nSes(iSymA)
        !iAtot = iA+nIsh(iSymA)+nAsh(iSymA)
        do iB=1,nSshB
          iBabs = iB+nSes(iSymB)
          if (iAabs < iBabs) cycle
          !iBtot = iB+nIsh(iSymB)+nAsh(iSymB)

          iASP = kTgeU(iTabs,iUabs)-nTgeUes(iSym)
          iISP = kAgeB(iAabs,iBabs)-nAgeBes(iSym)
          ValFP = GA_Arrays(ipTCP)%A(iASP+nASP*(iISP-1))
          if (iIabs == iJabs) ValFP = ValFP*Half
          ValFM = Zero
          if (iAabs /= iBabs) then
            if (iTabs /= iUabs) then
              iASM = kTgtU(iTabs,iUabs)-nTgtUes(iSym)
              iISM = kAgtB(iAabs,iBabs)-nAgtBes(iSym)
              ValFM = GA_Arrays(ipTCM)%A(iASM+nASM*(iISM-1))
            end if
          else
            ValFP = ValFP*SQI2
          end if
          VALFM = -VALFM !! why?

          AmpL1(iA,iB) = AmpL1(iA,iB)+ValFP+ValFM
          AmpL1(iB,iA) = AmpL1(iB,iA)+ValFP-ValFM
        end do
      end do

      AmpL1(:,:) = AmpL1(:,:)*Fac

      !! Calculate the actual contributions
      !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
      call OLagNS_post1(2,2,ERI1,AmpL1)

      !! Prepare for implicit (VV|VO) integrals
      !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
      !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu}
      call OLagNS_post2(nSshA,nSshB,nOccA,nOccB,AmpL1,WRK2)
      !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
      call OLagNS_post3(iIabs,iJabs,T2AO,WRK2)
    end do
  end do

end subroutine OLagNS_F

subroutine OLagNS_G(AmpL1)

  real(kind=wp), intent(out) :: AmpL1(nSshA,nSshB)
  integer(kind=iwp) :: iA, iB, iI, iJ

  if (nAshI*nIshJ*nSshA*nSshB == 0) return

  !nJ = nIshJ
  do iI=1,nAshI
    iIabs = iI+nIshI+nAes(iSymI)
    !iItot = iI + nCorI
    !if (iSymI == iSymJ) nJ = iI
    do iJ=1,nIshI
      iJabs = iJ+nIes(iSymJ)
      !iJtot = iJ
      !if (iIabs < iJabs) cycle
      Fac = One
      !if ((iI /= iJ) .and. (iSymI == iSymJ)) Fac = Two
      if (iSymI /= iSymJ) Fac = Two

      call Exch(iSymA,iSymI,iSymB,iSymJ,iI+nCorI,iJ+nFroJ,ERI1,Scr)

      AmpL1(:,:) = Zero

      do iA=1,nSshA
        iAabs = iA+nSes(iSymA)
        !iAtot = iA+nIsh(iSymA)+nAsh(iSymA)
        do iB=1,nSshB
          iBabs = iB+nSes(iSymB)
          if (iAabs < iBabs) cycle
          !iBtot = iB+nIsh(iSymB)+nAsh(iSymB)

          iAgeB = kAgeB(iAabs,iBabs)-nAgeBes(iSym) !! iSymAB
          iVjP = iJ+nIsh(iSymJ)*(iAgeB-1)+IOFF1(iSymJ)
          ValGP = GA_Arrays(ipTCP)%A(iI+nASP*(iVjP-1))
          ValGM = Zero
          if (iAabs /= iBabs) then
            ValGP = ValGP*SQ2
            iAgtB = kAgtB(iAabs,iBabs)-nAgtBes(iSym) !! iSymAB
            iVjM = iJ+nIsh(iSymJ)*(iAgtB-1)+IOFF2(iSymJ)
            ValGM = GA_Arrays(ipTCM)%A(iI+nASM*(iVjM-1))*SQ2*SQ3
          end if

          AmpL1(iA,iB) = AmpL1(iA,iB)+ValGP+ValGM
          AmpL1(iB,iA) = AmpL1(iB,iA)+ValGP-ValGM
        end do
      end do

      AmpL1(:,:) = AmpL1(:,:)*Fac

      !! Calculate the actual contributions
      !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
      call OLagNS_post1(2,2,ERI1,AmpL1)
      !call DGEMM_('N','T',nOrbA,nSshA,nSshB,One,ERI1(1+nOrbA*nOccB),nOrbA,AmpL1,nSshA,One,OLAG(nOrbA*nOccB+1),nOrbA)
      !call DGEMM_('T','N',nOrbA,nSshA,nSshB,One,ERI1(nOccA+1),nOrbA,AmpL1,nSshA,One,OLAG(nOrbA*nOccB+1),nOrbA)

      !! Prepare for implicit (VV|VO) integrals
      !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
      !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu} for all ij
      call OLagNS_post2(nSshA,nSshB,nOccA,nOccB,AmpL1,WRK2)
      !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
      call OLagNS_post3(iIabs,iJabs,T2AO,WRK2)
    end do
  end do

end subroutine OLagNS_G

subroutine OLagNS_H(AmpL1)

  real(kind=wp), intent(out) :: AmpL1(nSshA,nSshB)
  integer(kind=iwp) :: iA, iB, iI, iJ

  if (nIshI*nIshJ*nSshA*nSshB == 0) return

  nJ = nIshJ
  do iI=1,nIshI
    iIabs = iI+nIes(iSymI)
    !iItot = iI
    if (iSymI == iSymJ) nJ = iI
    do iJ=1,nJ
      iJabs = iJ+nIes(iSymJ)
      !iJtot = iJ
      if (iIabs < iJabs) cycle
      Fac = One
      !if ((iI /= iJ) .and. (iSymI == iSymJ)) Fac = Two
      if (iSymI /= iSymJ) Fac = Two

      call Exch(iSymA,iSymI,iSymB,iSymJ,iI+nFroI,iJ+nFroJ,ERI1,Scr)

      AmpL1(:,:) = Zero

      iViHP0 = kIgeJ(iIabs,iJabs)-nIgeJes(iSym)
      iViHP = nAgeB(iSym)*(iViHP0-1)
      iViHM0 = kIgtJ(iIabs,iJabs)-nIgtJes(iSym)
      iViHM = nAgtB(iSym)*(iViHM0-1)
      do iA=1,nSshA
        iAabs = iA+nSes(iSymA)
        !iAtot = iA+nIsh(iSymA)+nAsh(iSymA)
        do iB=1,nSshB
          iBabs = iB+nSes(iSymB)
          if (iAabs < iBabs) cycle
          !iBtot = iB+nIsh(iSymB)+nAsh(iSymB)
          iVaHP = kAgeB(iAabs,iBabs)-nAgeBes(iSym)
          iVHP = iVaHP+iViHP !! nAgeB(iSym)*(iViP-1)

          ValHP = GA_Arrays(ipTCP)%A(iVHP)
          ValHM = Zero
          if (iIabs /= iJabs) then
            if (iAabs /= iBabs) then
              ValHP = ValHP*Two
              iVaHM = kAgtB(iAabs,iBabs)-nAgtBes(iSym)
              iVHM = iVaHM+iViHM !! nAgtB(iSym)*(iViM-1)
              ValHM = GA_Arrays(ipTCM)%A(iVHM)*Two*SQ3
            else
              ValHP = ValHP*SQ2
            end if
          else
            if (iAabs /= iBabs) ValHP = ValHP*SQ2
          end if

          !write(u6,'(2i3,2f20.10)') ia,ib,valhp,valhm
          AmpL1(iA,iB) = AmpL1(iA,iB)+ValHP+ValHM
          AmpL1(iB,iA) = AmpL1(iB,iA)+ValHP-ValHM
        end do
      end do

      AmpL1(:,:) = AmpL1(:,:)*Fac

      !! Calculate the actual contributions
      !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
      call OLagNS_post1(2,2,ERI1,AmpL1)
      !call DGEMM_('N','T',nOrbA,nSshA,nSshB,One,ERI1(1+nOrbA*nOccB),nOrbA,AmpL1,nSshA,One,OLAG(nOrbA*nOccB+1),nOrbA)
      !call DGEMM_('T','N',nOrbA,nSshA,nSshB,One,ERI1(nOccA+1),nOrbA,AmpL1,nSshA,One,OLAG(nOrbA*nOccB+1),nOrbA)

      !! Prepare for implicit (VV|VO) integrals
      !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
      !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu} for all ij
      call OLagNS_post2(nSshA,nSshB,nOccA,nOccB,AmpL1,WRK2)
      !call DGEMM_('N','N',nBasT,nSshA,nSshB,One,CMOPT2(1+nBasT*nOccB),nBasT,AmpL1,nSshA,Zero,WRK1,nBasT)
      !call DGEMM_('N','T',nBasT,nBasT,nSshB,One,WRK1,nBasT,CMOPT2(1+nBasT*nOccA),nBasT,Zero,WRK2,nBasT)
      !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
      call OLagNS_post3(iIabs,iJabs,T2AO,WRK2)
      !do iBas=1,nBasT
      !  do jBas=1,nBasT
      !    loc1 = iJ-1+(jBas-1)*nOccA2+(iI-1)*nOccA2*nBasT+(iBas-1)*nOccA2*nBasT*nOccA2
      !    loc2 = iI-1+(jBas-1)*nOccA2+(iJ-1)*nOccA2*nBasT+(iBas-1)*nOccA2*nBasT*nOccA2
      !    loc3 = iBas-1+(jBas-1)*nBasT
      !    loc4 = jBas-1+(iBas-1)*nBasT
      !    T2AO(1+loc1) = T2AO(1+loc1)+WRK2(1+loc3)
      !    T2AO(1+loc2) = T2AO(1+loc2)+WRK2(1+loc4)
      !  end do
      !end do
    end do
  end do

end subroutine OLagNS_H

subroutine OLagNS_post1(iLeft,iRight,ERI,AmpMO)

  integer(kind=iwp), intent(in) :: iLeft, iRight
  real(kind=wp), intent(in) :: ERI(NBSQT), AmpMO(NBSQT)
  integer(kind=iwp) :: nDimA, nDimB, nSkpA, nSkpB

  nSkpA = 0
  nSkpB = 0
  nDimA = 0
  nDimB = 0
  if (iLeft == 1) then
    nDimA = nAshA
    nSkpA = nCorA
  else if (iLeft == 2) then
    nDimA = nSshA
    nSkpA = nOccA
  else if (iLeft == 3) then
    nDimA = nAshA+nSshA
    nSkpA = nCorA
  end if
  if (iRight == 1) then
    nDimB = nAshB
    nSkpB = nCorB
  else if (iRight == 2) then
    nDimB = nSshB
    nSkpB = nOccB
  else if (iRight == 3) then
    nDimB = nAshB+nSshB
    nSkpB = nCorB
  end if

  call DGEMM_('N','T',nOrbA,nDimA,nDimB,One,ERI(1+nOrbA*nSkpB),nOrbA,AmpMO,nDimA,One,OLAG(nOrbA*nSkpB+1),nOrbA)
  call DGEMM_('T','N',nOrbA,nDimA,nDimB,One,ERI(nSkpA+1),nOrbA,AmpMO,nDimA,One,OLAG(nOrbA*nSkpB+1),nOrbA)

end subroutine OLagNS_post1

subroutine OLagNS_post2(nDimA,nDimB,nSkpA,nSkpB,AmpMO,AmpAO)

  use caspt2_global, only: CMOPT2

  integer(kind=iwp), intent(in) :: nDimA, nDimB, nSkpA, nSkpB
  real(kind=wp), intent(in) :: AmpMO(nDimA,nDimB)
  real(kind=wp), intent(out) :: AmpAO(nBasA,nBasB)

  call DGEMM_('N','N',nBasA,nDimB,nDimA,One,CMOPT2(1+nBasA*nSkpA),nBasA,AmpMO,nDimA,Zero,WRK1,nBasA)
  call DGEMM_('N','T',nBasA,nBasB,nDimB,One,WRK1,nBasA,CMOPT2(1+nBasB*nSkpB),nBasA,Zero,AmpAO,nBasA)

end subroutine OLagNS_post2

subroutine OLagNS_post3(iIabs,iJabs,TampAO,TampIJ)

  integer(kind=iwp), intent(in) :: iIabs, iJabs
  real(kind=wp), intent(inout) :: TampAO(nOccA2,nBasI,nOccB2,nBasJ)
  real(kind=wp), intent(in) :: TampIJ(nBasI,nBasJ)
  integer(kind=iwp) :: iBas

  do iBas=1,nBasI
    ! FIXME: this can't work if iBas /= jBas
    TampAO(iJabs,:,iIabs,iBas) = TampAO(iJabs,:,iIabs,iBas)+TampIJ(iBas,:)
    TampAO(iIabs,:,iJabs,iBas) = TampAO(iIabs,:,iJabs,iBas)+TampIJ(:,iBas)
  end do

end subroutine OLagNS_post3

subroutine OLagNS_AAAA(AmpL1)

  use BDerNEV, only: Gder

  real(kind=wp), intent(inout) :: AmpL1(nAshA,nAshB)
  integer(kind=iwp) :: iA, iB, iI, iJ

  if (nAshI*nAshJ*nAshA*nAshB == 0) return

  do iI=1,nAshI
    iIabs = iI+nIshI+nAes(iSymI)
    do iJ=1,nAshJ
      iJabs = iJ+nIshJ+nAes(iSymJ)
      !Fac = One

      call Exch(iSymA,iSymI,iSymB,iSymJ,iI+nCorI,iJ+nCorJ,ERI1,Scr)
      AmpL1(:,:) = Zero

      do iA=1,nAshA
        iAabs = iA+nAes(iSymA)
        do iB=1,nAshB
          iBabs = iB+nAes(iSymB)
          ValA = Gder(iAabs,iI,iBabs,iJ)
          AmpL1(iA,iB) = AmpL1(iA,iB)+ValA
        end do
      end do

      !AmpL1(:,:) = Fac*AmpL1(:,:)
      call OLagNS_post1(1,1,ERI1,AmpL1)
      call OLagNS_post2(nAshA,nAshB,nCorA,nCorB,AmpL1,WRK2)
      call OLagNS_post3(iIabs,iJabs,T2AO,WRK2)
    end do
  end do

end subroutine OLagNS_AAAA

end subroutine OLagNS_Hel2
