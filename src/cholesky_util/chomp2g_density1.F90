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
! Copyright (C) 2010, Jonas Bostrom                                    *
!***********************************************************************

subroutine ChoMP2g_density1(irc,EOcc,EVir,EFro,Wrk,lWrk)
!
! Jonas Bostrom, Feb 2010
!
! Purpose: To Compute MP2 density from Cholesky MO-vectors and
!          decomposed MP2 amplitudes.

use Symmetry_Info, only: Mul
use Cholesky, only: nSym, NumCho
use ChoMP2, only: AdrR1, AdrR2, iAdrOff, iFro, iOcc, iT1am, kFLagr, kLagr, kPab, kPai, kPaK, kPij, kPik, kWab, kWai, kWaK, kWij, &
                  kWiK, kWJK, lFLagr, lLagr, lPab, lPai, lPaK, lPij, lPiK, lUnit_F, LuRInv, LuUVec, LuVVec, LuWVec, lWab, lWai, &
                  lWaK, lWij, lWiK, lWJK, nFro, nMoMo, nMP2Vec, nOcc, nVir
use Constants, only: Zero, One, Two, Four, Eight, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(in) :: EOcc(*), EVir(*), EFro(*)
integer(kind=iwp), intent(in) :: lWrk
real(kind=wp), intent(out) :: Wrk(lWrk)
integer(kind=iwp) :: i, iAdr, iB, iB1, iBBlock, iBRel, iClos, iI, iJ, iK, iOff, iOff1, iOffAmp, iOffL, iOffL1, iOffLic, iOffRic, &
                     iOffU, iOffU1, iOffV, iOffV1, iOffX, iOpt, iSeed, iSym, iSymA, iSymB, iSymC, iSymI, iSymJ, iSymJC, iSymK, &
                     iTypL, iTypR, iVecFV, iVecOF, iVecOO, iVecOV, iVecVV, jBat, jSym, jVec, jVec1, kAmp, kBat, kEndAmp, &
                     kEndFLagr, kEndLab, kEndLagr, kEndLic, kEndLij, kEndliK, kEndLji, kEndLKa, kEndPab, kEndPai, kEndPaK, &
                     kEndPij, kEndPiK, kEndRia, kEndRia2, kEndRib, kEndRic, kEndRjc, kEndU, kEndWab, kEndWai, kEndWaK, kEndWij, &
                     kEndWiK, kEndWJK, kEndX, kLab, kLia, kLic, kLij, kLiK, kLji, kLKa, kRia, kRia2, kRib, kRic, kRjc, kU, kV, &
                     kVec, kVec1, kX, lAmp, lLab, lLic, lLij, lLiK, lLji, lLKa, lRia, lRia2, lRib, lRic, lRjc, lScr, lTot, lU, lV, &
                     lX, maxvalue, nB, nBatL, nBatMax, nBatR, nBBlock, nFroMax, nMoMoMax, nMP2VecMax, nOccMax, NumB, NumChoMax, &
                     NumIC, NumLVecJ, NumRVecJ, NumVecJ, NumVecK, nVec, nVirMax
real(kind=wp) :: E_i, E_K, Fac
character(len=5) :: Fname
real(kind=wp), parameter :: X(0:1) = [Zero,One]
character(len=*), parameter :: SecNam = 'ChoMP2g_density1'
integer(kind=iwp), external :: IsFreeUnit

#include "macros.fh"
unused_var(EVir(1))

maxvalue = 200

! Do not delete vectors
! ---------------------
iClos = 2

! Set type of Choleskyvectors
! ---------------------------
iTypL = 1
iTypR = 2
iVecOV = 6
iVecOF = 4
iVecOO = 5
iVecFV = 3
iVecVV = 9

! Calc max number of cholesky vectors in a specific sym
! -----------------------------------------------------
nMP2VecMax = 0
NumChoMax = 0
do i=1,nSym
  nMP2VecMax = max(nMP2VecMax,nMP2Vec(i))
  NumChoMax = max(NumChoMax,NumCho(i))
end do

nMoMoMax = 0
do iSym=1,nSym
  nMoMoMax = max(nMoMo(iSym,iVecOV),nMoMoMax)
end do

! Allocate Memory for Pab-Vector
! ------------------------------
lPab = nVir(1)*nVir(1)
kPab(1) = 1
do iSym=2,nSym
  kPab(iSym) = kPab(iSym-1)+nVir(iSym-1)*nVir(iSym-1)
  lPab = lPab+nVir(iSym)*nVir(iSym)
end do
kEndPab = kPab(1)+lPab
Wrk(kPab(1):kEndPab-1) = Zero

! Allocate Memory for Wab-Vector
! ------------------------------
lWab = nVir(1)*nVir(1)
kWab(1) = kEndPab
do iSym=2,nSym
  kWab(iSym) = kWab(iSym-1)+nVir(iSym-1)*nVir(iSym-1)
  lWab = lWab+nVir(iSym)*nVir(iSym)
end do
kEndWab = kWab(1)+lWab
Wrk(kWab(1):kEndWab-1) = Zero

! Allocate memory for Pij-Vector
! ------------------------------

lPij = nOcc(1)*nOcc(1)
kPij(1) = kEndWab
do iSym=2,nSym
  kPij(iSym) = kPij(iSym-1)+nOcc(iSym-1)*nOcc(iSym-1)
  lPij = lPij+nOcc(iSym)*nOcc(iSym)
end do
kEndPij = kPij(1)+lPij
Wrk(kPij(1):kEndPij-1) = Zero

! Allocate memory for Wij-Vector
! ------------------------------

lWij = nOcc(1)*nOcc(1)
kWij(1) = kEndPij
do iSym=2,nSym
  kWij(iSym) = kWij(iSym-1)+nOcc(iSym-1)*nOcc(iSym-1)
  lWij = lWij+nOcc(iSym)*nOcc(iSym)
end do
kEndWij = kWij(1)+lWij
Wrk(kWij(1):kEndWij-1) = Zero

! Allocate memory for Pai-Vector
! ------------------------------

lPai = nVir(1)*nOcc(1)
kPai(1) = kEndWij
do iSym=2,nSym
  kPai(iSym) = kPai(iSym-1)+nVir(iSym-1)*nOcc(iSym-1)
  lPai = lPai+nVir(iSym)*nOcc(iSym)
end do
kEndPai = kPai(1)+lPai
Wrk(kPai(1):kEndPai-1) = Zero

! Allocate memory for Wai-Vector
! ------------------------------

lWai = nVir(1)*nOcc(1)
kWai(1) = kEndPai
do iSym=2,nSym
  kWai(iSym) = kWai(iSym-1)+nVir(iSym-1)*nOcc(iSym-1)
  lWai = lWai+nVir(iSym)*nOcc(iSym)
end do
kEndWai = kWai(1)+lWai
Wrk(kWai(1):kEndWai-1) = Zero

! Allocate memory for PaK-Vector
! ------------------------------

lPaK = nVir(1)*nFro(1)
kPaK(1) = kEndWai
do iSym=2,nSym
  kPaK(iSym) = kPaK(iSym-1)+nVir(iSym-1)*nFro(iSym-1)
  lPaK = lPaK+nVir(iSym)*nFro(iSym)
end do
kEndPaK = kPaK(1)+lPaK
Wrk(kPaK(1):kEndPaK-1) = Zero

! Allocate memory for WaK-Vector
! ------------------------------

lWaK = nVir(1)*nFro(1)
kWaK(1) = kEndPaK
do iSym=2,nSym
  kWaK(iSym) = kWaK(iSym-1)+nVir(iSym-1)*nOcc(iSym-1)
  lWaK = lWaK+nVir(iSym)*nOcc(iSym)
end do
kEndWaK = kWaK(1)+lWaK
Wrk(kWaK(1):kEndWaK-1) = Zero

! Allocate memory for PiK-vector (occ-fro)
! ----------------------------------------

lPiK = nOcc(1)*nFro(1)
kPiK(1) = kEndWaK
do iSym=2,nSym
  kPiK(iSym) = kPiK(iSym-1)+nOcc(iSym-1)*nFro(iSym-1)
  lPiK = lPiK+nOcc(iSym)*nFro(iSym)
end do
kEndPiK = kPiK(1)+lPiK
Wrk(kPiK(1):kEndPiK-1) = Zero

! Allocate memory for WiK-vector (occ-fro)
! ----------------------------------------

lWiK = nOcc(1)*nFro(1)
kWiK(1) = kEndPiK
do iSym=2,nSym
  kWiK(iSym) = kWiK(iSym-1)+nOcc(iSym-1)*nFro(iSym-1)
  lWiK = lWiK+nOcc(iSym)*nFro(iSym)
end do
kEndWiK = kWiK(1)+lWiK
Wrk(kWiK(1):kEndWiK-1) = Zero

! Allocate memory for WJK-vector (fro-fro)
! ----------------------------------------

lWJK = nFro(1)*nFro(1)
kWJK(1) = kEndWiK
do iSym=2,nSym
  kWJK(iSym) = kWJK(iSym-1)+nFro(iSym-1)*nFro(iSym-1)
  lWJK = lWJK+nFro(iSym)*nFro(iSym)
end do
kEndWJK = kWJK(1)+lWJK
Wrk(kWJK(1):kEndWJK-1) = Zero

! Allocate memory for Lagr-vector
! ------------------------------

lLagr = nOcc(1)*nVir(1)
kLagr(1) = kEndWJK
do iSym=2,nSym
  kLagr(iSym) = kLagr(iSym-1)+nOcc(iSym-1)*nVir(iSym-1)
  lLagr = lLagr+nOcc(iSym)*nVir(iSym)
end do
kEndLagr = kLagr(1)+lLagr
Wrk(kLagr(1):kEndLagr-1) = Zero

! Allocate memory for FrozenLagr-vector
! -------------------------------------

lFLagr = nFro(1)*nVir(1)
kFLagr(1) = kEndLagr
do iSym=2,nSym
  kFLagr(iSym) = kFLagr(iSym-1)+nFro(iSym-1)*nVir(iSym-1)
  lFLagr = lFLagr+nFro(iSym)*nVir(iSym)
end do
kEndFLagr = kFLagr(1)+lFLagr
Wrk(kFlagr(1):kEndFlagr-1) = Zero

nVec = min(maxvalue,max(nMP2VecMax,NumChoMax))
if (nVec < 1) call SysAbendMsg(SecNam,'Insufficient memory','[1]')

! Allocate memory for X^KJ-vector
! -------------------------------

lX = nVec*nMp2VecMax
kX = kEndFLagr
kEndX = kX+lX

do iSym=1,nSym

  ! Allocate memory for Ria-vectors
  ! -------------------------------
  lRia = nMoMo(iSym,iVecOV)*nVec
  kRia = kEndX
  kEndRia = kRia+lRia
  kLia = kRia

  lRia2 = nMoMo(iSym,iVecOV)*nVec
  kRia2 = kEndRia
  kEndRia2 = kRia2+lRia2

  ! Allocate memory for L-vectors
  ! -----------------------------
  lLKa = nMoMo(iSym,iVecFV)*nVec
  kLKa = kEndRia2
  kEndLKa = kLKa+lLKa

  lLab = nMoMo(iSym,iVecVV)*nVec
  kLab = kEndLKa
  kEndLab = kLab+lLab

  lLij = nMoMo(iSym,iVecOO)*nVec
  kLij = kEndLab
  kEndLij = kLij+lLij

  lLiK = nMoMo(iSym,iVecOF)*nVec
  kLiK = kEndLij
  kEndLiK = kLiK+lLiK

  ! Allocate memory for U-vector
  ! ----------------------------

  lU = nMoMo(iSym,iVecOV)*nVec
  kU = kEndLiK
  kEndU = kU+lU

  ! Setup batch over amplitude vectors.
  ! -----------------------------------
  nBatR = (nMP2Vec(iSym)-1)/nVec+1
  nBatL = (NumCho(iSym)-1)/nVec+1
  if ((nMoMo(iSym,iVecOV) > 0) .and. (nMp2Vec(iSym) > 0)) then

    ! Open the File for reordered amplitude vectors
    ! ---------------------------------------------
    call ChoMP2_OpenF(1,iTypR,iSym)
    call ChoMP2_OpenF(1,iTypL,iSym)

    iSeed = 7
    LuUVec = IsFreeUnit(iSeed)
    write(Fname,'(A4,I1)') 'TMPV',1
    call DaName_MF_WA(LuUVec,Fname)

    iSeed = 7
    LuVVec = IsFreeUnit(iSeed)
    write(Fname,'(A4,I1)') 'TMPV',2
    call DaName_MF_WA(LuVVec,Fname)

    ! Calculate Intermediate vectors U
    ! --------------------------------
    do kBat=1,nBatR
      Wrk(kX:kX+lX-1) = Zero
      if (kBat == nBatR) then
        NumVecK = nMP2Vec(iSym)-nVec*(nBatR-1)
      else
        NumVecK = nVec
      end if
      kVec = nVec*(kBat-1)+1

      ! Read Amplitude vectors from kBat
      ! --------------------------------
      iOpt = 2
      lTot = nMoMo(iSym,iVecOV)*NumVecK
      iAdr = nMoMo(iSym,iVecOV)*(kVec-1)+1
      call dDaFile(lUnit_F(iSym,iTypR),iOpt,Wrk(kRia),lTot,iAdr)

      do jBat=1,nBatR
        if (jBat == nBatR) then
          NumVecJ = nMP2Vec(iSym)-nVec*(nBatR-1)
        else
          NumVecJ = nVec
        end if
        jVec = nVec*(jBat-1)+1

        ! Read Amplitude vectors from jBat
        ! --------------------------------
        iOpt = 2
        lTot = nMoMo(iSym,iVecOV)*NumVecJ
        iAdr = nMoMo(iSym,iVecOV)*(jVec-1)+1
        call dDaFile(lUnit_F(iSym,iTypR),iOpt,Wrk(kRia2),lTot,iAdr)

        ! Construct X^JK-vector
        ! ---------------------
        iOffX = NumVecK*(jVec-1)
        call dGemm_('T','N',NumVecK,NumVecJ,nMoMo(iSym,iVecOV),One,Wrk(kRia),nMoMo(iSym,iVecOV),Wrk(kRia2),nMoMo(iSym,iVecOV), &
                    Zero,Wrk(kX+iOffX),NumVecK)

      end do

      iOpt = 1
      lTot = nMP2Vec(iSym)*NumVecK
      iAdr = 1+nMP2Vec(iSym)*(kVec-1)
      call dDaFile(LuVVec,iOpt,Wrk(kX),lTot,iAdr)
    end do
    do kBat=1,nBatR
      if (kBat == nBatR) then
        NumVecK = nMP2Vec(iSym)-nVec*(nBatR-1)
      else
        NumVecK = nVec
      end if
      kVec = nVec*(kBat-1)+1

      iOpt = 2
      lTot = nMP2Vec(iSym)*NumVecK
      iAdr = 1+nMP2Vec(iSym)*(kVec-1)
      call dDaFile(LuVVec,iOpt,Wrk(kX),lTot,iAdr)

      do jBat=1,nBatR
        if (jBat == nBatR) then
          NumVecJ = nMP2Vec(iSym)-nVec*(nBatR-1)
        else
          NumVecJ = nVec
        end if
        jVec = nVec*(jBat-1)+1

        ! Read Amplitude vectors from jBat
        ! --------------------------------
        iOpt = 2
        lTot = nMoMo(iSym,iVecOV)*NumVecJ
        iAdr = nMoMo(iSym,iVecOV)*(jVec-1)+1
        call dDaFile(lUnit_F(iSym,iTypR),iOpt,Wrk(kRia),lTot,iAdr)

        iOffX = NumVecK*(jVec-1)
        Fac = X(min((jBat-1),1))
        call dGemm_('N','T',nMoMo(iSym,iVecOV),NumVecK,NumVecJ,One,Wrk(kRia),nMoMo(iSym,iVecOV),Wrk(kX+iOffX),NumVecK,Fac, &
                    Wrk(kU),nMoMo(iSym,iVecOV))
      end do

      iOpt = 1
      lTot = nMoMo(iSym,iVecOV)*NumVecK
      iAdr = 1+nMoMo(iSym,iVecOV)*(kVec-1)
      call dDaFile(LuUVec,iOpt,Wrk(kU),lTot,iAdr)
    end do

    ! Calculate "Coulomb"-contributions to Densities from U
    ! -----------------------------------------------------
    do kBat=1,nBatR
      if (kBat == nBatR) then
        NumVecK = nMP2Vec(iSym)-nVec*(nBatR-1)
      else
        NumVecK = nVec
      end if
      kVec = nVec*(kBat-1)+1

      ! Read Amplitude vectors from kBat
      ! --------------------------------
      iOpt = 2
      lTot = nMoMo(iSym,iVecOV)*NumVecK
      iAdr = nMoMo(iSym,iVecOV)*(kVec-1)+1
      call dDaFile(lUnit_F(iSym,iTypR),iOpt,Wrk(kRia),lTot,iAdr)

      ! Load intermediate vectors U^K_ib
      ! --------------------------------
      iOpt = 2
      lTot = nMoMo(iSym,iVecOV)*NumVecK
      iAdr = nMoMo(iSym,iVecOv)*(kVec-1)+1
      call dDaFile(LuUVec,iOpt,Wrk(kU),lTot,iAdr)

      ! Calculate the "Coulomb" Contribution to Pab
      ! -------------------------------------------
      iOff1 = 0
      do iSymI=1,nSym
        iSymA = Mul(iSym,iSymI)
        if (nOcc(iSymI)*nVir(iSymA) /= 0) then
          do iI=1,nOcc(iSymI)
            iOff = (iI-1)*nVir(iSymA)+iOff1
            call dGemm_('N','T',nVir(iSymA),nVir(iSymA),NumVecK,Four,Wrk(kRia+iOff),nMoMo(iSym,iVecOV),Wrk(kU+iOff), &
                        nMoMo(iSym,iVecOV),One,Wrk(kPab(iSymA)),nVir(iSymA))
          end do
        end if
        iOff1 = iOff1+nOcc(iSymI)*nVir(iSymA)
      end do
      ! Calculate the "Coulomb" Contribution to Pik
      ! -------------------------------------------

      do kVec1=1,NumVecK
        iOff1 = 0
        do iSymK=1,nSym
          iSymI = iSymK
          iSymB = Mul(iSymK,iSym)
          if (nOcc(iSymK)*nVir(iSymB) /= 0) then
            iOff = nMoMo(iSym,iVecOV)*(kVec1-1)+iOff1
            call dGemm_('T','N',nOcc(iSymI),nOcc(iSymK),nVir(iSymB),-Four,Wrk(kU+iOff),nVir(iSymB),Wrk(kRia+iOff),nVir(iSymB), &
                        One,Wrk(kPij(iSymK)),nOcc(iSymI))
          end if
          iOff1 = iOff1+nOcc(iSymI)*nVir(iSymB)
        end do
      end do

    end do

    ! Calculate Intermediate vectors U*
    ! ---------------------------------
    do kBat=1,nBatL
      Wrk(kX:kX+lX-1) = Zero
      if (kBat == nBatL) then
        NumVecK = NumCho(iSym)-nVec*(nBatL-1)
      else
        NumVecK = nVec
      end if
      kVec = nVec*(kBat-1)+1

      ! Read Integral vectors from kBat
      ! --------------------------------
      iOpt = 2
      lTot = nMoMo(iSym,iVecOV)*NumVecK
      iAdr = 1+nMoMo(iSym,iVecOV)*(kVec-1)+iAdrOff(iSym,iVecOV)
      call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLia),lTot,iAdr)

      do jBat=1,nBatR
        if (jBat == nBatR) then
          NumVecJ = nMP2Vec(iSym)-nVec*(nBatR-1)
        else
          NumVecJ = nVec
        end if
        jVec = nVec*(jBat-1)+1

        ! Read Amplitude vectors from jBat
        ! --------------------------------
        iOpt = 2
        lTot = nMoMo(iSym,iVecOV)*NumVecJ
        iAdr = nMoMo(iSym,iVecOV)*(jVec-1)+1
        call dDaFile(lUnit_F(iSym,iTypR),iOpt,Wrk(kRia2),lTot,iAdr)

        ! Construct X^JK-vector
        ! ---------------------
        iOffX = NumVecK*(jVec-1)
        call dGemm_('T','N',NumVecK,NumVecJ,nMoMo(iSym,iVecOV),One,Wrk(kLia),nMoMo(iSym,iVecOV),Wrk(kRia2),nMoMo(iSym,iVecOV), &
                    Zero,Wrk(kX+iOffX),NumVecK)
      end do

      iOpt = 1
      lTot = nMP2Vec(iSym)*NumVecK
      iAdr = 1+nMP2Vec(iSym)*(kVec-1)
      call dDaFile(LuVVec,iOpt,Wrk(kX),lTot,iAdr)
    end do
    do kBat=1,nBatL
      if (kBat == nBatL) then
        NumVecK = NumCho(iSym)-nVec*(nBatL-1)
      else
        NumVecK = nVec
      end if
      kVec = nVec*(kBat-1)+1

      iOpt = 2
      lTot = nMP2Vec(iSym)*NumVecK
      iAdr = 1+nMP2Vec(iSym)*(kVec-1)
      call dDaFile(LuVVec,iOpt,Wrk(kX),lTot,iAdr)

      do jBat=1,nBatR
        if (jBat == nBatR) then
          NumVecJ = nMP2Vec(iSym)-nVec*(nBatR-1)
        else
          NumVecJ = nVec
        end if
        jVec = nVec*(jBat-1)+1

        ! Read Amplitude vectors from jBat
        ! --------------------------------
        iOpt = 2
        lTot = nMoMo(iSym,iVecOV)*NumVecJ
        iAdr = nMoMo(iSym,iVecOV)*(jVec-1)+1
        call dDaFile(lUnit_F(iSym,iTypR),iOpt,Wrk(kRia),lTot,iAdr)

        iOffX = NumVecK*(jVec-1)
        Fac = X(min((jBat-1),1))
        call dGemm_('N','T',nMoMo(iSym,iVecOV),NumVecK,NumVecJ,One,Wrk(kRia),nMoMo(iSym,iVecOV),Wrk(kX+iOffX),NumVecK,Fac,Wrk(kU), &
                    nMoMo(iSym,iVecOV))
      end do

      iOpt = 1
      lTot = nMoMo(iSym,iVecOV)*NumVecK
      iAdr = 1+nMoMo(iSym,iVecOV)*(kVec-1)
      call dDaFile(LuUVec,iOpt,Wrk(kU),lTot,iAdr)
    end do

    ! Calculate "Coulomb"-contributions to Densities and
    ! the Lagrangian from U*
    ! -----------------------------------------------------
    do kBat=1,nBatL

      if (kBat == nBatL) then
        NumVecK = NumCho(iSym)-nVec*(nBatL-1)
      else
        NumVecK = nVec
      end if
      kVec = nVec*(kBat-1)+1

      ! Read Integral vectors(FV) from kBat
      ! --------------------------------
      iOpt = 2
      lTot = nMoMo(iSym,iVecFV)*NumVecK
      iAdr = 1+nMoMo(iSym,iVecFV)*(kVec-1)+iAdrOff(iSym,iVecFV)
      call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLKa),lTot,iAdr)

      ! Read Integral vectors(VV) from kBat
      ! --------------------------------
      iOpt = 2
      lTot = nMoMo(iSym,iVecVV)*NumVecK
      iAdr = 1+nMoMo(iSym,iVecVV)*(kVec-1)+iAdrOff(iSym,iVecVV)
      call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLab),lTot,iAdr)

      ! Read Integral vectors(OO) from kBat
      ! --------------------------------
      iOpt = 2
      lTot = nMoMo(iSym,iVecOO)*NumVecK
      iAdr = 1+nMoMo(iSym,iVecOO)*(kVec-1)+iAdrOff(iSym,iVecOO)
      call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLij),lTot,iAdr)

      ! Read Integral vectors(OF) from kBat
      ! --------------------------------
      iOpt = 2
      lTot = nMoMo(iSym,iVecOF)*NumVecK
      iAdr = 1+nMoMo(iSym,iVecOF)*(kVec-1)+iAdrOff(iSym,iVecOF)
      call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLiK),lTot,iAdr)

      ! Read Integral vectors(OV) from kBat
      ! --------------------------------
      iOpt = 2
      lTot = nMoMo(iSym,iVecOV)*NumVecK
      iAdr = 1+nMoMo(iSym,iVecOV)*(kVec-1)+iAdrOff(iSym,iVecOV)
      call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLia),lTot,iAdr)

      ! Load intermediate vectors
      ! -------------------------
      iOpt = 2
      lTot = nMoMo(iSym,iVecOV)*NumVecK
      iAdr = nMoMo(iSym,iVecOv)*(kVec-1)+1
      call dDaFile(LuUVec,iOpt,Wrk(kU),lTot,iAdr)

      ! Calculate the "Coulomb" Contribution to Froz-occupied density PiK
      ! -----------------------------------------------------------------

      do kVec1=1,NumVecK
        iOffU1 = 0
        iOffL1 = 0
        do iSymK=1,nSym
          iSymI = iSymK
          iSymB = Mul(iSymK,iSym)
          if (nOcc(iSymI)*nFro(iSymK)*nVir(iSymB) /= 0) then
            iOffU = nMoMo(iSym,iVecOV)*(kVec1-1)+iOffU1
            iOffL = nMoMo(iSym,iVecFV)*(kVec1-1)+iOffL1
            call dGemm_('T','N',nOcc(iSymI),nFro(iSymK),nVir(iSymB),-Four,Wrk(kU+iOffU),nVir(iSymB),Wrk(kLKa+iOffL),nVir(iSymB), &
                        One,Wrk(kPiK(iSymK)),nOcc(iSymI))
          end if
          iOffU1 = iOffU1+nOcc(iSymI)*nVir(iSymB)
          iOffL1 = iOffL1+nFro(iSymK)*nVir(iSymB)
        end do
      end do

      ! Calculate the "Coulomb" contribution to Wij (froz-occ)
      ! ------------------------------------------------------

      do kVec1=1,NumVecK
        iOffU1 = 0
        iOffL1 = 0
        do iSymK=1,nSym
          iSymI = iSymK
          iSymB = Mul(iSymK,iSym)
          if (nOcc(iSymI)*nFro(iSymK)*nVir(iSymB) /= 0) then
            iOffU = nMoMo(iSym,iVecOV)*(kVec1-1)+iOffU1
            iOffL = nMoMo(iSym,iVecFV)*(kVec1-1)+iOffL1
            call dGemm_('T','N',nOcc(iSymI),nFro(iSymK),nVir(iSymB),Four,Wrk(kU+iOffU),nVir(iSymB),Wrk(kLKa+iOffL),nVir(iSymB), &
                        One,Wrk(kWiK(iSymK)),nOcc(iSymI))
          end if
          iOffU1 = iOffU1+nOcc(iSymI)*nVir(iSymB)
          iOffL1 = iOffL1+nFro(iSymK)*nVir(iSymB)
        end do
      end do

      ! Calculate "Coulomb"-contribution to Lagr(IV)
      ! --------------------------------------------
      do kVec1=1,NumVecK
        iOffU1 = 0
        iOffL1 = 0
        do iSymA=1,nSym
          iSymI = iSymA
          iSymB = Mul(iSymA,iSym)
          if (nOcc(iSymI)*nVir(iSymA)*nVir(iSymB) /= 0) then
            iOffU = nMoMo(iSym,iVecOV)*(kVec1-1)+iOffU1
            iOffL = nMoMo(iSym,iVecVV)*(kVec1-1)+iOffL1
            call dGemm_('T','N',nVir(iSymA),nOcc(iSymI),nVir(iSymB),Four,Wrk(kLab+iOffL),nVir(iSymB),Wrk(kU+iOffU),nVir(iSymB), &
                        One,Wrk(kLagr(iSymA)),nVir(iSymA))
          end if
          iOffU1 = iOffU1+nOcc(iSymI)*nVir(iSymB)
          iOffL1 = iOffL1+nVir(iSymA)*nVir(iSymB)
        end do
      end do

      ! Calculate "Coulomb"-contribution to Lagr(III)
      ! ---------------------------------------------
      iOffU1 = 0
      iOffL1 = 0
      do iSymJ=1,nSym
        iSymA = Mul(iSym,iSymJ)
        iSymI = iSymA
        if (nOcc(iSymI)*nVir(iSymA)*nOcc(iSymJ) /= 0) then
          do iJ=1,nOcc(iSymJ)
            iOffU = (iJ-1)*nVir(iSymA)+iOffU1
            iOffL = (iJ-1)*nOcc(iSymI)+iOffL1
            call dGemm_('N','T',nVir(iSymA),nOcc(iSymI),NumVecK,-Four,Wrk(kU+iOffU),nMoMo(iSym,iVecOV),Wrk(kLij+iOffL), &
                        nMoMo(iSym,iVecOO),One,Wrk(kLagr(iSymA)),nVir(iSymA))
          end do
        end if
        iOffU1 = iOffU1+nOcc(iSymJ)*nVir(iSymA)
        iOffL1 = iOffL1+nOcc(iSymJ)*nOcc(iSymI)
      end do

      ! Calculate "Coulomb"-contribution to Wai
      ! ---------------------------------------
      iOffU1 = 0
      iOffL1 = 0
      do iSymJ=1,nSym
        iSymA = Mul(iSym,iSymJ)
        iSymI = iSymA
        if (nOcc(iSymI)*nVir(iSymA)*nOcc(iSymJ) /= 0) then
          do iJ=1,nOcc(iSymJ)
            iOffU = (iJ-1)*nVir(iSymA)+iOffU1
            iOffL = (iJ-1)*nOcc(iSymI)+iOffL1
            call dGemm_('N','T',nVir(iSymA),nOcc(iSymI),NumVecK,Eight,Wrk(kU+iOffU),nMoMo(iSym,iVecOV),Wrk(kLij+iOffL), &
                        nMoMo(iSym,iVecOO),One,Wrk(kWai(iSymA)),nVir(iSymA))
          end do
        end if
        iOffU1 = iOffU1+nOcc(iSymJ)*nVir(iSymA)
        iOffL1 = iOffL1+nOcc(iSymJ)*nOcc(iSymI)
      end do

      ! Calculate "Coulomb"-contribution to Wab
      ! ---------------------------------------
      iOff1 = 0
      do iSymI=1,nSym
        iSymA = Mul(iSym,iSymI)
        if (nOcc(iSymI)*nVir(iSymA) /= 0) then
          do iI=1,nOcc(iSymI)
            iOff = (iI-1)*nVir(iSymA)+iOff1
            call dGemm_('N','T',nVir(iSymA),nVir(iSymA),NumVecK,Four,Wrk(kLia+iOff),nMoMo(iSym,iVecOV),Wrk(kU+iOff), &
                        nMoMo(iSym,iVecOV),One,Wrk(kWab(iSymA)),nVir(iSymA))
          end do
        end if
        iOff1 = iOff1+nOcc(iSymI)*nVir(iSymA)
      end do

      ! Calculate "Coulomb"-contribution to Wij (Wik)
      ! ---------------------------------------------
      do kVec1=1,NumVecK
        iOff1 = 0
        do iSymK=1,nSym
          iSymI = iSymK
          iSymB = Mul(iSymK,iSym)
          if (nOcc(iSymK)*nVir(iSymB) /= 0) then
            iOff = nMoMo(iSym,iVecOV)*(kVec1-1)+iOff1
            call dGemm_('T','N',nOcc(iSymI),nOcc(iSymK),nVir(iSymB),Four,Wrk(kU+iOff),nVir(iSymB),Wrk(kLia+iOff),nVir(iSymB),One, &
                        Wrk(kWij(iSymK)),nOcc(iSymI))
          end if
        end do
      end do

      ! Calculate "Coulomb"-contribution to the Frozen Lagr(III)
      ! --------------------------------------------------------
      iOffU1 = 0
      iOffL1 = 0
      do iSymJ=1,nSym
        iSymA = Mul(iSym,iSymJ)
        iSymI = iSymA
        if (nFro(iSymI)*nVir(iSymA)*nOcc(iSymJ) /= 0) then
          do iJ=1,nOcc(iSymJ)
            iOffU = (iJ-1)*nVir(iSymA)+iOffU1
            iOffL = (iJ-1)*nFro(iSymI)+iOffL1
            call dGemm_('N','T',nVir(iSymA),nFro(iSymI),NumVecK,-Four,Wrk(kU+iOffU),nMoMo(iSym,iVecOV),Wrk(kLiK+iOffL), &
                        nMoMo(iSym,iVecOF),One,Wrk(kFLagr(iSymA)),nVir(iSymA))
          end do
        end if
        iOffU1 = iOffU1+nOcc(iSymJ)*nVir(iSymA)
        iOffL1 = iOffL1+nOcc(iSymJ)*nFro(iSymI)
      end do

      ! Calculate "Coulomb"-contribution to the WaK (vir-froz)
      ! --------------------------------------------------------
      iOffU1 = 0
      iOffL1 = 0
      do iSymJ=1,nSym
        iSymA = Mul(iSym,iSymJ)
        iSymI = iSymA
        if (nFro(iSymI)*nVir(iSymA)*nOcc(iSymJ) /= 0) then
          do iJ=1,nOcc(iSymJ)
            iOffU = (iJ-1)*nVir(iSymA)+iOffU1
            iOffL = (iJ-1)*nFro(iSymI)+iOffL1
            call dGemm_('N','T',nVir(iSymA),nFro(iSymI),NumVecK,Eight,Wrk(kU+iOffU),nMoMo(iSym,iVecOV),Wrk(kLiK+iOffL), &
                        nMoMo(iSym,iVecOF),One,Wrk(kWaK(iSymA)),nVir(iSymA))
          end do
        end if
        iOffU1 = iOffU1+nOcc(iSymJ)*nVir(iSymA)
        iOffL1 = iOffL1+nOcc(iSymJ)*nFro(iSymI)
      end do
    end do

    call ChoMP2_OpenF(iClos,iTypR,iSym)
    call ChoMP2_OpenF(iClos,iTypL,iSym)
    call DaClos(LuUVec)
    call DaClos(LuVVec)

  end if                 ! vector check
end do !iSym

!******************************************************
! --------------------------------------------------  *
! Calculate the exchange part needing U_ic^K-vectors  *
! --------------------------------------------------  *
!******************************************************
! Calculate the max number of orbitals of a specific type
! -------------------------------------------------------
nFroMax = 0
nVirMax = 0
nOccMax = 0
do iSym=1,nSym
  nFroMax = max(nFroMax,nFro(iSym))
  nVirMax = max(nVirMax,nVir(iSym))
  nOccMax = max(nOccMax,nOcc(iSym))
end do

if (nVirMax > 1) then
  nB = min(nVirMax,maxvalue)
else
  nB = nVirMax
end if

! Allocate memory for T_[i]j^[b]c
! -------------------------------
lAmp = nVirMax*nB*nOccMax
kAmp = kEndFLagr
kEndAmp = kAmp+lAmp
Wrk(kAmp:kEndAmp-1) = Zero

! Allocate memory for Rjc^K
! -------------------------
lRjc = nVirMax*nVec
kRjc = kEndAmp
kEndRjc = kRjc+lRjc

! Allocate memory for Rib^K
! -------------------------
lRib = nOccMax*nVec
kRib = kEndRjc
kEndRib = kRib+lRib

! Open files for Reordered R-vectors
! ----------------------------------
iSeed = 7
LuRInv(1) = IsFreeUnit(iSeed)
write(Fname,'(A4,I1)') 'TMPV',4
call DaName_MF_WA(LuRInv(1),Fname)

iSeed = 8
LuRInv(2) = IsFreeUnit(iSeed)
write(Fname,'(A4,I1)') 'TMPV',5
call DaName_MF_WA(LuRInv(2),Fname)

do jSym=1,nSym
  if ((nMoMo(jSym,iVecOV) <= 0) .or. (nMp2Vec(jSym)*NumCho(jSym) <= 0)) cycle

  ! Allocate memory for Ric^J
  ! -------------------------
  lRic = nMoMo(jSym,iVecOV)*nVec
  kRic = kEndRib
  kEndRic = kRic+lRic

  ! Allocate memory for Lic^J
  ! -------------------------
  lLic = nMoMo(jSym,iVecOV)*nVec
  kLic = kEndRic
  kEndLic = kLic+lLic

  ! Allocate memory for Lab^J
  ! -------------------------
  lLab = nMoMo(jSym,iVecVV)*nVec
  kLab = kEndLic
  kEndLab = kLab+lLab

  ! Allocate memory for Lji^J
  ! -------------------------
  lLji = nMoMo(jSym,iVecOO)*nVec
  kLji = kEndLab
  kEndLji = kLji+lLji

  ! Allocate memory for LKa^J
  ! -------------------------
  lLKa = nMoMo(jSym,iVecFV)*nVec
  kLKa = kEndLji
  kEndLKa = kLKa+lLKa

  ! Allocate memory for LiK^J
  ! -------------------------
  lLiK = nMoMo(jSym,iVecOF)*nVec
  kLiK = kEndLKa
  kEndLiK = kLiK+lLiK

  ! Allocate memory for Ujb^J
  ! -------------------------
  lU = nMoMo(jSym,iVecOV)*nVec
  kU = kEndLiK
  kEndU = kU+lU

  ! Allocate memory for Vjb^J
  ! -------------------------
  lV = nMoMo(jSym,iVecOV)*nVec
  kV = kEndU

  nBatR = (nMP2Vec(jSym)-1)/nVec+1
  nBatL = (NumCho(jSym)-1)/nVec+1

  ! Open the File for reordered amplitude vectors
  ! ---------------------------------------------
  call ChoMP2_OpenF(1,iTypR,jSym)
  call ChoMP2_OpenF(1,iTypL,jSym)

  iSeed = 7
  LuUVec = IsFreeUnit(iSeed)
  write(Fname,'(A4,I1)') 'TMPV',1
  call DaName_MF_WA(LuUVec,Fname)

  iSeed = 7
  LuVVec = IsFreeUnit(iSeed)
  write(Fname,'(A4,I1)') 'TMPV',2
  call DaName_MF_WA(LuVVec,Fname)

  nBatMax = max(nBatR,nBatL)
  do jBat=1,nBatMax
    Wrk(kU:kU+lU-1) = Zero
    Wrk(kV:kV+lV-1) = Zero
    if (jBat == nBatR) then
      NumRVecJ = nMP2Vec(jSym)-nVec*(nBatR-1)
    else if (jBat > nBatR) then
      NumRVecJ = 0
    else
      NumRVecJ = nVec
    end if
    if (jBat == nBatL) then
      NumLVecJ = NumCho(jSym)-nVec*(nBatL-1)
    else if (jBat > nBatL) then
      NumLVecJ = 0
    else
      NumLVecJ = nVec
    end if
    jVec = nVec*(jBat-1)+1

    ! Read Amplitude vectors from kBat
    ! --------------------------------
    if (NumRVecJ > 0) then
      iOpt = 2
      lTot = nMoMo(jSym,iVecOV)*NumRVecJ
      iAdr = nMoMo(jSym,iVecOV)*(jVec-1)+1
      call dDaFile(lUnit_F(jSym,iTypR),iOpt,Wrk(kRic),lTot,iAdr)
    end if
    if (NumLVecJ > 0) then
      iOpt = 2
      lTot = nMoMo(jSym,iVecOV)*NumLVecJ
      iAdr = nMoMo(jSym,iVecOV)*(jVec-1)+1
      call dDaFile(lUnit_F(jSym,iTypL),iOpt,Wrk(kLic),lTot,iAdr)
    end if

    do iSymJ=1,nSym

      iSymB = Mul(jSym,iSymJ)
      if (nVir(iSymB) == 0) cycle
      nBBlock = (nVir(iSymB)-1)/nB+1

      do iSymI=1,nSym
        iSymC = Mul(jSym,iSymI)
        iSymJC = Mul(iSymJ,iSymC)
        NumIC = nOcc(iSymI)*nVir(iSymC)
        if (nOcc(iSymI)*nVir(iSymC)*(NumRVecJ+NumLVecJ) == 0) cycle
        do iJ=1,nOcc(iSymJ)
          do iBBlock=1,nBBlock
            if (iBBlock == nBBlock) then
              NumB = nVir(iSymB)-nB*(nBBlock-1)
            else
              NumB = nB
            end if
            iB1 = nB*(iBBlock-1)+1

            do kBat=1,nBatR
              if (kBat == nBatR) then
                NumVecK = nMP2Vec(iSymJC)-nVec*(nBatR-1)
              else
                NumVecK = nVec
              end if
              kVec = nVec*(kBat-1)+1

              iOpt = 2
              lTot = NumVecK*nVir(iSymC)
              iAdr = AdrR1(iSymC,iSymJ,iJ)+(kVec-1)*nVir(iSymC)
              call dDaFile(LuRInv(1),iOpt,Wrk(kRjc),lTot,iAdr)

              do iBRel=1,NumB
                iB = iBrel+(iB1-1)
                iOpt = 2
                lTot = NumVecK*nOcc(iSymI)
                iAdr = AdrR2(iSymB,iSymI,iB)+(kVec-1)*nOcc(iSymI)

                call dDaFile(LuRInv(2),iOpt,Wrk(kRib),lTot,iAdr)
                iOffAmp = (iBrel-1)*nVir(iSymC)*nOcc(iSymI)
                Fac = X(min((kBat-1),1))
                call dGemm_('N','T',nVir(iSymC),nOcc(iSymI),NumVecK,One,Wrk(kRjc),nVir(iSymC),Wrk(kRib),nOcc(iSymI),Fac, &
                            Wrk(kAmp+iOffAmp),nVir(iSymC))
              end do
            end do

            ! Calculate the contribution to U^J_jb
            ! ------------------------------------
            if (NumRVecJ > 0) then
              iOffRic = iT1am(iSymC,iSymI)
              iOffU = iT1am(iSymB,iSymJ)+(iJ-1)*nVir(iSymB)+(iB1-1)
              call dGemm_('T','N',NumB,NumRVecJ,NumIC,One,Wrk(kAmp),NumIC,Wrk(kRic+iOffRic),nMoMo(jSym,iVecOV),One,Wrk(kU+iOffU), &
                          nMoMo(jSym,iVecOV))
            end if
            ! Calculate the contribution to V^J_jb
            ! ------------------------------------
            if (NumLVecJ > 0) then
              iOffLic = iT1am(iSymC,iSymI)
              iOffV = iT1am(iSymB,iSymJ)+(iJ-1)*nVir(iSymB)+(iB1-1)
              call dGemm_('T','N',NumB,NumLVecJ,NumIC,One,Wrk(kAmp),NumIC,Wrk(kLic+iOffLic),nMoMo(jSym,iVecOV),One,Wrk(kV+iOffV), &
                          nMoMo(jSym,iVecOV))
            end if

          end do

        end do        !iJ
      end do          !iSymI
    end do            !iSymJ

    iOpt = 1
    lTot = nMoMo(jSym,iVecOV)*NumRVecJ
    iAdr = 1+nMoMo(jSym,iVecOV)*(jVec-1)
    call dDaFile(LuUVec,iOpt,Wrk(kU),lTot,iAdr)

    iOpt = 1
    lTot = nMoMo(jSym,iVecOV)*NumLVecJ
    iAdr = 1+nMoMo(jSym,iVecOV)*(jVec-1)
    call dDaFile(LuVVec,iOpt,Wrk(kV),lTot,iAdr)
  end do              !jBat

  do jBat=1,nBatMax
    if (jBat == nBatR) then
      NumRVecJ = nMP2Vec(jSym)-nVec*(nBatR-1)
    else if (jBat > nBatR) then
      NumRVecJ = 0
    else
      NumRVecJ = nVec
    end if
    if (jBat == nBatL) then
      NumLVecJ = NumCho(jSym)-nVec*(nBatL-1)
    else if (jBat > nBatL) then
      NumLVecJ = 0
    else
      NumLVecJ = nVec
    end if
    jVec = nVec*(jBat-1)+1

    ! Read Amplitude Rja
    ! --------------------------------
    if (NumRVecJ > 0) then
      iOpt = 2
      lTot = nMoMo(jSym,iVecOV)*NumRVecJ
      iAdr = nMoMo(jSym,iVecOV)*(jVec-1)+1
      call dDaFile(lUnit_F(jSym,iTypR),iOpt,Wrk(kRic),lTot,iAdr)
    end if
    ! Read Integrals R_ab^J from disk
    ! -------------------------------
    if (NumLVecJ > 0) then
      iOpt = 2
      lTot = nMoMo(jSym,iVecVV)*NumLVecJ
      iAdr = 1+nMoMo(jSym,iVecVV)*(jVec-1)+iAdrOff(jSym,iVecVV)
      call dDaFile(lUnit_F(jSym,iTypL),iOpt,Wrk(kLab),lTot,iAdr)
    end if

    ! Read Integrals L_ab^J from disk
    ! -------------------------------
    if (NumLVecJ > 0) then
      iOpt = 2
      lTot = nMoMo(jSym,iVecVV)*NumLVecJ
      iAdr = 1+nMoMo(jSym,iVecVV)*(jVec-1)+iAdrOff(jSym,iVecVV)
      call dDaFile(lUnit_F(jSym,iTypL),iOpt,Wrk(kLab),lTot,iAdr)
    end if

    ! Read L_ji^J-vectors from disk
    ! -----------------------------
    if (NumLVecJ > 0) then
      iOpt = 2
      lTot = nMoMo(jSym,iVecOO)*NumLVecJ
      iAdr = 1+nMoMo(jSym,iVecOO)*(jVec-1)+iAdrOff(jSym,iVecOO)
      call dDaFile(lUnit_F(jSym,iTypL),iOpt,Wrk(kLji),lTot,iAdr)
    end if

    ! Read L_ic^J-vectors from disk
    ! -----------------------------
    if (NumLVecJ > 0) then
      iOpt = 2
      lTot = nMoMo(jSym,iVecOV)*NumLVecJ
      iAdr = nMoMo(jSym,iVecOV)*(jVec-1)+1
      call dDaFile(lUnit_F(jSym,iTypL),iOpt,Wrk(kLic),lTot,iAdr)
    end if

    ! Read L_Ka^J-vectors from disk
    ! -----------------------------
    if (NumLVecJ > 0) then
      iOpt = 2
      lTot = nMoMo(jSym,iVecFV)*NumLVecJ
      iAdr = 1+nMoMo(jSym,iVecFV)*(jVec-1)+iAdrOff(jSym,iVecFV)
      call dDaFile(lUnit_F(jSym,iTypL),iOpt,Wrk(kLKa),lTot,iAdr)
    end if

    ! Read L_Ka^J-vectors from disk
    ! -----------------------------
    if (NumLVecJ > 0) then
      iOpt = 2
      lTot = nMoMo(jSym,iVecOF)*NumLVecJ
      iAdr = 1+nMoMo(jSym,iVecOF)*(jVec-1)+iAdrOff(jSym,iVecOF)
      call dDaFile(lUnit_F(jSym,iTypL),iOpt,Wrk(kLiK),lTot,iAdr)
    end if

    ! Read U_jb^J-vectors from disk
    ! -----------------------------
    if (NumRVecJ > 0) then
      iOpt = 2
      lTot = nMoMo(jSym,iVecOV)*NumRVecJ
      iAdr = nMoMo(jSym,iVecOv)*(jVec-1)+1
      call dDaFile(LuUVec,iOpt,Wrk(kU),lTot,iAdr)
    end if

    ! Read V_jb^J-vectors from disk
    ! -----------------------------
    if (NumLVecJ > 0) then
      iOpt = 2
      lTot = nMoMo(jSym,iVecOV)*NumLVecJ
      iAdr = nMoMo(jSym,iVecOv)*(jVec-1)+1
      call dDaFile(LuVVec,iOpt,Wrk(kV),lTot,iAdr)
    end if

    ! Calculate the "Exchange" Contribution to Pab
    ! -------------------------------------------

    iOff1 = 0
    do iSymJ=1,nSym
      iSymA = Mul(jSym,iSymJ)
      if (nOcc(iSymJ)*nVir(iSymA)*NumRVecJ /= 0) then
        do iI=1,nOcc(iSymJ)
          iOff = (iI-1)*nVir(iSymA)+iOff1
          call dGemm_('N','T',nVir(iSymA),nVir(iSymA),NumRVecJ,-Two,Wrk(kRic+iOff),nMoMo(jSym,iVecOV),Wrk(kU+iOff), &
                      nMoMo(jSym,iVecOV),One,Wrk(kPab(iSymA)),nVir(iSymA))
        end do
      end if
      iOff1 = iOff1+nOcc(iSymJ)*nVir(iSymA)
    end do

    ! Calculate the "Exchange" Contribution to Pij
    ! -------------------------------------------

    do jVec1=1,NumRVecJ
      iOff1 = 0
      do iSymJ=1,nSym
        iSymI = iSymJ
        iSymB = Mul(iSymJ,jSym)
        if (nOcc(iSymJ)*nVir(iSymB) /= 0) then
          iOff = nMoMo(jSym,iVecOV)*(jVec1-1)+iOff1
          call dGemm_('T','N',nOcc(iSymJ),nOcc(iSymI),nVir(iSymB),Two,Wrk(kU+iOff),nVir(iSymB),Wrk(kRic+iOff),nVir(iSymB),One, &
                      Wrk(kPij(iSymJ)),nOcc(iSymJ))
        end if
        iOff1 = iOff1+nOcc(iSymJ)*nVir(iSymB)
      end do
    end do

    ! Calculate the "Exchange"-contribution to PiK (Froz-occ)
    ! -------------------------------------------------------

    do jVec1=1,NumLVecJ
      iOffV1 = 0
      iOffL1 = 0
      do iSymK=1,nSym
        iSymI = iSymK
        iSymA = Mul(iSymK,jSym)
        if (nOcc(iSymI)*nFro(iSymK)*nVir(iSymA) /= 0) then
          iOffV = nMoMo(jSym,iVecOV)*(jVec1-1)+iOffV1
          iOffL = nMoMo(jSym,iVecFV)*(jVec1-1)+iOffL1
          call dGemm_('T','N',nOcc(iSymI),nFro(iSymK),nVir(iSymA),Two,Wrk(kV+iOffV),nVir(iSymA),Wrk(kLKa+iOffL),nVir(iSymA),One, &
                      Wrk(kPiK(iSymK)),nOcc(iSymI))
        end if
        iOffV1 = iOffV1+nOcc(iSymI)*nVir(iSymA)
        iOffL1 = iOffL1+nFro(iSymK)*nVir(iSymA)
      end do
    end do

    ! Calculate the "Exchange"-contribution to WiK (Froz-occ)
    ! -------------------------------------------------------

    do jVec1=1,NumLVecJ
      iOffV1 = 0
      iOffL1 = 0
      do iSymK=1,nSym
        iSymI = iSymK
        iSymA = Mul(iSymK,jSym)
        if (nOcc(iSymI)*nFro(iSymK)*nVir(iSymA) /= 0) then
          iOffV = nMoMo(jSym,iVecOV)*(jVec1-1)+iOffV1
          iOffL = nMoMo(jSym,iVecFV)*(jVec1-1)+iOffL1
          call dGemm_('T','N',nOcc(iSymI),nFro(iSymK),nVir(iSymA),-Two,Wrk(kV+iOffV),nVir(iSymA),Wrk(kLKa+iOffL),nVir(iSymA),One, &
                      Wrk(kWiK(iSymK)),nOcc(iSymI))
        end if
        iOffV1 = iOffV1+nOcc(iSymI)*nVir(iSymA)
        iOffL1 = iOffL1+nFro(iSymK)*nVir(iSymA)
      end do
    end do

    ! Calculate the "Exchange"-contribution to Lagr(IV)
    ! -------------------------------------------------

    do jVec1=1,NumLVecJ
      iOffV1 = 0
      iOffL1 = 0
      do iSymI=1,nSym
        iSymA = iSymI
        iSymB = Mul(iSymI,jSym)
        if (nOcc(iSymI)*nVir(iSymB)*nVir(iSymA) /= 0) then
          iOffV = nMoMo(jSym,iVecOV)*(jVec1-1)+iOffV1
          iOffL = nMoMo(jSym,iVecVV)*(jVec1-1)+iOffL1
          call dGemm_('T','N',nVir(iSymA),nOcc(iSymI),nVir(iSymB),-Two,Wrk(kLab+iOffL),nVir(iSymB),Wrk(kV+iOffV),nVir(iSymB),One, &
                      Wrk(kLagr(iSymA)),nVir(iSymA))
        end if
        iOffL1 = iOffL1+nVir(iSymA)*nVir(iSymB)
        iOffV1 = iOffV1+nOcc(iSymI)*nVir(iSymB)
      end do
    end do

    ! Calculate the "Exchange"-contribution to Lagr(III)
    ! --------------------------------------------------

    iOffV1 = 0
    iOffL1 = 0
    do iSymJ=1,nSym
      iSymA = Mul(jSym,iSymJ)
      iSymI = iSymA
      if (nOcc(iSymI)*nVir(iSymA)*nOcc(iSymJ) /= 0) then
        do iJ=1,nOcc(iSymJ)
          iOffV = (iJ-1)*nVir(iSymA)+iOffV1
          iOffL = (iJ-1)*nOcc(iSymI)+iOffL1
          call dGemm_('N','T',nVir(iSymA),nOcc(iSymI),NumLVecJ,Two,Wrk(kV+iOffV),nMoMo(jSym,iVecOV),Wrk(kLji+iOffL), &
                      nMoMo(jSym,iVecOO),One,Wrk(kLagr(iSymA)),nVir(iSymA))
        end do
      end if
      iOffV1 = iOffV1+nOcc(iSymJ)*nVir(iSymA)
      iOffL1 = iOffL1+nOcc(iSymJ)*nOcc(iSymI)
    end do

    ! Calculate the "Exchange"-contribution to Wai
    ! --------------------------------------------------

    iOffV1 = 0
    iOffL1 = 0
    do iSymJ=1,nSym
      iSymA = Mul(jSym,iSymJ)
      iSymI = iSymA
      if (nOcc(iSymI)*nVir(iSymA)*nOcc(iSymJ)*NumLVecJ /= 0) then
        do iJ=1,nOcc(iSymJ)
          iOffV = (iJ-1)*nVir(iSymA)+iOffV1
          iOffL = (iJ-1)*nOcc(iSymI)+iOffL1
          call dGemm_('N','T',nVir(iSymA),nOcc(iSymI),NumLVecJ,-Four,Wrk(kV+iOffV),nMoMo(jSym,iVecOV),Wrk(kLji+iOffL), &
                      nMoMo(jSym,iVecOO),One,Wrk(kWai(iSymA)),nVir(iSymA))
        end do
      end if
      iOffV1 = iOffV1+nOcc(iSymJ)*nVir(iSymA)
      iOffL1 = iOffL1+nOcc(iSymJ)*nOcc(iSymI)
    end do

    ! Calculate the "Exchange"-contribution to Wab
    ! ---------------------------------------------
    iOff1 = 0
    do iSymJ=1,nSym
      iSymA = Mul(jSym,iSymJ)
      if (nOcc(iSymJ)*nVir(iSymA)*NumLVecJ /= 0) then
        do iI=1,nOcc(iSymJ)
          iOff = (iI-1)*nVir(iSymA)+iOff1
          call dGemm_('N','T',nVir(iSymA),nVir(iSymA),NumLVecJ,-Two,Wrk(kLic+iOff),nMoMo(jSym,iVecOV),Wrk(kV+iOff), &
                      nMoMo(jSym,iVecOV),One,Wrk(kWab(iSymA)),nVir(iSymA))
        end do
      end if
      iOff1 = iOff1+nOcc(iSymJ)*nVir(iSymA)
    end do

    ! Calculate the "Exchange"-contribution to Wij
    ! --------------------------------------------

    do jVec1=1,NumLVecJ
      iOff1 = 0
      do iSymJ=1,nSym
        iSymI = iSymJ
        iSymB = Mul(iSymJ,jSym)
        if (nOcc(iSymJ)*nVir(iSymB) /= 0) then
          iOff = nMoMo(jSym,iVecOV)*(jVec1-1)+iOff1
          call dGemm_('T','N',nOcc(iSymJ),nOcc(iSymI),nVir(iSymB),-Two,Wrk(kV+iOff),nVir(iSymB),Wrk(kLic+iOff),nVir(iSymB),One, &
                      Wrk(kWij(iSymJ)),nOcc(iSymJ))
        end if
        iOff1 = iOff1+nOcc(iSymJ)*nVir(iSymB)
      end do
    end do

    ! Calculate "Exchange"-contribution to the Frozen Lagr(III)
    ! --------------------------------------------------------
    iOffV1 = 0
    iOffL1 = 0
    do iSymJ=1,nSym
      iSymA = Mul(jSym,iSymJ)
      iSymI = iSymA
      if (nFro(iSymI)*nVir(iSymA)*nOcc(iSymJ) /= 0) then
        do iJ=1,nOcc(iSymJ)
          iOffV = (iJ-1)*nVir(iSymA)+iOffV1
          iOffL = (iJ-1)*nFro(iSymI)+iOffL1
          call dGemm_('N','T',nVir(iSymA),nFro(iSymI),NumLVecJ,Two,Wrk(kV+iOffV),nMoMo(jSym,iVecOV),Wrk(kLiK+iOffL), &
                      nMoMo(jSym,iVecOF),One,Wrk(kFLagr(iSymA)),nVir(iSymA))
        end do
      end if
      iOffV1 = iOffV1+nOcc(iSymJ)*nVir(iSymA)
      iOffL1 = iOffL1+nOcc(iSymJ)*nFro(iSymI)
    end do

    ! Calculate "Exchange"-contribution to the Wak(Vir-Froz)
    ! --------------------------------------------------------
    iOffV1 = 0
    iOffL1 = 0
    do iSymJ=1,nSym
      iSymA = Mul(jSym,iSymJ)
      iSymI = iSymA
      if (nFro(iSymI)*nVir(iSymA)*nOcc(iSymJ) /= 0) then
        do iJ=1,nOcc(iSymJ)
          iOffV = (iJ-1)*nVir(iSymA)+iOffV1
          iOffL = (iJ-1)*nFro(iSymI)+iOffL1
          call dGemm_('N','T',nVir(iSymA),nFro(iSymI),NumLVecJ,-Four,Wrk(kV+iOffV),nMoMo(jSym,iVecOV),Wrk(kLiK+iOffL), &
                      nMoMo(jSym,iVecOF),One,Wrk(kWaK(iSymA)),nVir(iSymA))
        end do
      end if
      iOffV1 = iOffV1+nOcc(iSymJ)*nVir(iSymA)
      iOffL1 = iOffL1+nOcc(iSymJ)*nFro(iSymI)
    end do

  end do

  call DaClos(LuUVec)
  call DaClos(LuVVec)
  call ChoMP2_OpenF(iClos,iTypR,jSym)
  call ChoMP2_OpenF(iClos,iTypL,jSym)
end do            ! jSym

call DaClos(LuRInv(1))
call DaClos(LuRInv(2))

! Scale the Frozen-Occupied part of the density
! ---------------------------------------------
do iSym=1,nSym
  do iI=1,nOcc(iSym)
    E_i = EOcc(iOcc(iSym)+iI)
    do iK=1,nFro(iSym)
      E_K = EFro(iFro(iSym)+iK)
      i = iI-1+(iK-1)*nOcc(iSym)
      Wrk(kPik(iSym)+i) = Wrk(kPik(iSym)+i)/(E_i-E_k)
    end do
  end do
end do

do iSym=1,nSym
  if (NumCho(iSym) > 0) then

    call ChoMP2_OpenF(1,1,iSym)

    ! I doubt the intermediate files here need to be
    ! Multifiles, probably they will never be used again after the
    ! end of the iSym-loop (or even IB-block loop).
    iSeed = 7
    LuUVec = IsFreeUnit(iSeed)
    write(Fname,'(A4,I1)') 'TMPV',1
    call DaName_MF_WA(LuUVec,Fname)

    iSeed = 7
    LuVVec = IsFreeUnit(iSeed)
    write(Fname,'(A4,I1)') 'TMPV',2
    call DaName_MF_WA(LuVVec,Fname)

    iSeed = 7
    LuWVec = IsFreeUnit(iSeed)
    write(Fname,'(A4,I1)') 'TMPV',3
    call DaName_MF_WA(LuWVec,Fname)

    lScr = lWrk-kEndFLagr
    ! Construct Lagr(i)
    ! -----------------
    call ChoMP2g_ConstrAP(irc,Wrk(kEndFLagr),lScr,'oovo',iSym,nVec,Wrk(kLagr(1)),lLagr,Wrk(kPij(1)),lPij,-Half)

    call ChoMP2g_ConstrAP(irc,Wrk(kEndFLagr),lScr,'fovo',iSym,nVec,Wrk(kLagr(1)),lLagr,Wrk(kPiK(1)),lPiK,-One)
    ! Construct FLagr(i)
    ! ------------------
    call ChoMP2g_ConstrAP(irc,Wrk(kEndFLagr),lScr,'oovf',iSym,nVec,Wrk(kFLagr(1)),lFLagr,Wrk(kPij(1)),lPij,-Half)

    call ChoMP2g_ConstrAP(irc,Wrk(kEndFLagr),lScr,'fovf',iSym,nVec,Wrk(kFLagr(1)),lFLagr,Wrk(kPiK(1)),lPiK,-One)
    ! Construct Lagr(ii)
    ! ------------------
    call ChoMP2g_ConstrAP(irc,Wrk(kEndFLagr),lScr,'vvvo',iSym,nVec,Wrk(kLagr(1)),lLagr,Wrk(kPab(1)),lPab,-Half)
    ! Construct FLagr(ii)
    ! -------------------
    call ChoMP2g_ConstrAP(irc,Wrk(kEndFLagr),lScr,'vvvf',iSym,nVec,Wrk(kFLagr(1)),lFLagr,Wrk(kPab(1)),lPab,-Half)

    call ChoMP2_OpenF(iClos,1,iSym)
    call DaClos(LuUVec)
    call DaClos(LuVVec)
    call DaClos(LuWVec)

  end if
end do

#ifdef _DEBUGPRINT_
write(u6,*) 'Pab'
do i=1,lPab
  write(u6,*) Wrk(kPab(1)+i-1)
end do
write(u6,*) 'lWab'
do i=1,lWab
  write(u6,*) Wrk(kWab(1)+i-1)
end do
write(u6,*) 'Pij'
do i=1,lPij
  write(u6,*) Wrk(kPij(1)+i-1)
end do
write(u6,*) 'Wij'
do i=1,lWij
  write(u6,*) Wrk(kWij(1)+i-1)
end do
write(u6,*) 'WiK'
do i=1,lWiK
  write(u6,*) Wrk(kWiK(1)+i-1)
end do
write(u6,*) 'WaK'
do i=1,lWaK
  write(u6,*) Wrk(kWaK(1)+i-1)
end do
write(u6,*) 'Wai'
do i=1,lWai
  write(u6,*) Wrk(kWai(1)+i-1)
end do
write(u6,*) 'PiK'
do i=1,lPiK
  write(u6,*) Wrk(kPiK(1)+i-1)
end do

write(u6,*) 'Lagr'
do i=1,lLagr
  write(u6,*) Wrk(kLagr(1)+i-1)
end do

write(u6,*) 'FLagr'
do i=1,lFLagr
  write(u6,*) Wrk(kFLagr(1)+i-1)
end do
#endif

end subroutine ChoMP2g_density1
