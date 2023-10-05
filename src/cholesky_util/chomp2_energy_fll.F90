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
! Copyright (C) 2004, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_Energy_Fll(irc,Delete,EMP2,EOcc,EVir,Wrk,lWrk)
!
! Thomas Bondo Pedersen, Dec. 2004.
!
! Purpose: compute MP2 energy contribution using original
!          Cholesky vectors on disk for the case nBatch=1.

use Symmetry_Info, only: Mul
use Index_Functions, only: iTri
use Cholesky, only: nSym, NumCho
use ChoMP2, only: ChoAlg, DecoMP2, iMatab, iT1am, LiMatij, lUnit_F, nBatch, nMatab, nMP2Vec, nOcc, nT1am, nVir, Wref
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
logical(kind=iwp), intent(in) :: Delete
integer(kind=iwp), intent(in) :: lWrk
real(kind=wp), intent(out) :: EMP2, Wrk(lWrk)
real(kind=wp), intent(in) :: EOcc(*), EVir(*)
integer(kind=iwp) :: i, iAdr, iBat, iClos, ij, iOpt, iSym, iSyma, iSymab, iSymb, iSymi, iSymij, iSymj, iTyp, iVaJi(8), iVec, &
                     iVec0, iVec1, j, kEnd0, kEnd1, kEnd2, kMabij, kXaibj, kOff1, kOff2, kOffi, kOffj, kOffM, kVec, kVecai, kXint, &
                     LiT2am(8), LnT2am, lTot, lWrk0, lWrk1, lWrk2, Nai, nBat, nEnrVec(8), NumV, NumVec, nVaJi, nVec
real(kind=wp) :: Fac
real(kind=wp), parameter :: X(0:1) = [Zero,One]
character(len=*), parameter :: SecNam = 'ChoMP2_Energy_Fll'

if (nBatch /= 1) then
  irc = -1
  return
end if

irc = 0

! Determine if vector files are to be deleted after use.
! ------------------------------------------------------

if (Delete) then
  iClos = 3
else
  iClos = 2
end if

! Set number and type of vectors.
! -------------------------------

if (DecoMP2) then
  iTyp = 2
  nEnrVec(1:nSym) = nMP2Vec(1:nSym)
else
  iTyp = 1
  nEnrVec(1:nSym) = NumCho(1:nSym)
end if

! Set (ai|bj) indices.
! --------------------

call ChoMP2_Energy_GetInd(LnT2am,LiT2am,1,1)

! Allocate memory for integrals.
! ------------------------------

kXaibj = 1
kEnd0 = kXaibj+LnT2am
lWrk0 = lWrk-kEnd0+1
if (lWrk0 < 0) call SysAbendMsg(SecNam,'insufficient memory','[0]')

! Initialize MP2 energy correction.
! ---------------------------------

EMP2 = Zero
Wref = Zero

! Special code for ChoAlg=2:
! compute M(ab,ij) = (ai|bj) with i<=j using level 3 BLAS.
! For ChoAlg=1: use strictly lower triangular storage (=>
! level 2 BLAS).
! --------------------------------------------------------

if (ChoAlg == 2) then ! level 3 BLAS algorithm

  kMabij = kXaibj ! rename pointer
  Wrk(kMabij:kMabij+LnT2am-1) = Zero ! initialize

  ! Loop over Cholesky symmetries.
  ! ------------------------------

  do iSym=1,nSym

    Nai = nT1am(iSym)
    if ((Nai > 0) .and. (nEnrVec(iSym) > 0)) then

      ! Reserve memory for reading a single vector.
      ! -------------------------------------------

      kVecai = kEnd0
      kEnd1 = kVecai+Nai
      lWrk1 = lWrk-kEnd1+1

      if (lWrk1 < Nai) call SysAbendMsg(SecNam,'Insufficient memory','[ChoAlg.2.1]')

      ! Set up batch over Cholesky vectors.
      ! -----------------------------------

      nVec = min(lWrk1/Nai,nEnrVec(iSym))
      if (nVec < 1) call SysAbendMsg(SecNam,'Insufficient memory','[ChoAlg.2.2]') ! should not happen
      nBat = (nEnrVec(iSym)-1)/nVec+1

      ! Open Cholesky vector files.
      ! ---------------------------

      call ChoMP2_OpenF(1,iTyp,iSym)

      ! Start vector batch loop.
      ! ------------------------

      do iBat=1,nBat

        if (iBat == nBat) then
          NumVec = nEnrVec(iSym)-nVec*(nBat-1)
        else
          NumVec = nVec
        end if
        iVec1 = nVec*(iBat-1)+1

        ! Set up index arrays for reordered vectors.
        ! ------------------------------------------

        nVaJi = 0
        do iSymi=1,nSym
          iSyma = Mul(iSymi,iSym)
          iVaJi(iSymi) = nVaJi
          nVaJi = nVaJi+nVir(iSyma)*NumVec*nOcc(iSymi)
        end do

        ! Pointer to reordered vectors: kVec.
        ! -----------------------------------

        kVec = kEnd1
        kEnd2 = kVec+nVaJi
        lWrk2 = lWrk-kEnd2+1
        if (lWrk2 < 0) call SysAbendMsg(SecNam,'Insufficient memory','[ChoAlg.2.3]') ! should not happen

        ! Read one vector at a time and reorder.
        ! --------------------------------------

        iVec0 = iVec1-1
        do iVec=1,NumVec

          iOpt = 2
          lTot = nT1am(iSym)
          iAdr = nT1am(iSym)*(iVec0+iVec-1)+1
          call ddaFile(lUnit_F(iSym,iTyp),iOpt,Wrk(kVecai),lTot,iAdr)

          do iSymi=1,nSym
            iSyma = Mul(iSymi,iSym)
            do i=1,nOcc(iSymi)
              kOff1 = kVecai+iT1am(iSyma,iSymi)+nVir(iSyma)*(i-1)
              kOff2 = kVec+iVaJi(iSymi)+nVir(iSyma)*NumVec*(i-1)+nVir(iSyma)*(iVec-1)
              Wrk(kOff2:kOff2+nVir(iSyma)-1) = Wrk(kOff1:kOff1+nVir(iSyma)-1)
            end do
          end do

        end do

        ! Compute M(ab,ij) for i<=j.
        ! First do iSymi=iSymj, then iSymi<iSymj.
        ! ---------------------------------------

        do iSymj=1,nSym

          iSymb = Mul(iSymj,iSym)

          if (nVir(iSymb) > 0) then

            do j=1,nOcc(iSymj)
              do i=1,j

                ij = LiMatij(iSymj,iSymj,1)+iTri(i,j)

                kOffi = kVec+iVaJi(iSymj)+nVir(iSymb)*NumVec*(i-1)
                kOffj = kVec+iVaJi(iSymj)+nVir(iSymb)*NumVec*(j-1)
                kOffM = kMabij+LiT2am(1)+nMatab(1)*(ij-1)+iMatab(iSymb,iSymb)

                call DGEMM_('N','T',nVir(iSymb),nVir(iSymb),NumVec,ONe,Wrk(kOffi),nVir(iSymb),Wrk(kOffj),nVir(iSymb),One, &
                            Wrk(kOffM),nVir(iSymb))

              end do
            end do

            do iSymi=1,iSymj-1

              iSyma = Mul(iSymi,iSym)
              iSymab = Mul(iSyma,iSymb)
              iSymij = iSymab

              if ((nOcc(iSymi) > 0) .and. (nVir(iSyma) > 0)) then

                do j=1,nOcc(iSymj)
                  do i=1,nOcc(iSymi)

                    ij = LiMatij(iSymi,iSymj,1)+nOcc(iSymi)*(j-1)+i

                    kOffi = kVec+iVaJi(iSymi)+nVir(iSyma)*NumVec*(i-1)
                    kOffj = kVec+iVaJi(iSymj)+nVir(iSymb)*NumVec*(j-1)
                    kOffM = kMabij+LiT2am(iSymij)+nMatab(iSymab)*(ij-1)+iMatab(iSyma,iSymb)

                    call DGEMM_('N','T',nVir(iSyma),nVir(iSymb),NumVec,One,Wrk(kOffi),nVir(iSyma),Wrk(kOffj),nVir(iSymb),One, &
                                Wrk(kOffM),nVir(iSyma))

                  end do
                end do

              end if

            end do
          end if

        end do

      end do

      ! Close (and possibly delete) Cholesky vector files.
      ! --------------------------------------------------

      call ChoMP2_OpenF(iClos,iTyp,iSym)

    end if

  end do ! iSym

else ! level 2 BLAS algorithm

  ! Loop over Cholesky symmetries.
  ! ------------------------------

  do iSym=1,nSym
    if ((nT1am(iSym) > 0) .and. (nEnrVec(iSym) > 0)) then

      ! Set up Cholesky batch.
      ! ----------------------

      NumVec = min(lWrk0/nT1am(iSym),nEnrVec(iSym))
      if (NumVec < 1) call SysAbendMsg(SecNam,'insufficient memory','[2]')
      nBat = (nEnrVec(iSym)-1)/NumVec+1

      ! Open Cholesky vector file.
      ! --------------------------

      call ChoMP2_OpenF(1,iTyp,iSym)

      ! Start batch loop.
      ! -----------------

      do iBat=1,nBat

        if (iBat == nBat) then
          NumV = nEnrVec(iSym)-NumVec*(nBat-1)
        else
          NumV = NumVec
        end if
        iVec1 = NumVec*(iBat-1)+1

        ! Read vectors.
        ! -------------

        kVec = kEnd0

        iOpt = 2
        lTot = nT1am(iSym)*NumV
        iAdr = nT1am(iSym)*(iVec1-1)+1
        call ddaFile(lUnit_F(iSym,iTyp),iOpt,Wrk(kVec),lTot,iAdr)

        ! Compute contribution.
        ! ---------------------

        Fac = X(min((iBat-1),1))
        kXint = kXaibj+LiT2am(iSym)
        call dGeMM_Tri('N','T',nT1am(iSym),nT1am(iSym),NumV,One,Wrk(kVec),nT1am(iSym),Wrk(kVec),nT1am(iSym),Fac,Wrk(kXint), &
                       nT1am(iSym))

      end do

      ! Close (and possibly delete) Cholesky vector file.
      ! -------------------------------------------------

      call ChoMP2_OpenF(iClos,iTyp,iSym)

    end if
  end do

end if ! ChoAlg

! Compute energy contribution.
! ----------------------------

call ChoMP2_Energy_Contr(EMP2,EOcc,EVir,Wrk(kXaibj),LnT2am,LiT2am,1,1)

! Change sign on energy.
! ----------------------

EMP2 = -EMP2

end subroutine ChoMP2_Energy_Fll
