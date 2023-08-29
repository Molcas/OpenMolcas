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
! Copyright (C) 2008, Francesco Aquilante                              *
!***********************************************************************

subroutine ChoMP2_fno_Fll(irc,Delete,P_ab,P_ii,EOcc,EVir,Wrk,lWrk)
!
!  F. Aquilante, Geneva May 2008  (snick to Pedersen's code)

use Symmetry_Info, only: Mul
use Index_Functions, only: iTri
use Cholesky, only: nSym, NumCho
use ChoMP2, only: ChoAlg, DecoMP2, DeMP2, iOcc, iMatab, iT1am, iVir, LiMatij, lUnit_F, MP2_small, nBatch, nMP2Vec, nMatab, nOcc, &
                  nT1am, nVir, shf
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
logical(kind=iwp), intent(in) :: Delete
real(kind=wp), intent(inout) :: P_ab(*), P_ii(*)
real(kind=wp), intent(in) :: EOcc(*), EVir(*)
integer(kind=iwp), intent(in) :: lWrk
real(kind=wp), intent(out) :: Wrk(lWrk)
integer(kind=iwp) :: i, iAdr, iBat, iClos, ij, iOpt, iS, iSym, iSyma, iSymb, iSymi, iSymj, iTyp, iVaJi(8), iVec, iVec0, iVec1, j, &
                     ja, jb, kEnd0, kEnd1, kEnd2, kMabij, kOff1, kOff2, kOffi, kOffj, kOffM, kOffMM, kP(8), kVec, kVecai, kXaibj, &
                     LiT2am(8), LnT2am, lP(8), lTot, lWrk0, lWrk1, lWrk2, Nai, nBat, nEnrVec(8), NumVec, nVaJi, nVec
real(kind=wp) :: Dnom, xsDnom
character(len=*), parameter :: SecNam = 'ChoMP2_fno_Fll'
real(kind=wp), external :: ddot_

if (nBatch /= 1) then
  irc = -1
  return
end if

irc = 0

kP(1) = 1
lP(1) = 0
do iS=2,nSym
  kP(iS) = kP(iS-1)+nVir(iS-1)**2
  lP(iS) = lP(iS-1)+nOcc(iS-1)
end do

DeMP2 = Zero

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

if ((ChoAlg == 2) .and. MP2_small) then ! level 3 BLAS algorithm

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

        ! Compute M(ab,ii) .
        !-------------------

        do iSymj=1,nSym

          iSymb = Mul(iSymj,iSym)

          if (nVir(iSymb) > 0) then

            do j=1,nOcc(iSymj)

              i = j

              ij = LiMatij(iSymj,iSymj,1)+iTri(i,j)

              kOffi = kVec+iVaJi(iSymj)+nVir(iSymb)*NumVec*(i-1)
              kOffj = kVec+iVaJi(iSymj)+nVir(iSymb)*NumVec*(j-1)
              kOffM = kMabij+LiT2am(1)+nMatab(1)*(ij-1)+iMatab(iSymb,iSymb)

              call dGeMM_('N','T',nVir(iSymb),nVir(iSymb),NumVec,One,Wrk(kOffi),nVir(iSymb),Wrk(kOffj),nVir(iSymb),One,Wrk(kOffM), &
                          nVir(iSymb))

            end do

          end if

        end do

      end do

      call Cho_GAdGOp(Wrk(kMabij),LnT2am,'+')

      ! Close (and possibly delete) Cholesky vector files.
      ! --------------------------------------------------

      call ChoMP2_OpenF(iClos,iTyp,iSym)

      ! Compute Energy contribution.
      ! ----------------------------

      do iSymj=1,nSym

        iSymb = Mul(iSymj,iSym)

        if (nVir(iSymb) > 0) then

          do j=1,nOcc(iSymj)

            i = j

            ij = LiMatij(iSymj,iSymj,1)+iTri(i,j)

            kOffM = kMabij+LiT2am(1)+nMatab(1)*(ij-1)+iMatab(iSymb,iSymb)

            do jb=1,nVir(iSymb)
              do ja=1,nVir(iSymb)
                Dnom = EVir(iVir(iSymb)+ja)+EVir(iVir(iSymb)+jb)-Two*EOcc(iOcc(iSymj)+j)
                xsDnom = Dnom/(Dnom**2+shf**2) !Reg shf
                kOffMM = kOffM+nVir(iSymb)*(jb-1)+ja-1
                !DeMP2 = DeMP2+Wrk(kOffMM)**2/Dnom
                DeMP2 = DeMP2+Wrk(kOffMM)**2*xsDnom
              end do
            end do

          end do

        end if

      end do

    end if

  end do ! iSym

else if (ChoAlg == 2) then ! level 3 BLAS algorithm

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

        ! Compute M(ab,ii) .
        ! ------------------

        do iSymj=1,nSym

          iSymb = Mul(iSymj,iSym)

          if (nVir(iSymb) > 0) then

            do j=1,nOcc(iSymj)

              i = j

              ij = LiMatij(iSymj,iSymj,1)+iTri(i,j)

              kOffi = kVec+iVaJi(iSymj)+nVir(iSymb)*NumVec*(i-1)
              kOffj = kVec+iVaJi(iSymj)+nVir(iSymb)*NumVec*(j-1)
              kOffM = kMabij+LiT2am(1)+nMatab(1)*(ij-1)+iMatab(iSymb,iSymb)

              call dGeMM_('N','T',nVir(iSymb),nVir(iSymb),NumVec,One,Wrk(kOffi),nVir(iSymb),Wrk(kOffj),nVir(iSymb),One,Wrk(kOffM), &
                          nVir(iSymb))

            end do

          end if

        end do

      end do

      call Cho_GAdGOp(Wrk(kMabij),LnT2am,'+')

      ! Close (and possibly delete) Cholesky vector files.
      ! --------------------------------------------------

      call ChoMP2_OpenF(iClos,iTyp,iSym)

      do iSymj=1,nSym

        iSymb = Mul(iSymj,iSym)

        if (nVir(iSymb) > 0) then

          do j=1,nOcc(iSymj)

            i = j

            ij = LiMatij(iSymj,iSymj,1)+iTri(i,j)

            kOffM = kMabij+LiT2am(1)+nMatab(1)*(ij-1)+iMatab(iSymb,iSymb)

            ! Compute T(a,b)[i] and energy contrib.
            ! -------------------------------------
            do jb=1,nVir(iSymb)
              do ja=1,nVir(iSymb)
                Dnom = EVir(iVir(iSymb)+ja)+EVir(iVir(iSymb)+jb)-Two*EOcc(iOcc(iSymj)+j)
                xsDnom = Dnom/(Dnom**2+shf**2) !Reg shf
                kOffMM = kOffM+nVir(iSymb)*(jb-1)+ja-1
                !DeMP2 = DeMP2+Wrk(kOffMM)**2/Dnom
                DeMP2 = DeMP2+Wrk(kOffMM)**2*xsDnom
                !Wrk(kOffMM) = Wrk(kOffMM)/Dnom
                Wrk(kOffMM) = Wrk(kOffMM)*xsDnom
              end do
            end do

            P_ii(lP(iSymj)+i) = P_ii(lP(iSymj)+i)+ddot_(nVir(iSymb)**2,Wrk(kOffM),1,Wrk(kOffM),1)

            ! Compute P(a,b) += sum_c T(a,c)*T(c,b)
            ! -------------------------------------
            call dGeMM_('N','N',nVir(iSymb),nVir(iSymb),nVir(iSymb),One,Wrk(kOffM),nVir(iSymb),Wrk(kOffM),nVir(iSymb),One, &
                        P_ab(kP(iSymb)),nVir(iSymb))

          end do

        end if

      end do

    end if

  end do ! iSym

end if ! ChoAlg

end subroutine ChoMP2_fno_Fll
