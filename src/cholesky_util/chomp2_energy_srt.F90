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
! Copyright (C) 2004,2005, Thomas Bondo Pedersen                       *
!***********************************************************************

subroutine ChoMP2_Energy_Srt(irc,Delete,EMP2,EOcc,EVir,Wrk,lWrk)
!
! Thomas Bondo Pedersen, Dec. 2004 / Feb. 2005.
!
! Purpose: compute MP2 energy contribution using presorted MO
!          Cholesky vectors on disk.

use Symmetry_Info, only: Mul
use Index_Functions, only: iTri
use Cholesky, only: nSym, NumCho
use ChoMP2, only: ChoAlg, DecoMP2, iMatab, LiMatij, LiT1am, LnOcc, LnT1am, lUnit, nBatch, nMatab, nMP2Vec, nVir, Verbose, Wref
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
logical(kind=iwp), intent(in) :: Delete
integer(kind=iwp), intent(in) :: lWrk
real(kind=wp), intent(out) :: EMP2, Wrk(lWrk)
real(kind=wp), intent(in) :: EOcc(*), EVir(*)
integer(kind=iwp) :: i, iAdr, iBat, iBatch, ij, iOpt, iSym, iSyma, iSymab, iSymb, iSymi, iSymij, iSymj, iVaJi(8), iVec, iVec0, &
                     iVec1, j, jBatch, kEnd0, kEnd1, kEnd2, kMabij, kOff1, kOff2, kOffi, kOffj, kOffM, kVai, kVbj, kVec, kVecai, &
                     kXaibj, kXint, LiT2am(8), LnT2am, lTot, lWrk0, lWrk1, lWrk2, MinMem, Nai, nBat, Nbj, nEnrVec(8), NumV, &
                     NumVec, nVaJi, nVec
real(kind=wp) :: Fac
integer(kind=iwp), parameter :: iDummy = -999999
real(kind=wp), parameter :: X(0:1) = [Zero,One]
character(len=*), parameter :: SecNam = 'ChoMP2_Energy_Srt'

irc = 0

! Set number of vectors.
! ----------------------

if (DecoMP2) then
  nEnrVec(1:nSym) = nMP2Vec(1:nSym)
else
  nEnrVec(1:nSym) = NumCho(1:nSym)
end if

! Initialize MP2 energy correction.
! ---------------------------------

EMP2 = Zero
Wref = Zero

! Print header of status table.
! -----------------------------

if (Verbose) call ChoMP2_Energy_Prt(SecNam,0,iDummy)

! Loop over occupied orbital batches.
! Special handling of diagonal batches for ChoAlg=2.
! --------------------------------------------------

do iBatch=1,nBatch
  if (Verbose) call ChoMP2_Energy_Prt(SecNam,1,iBatch)
  do jBatch=iBatch,nBatch

    call ChoMP2_Energy_GetInd(LnT2am,LiT2am,iBatch,jBatch)

    kXaibj = 1
    kEnd0 = kXaibj+LnT2am
    lWrk0 = lWrk-kEnd0+1
    if (lWrk0 < 1) call SysAbendMsg(SecNam,'insufficient memory','[0]')

    ! Special code for iBatch=jBatch and ChoAlg=2:
    ! compute M(ab,ij) = (ai|bj) with i<=j using level 3 BLAS.
    ! For ChoAlg=1: use strictly lower triangular storage (=>
    ! level 2 BLAS).
    ! --------------------------------------------------------

    if ((jBatch == iBatch) .and. (ChoAlg == 2)) then

      kMabij = kXaibj  ! rename pointer
      Wrk(kMabij:kMabij+LnT2am-1) = Zero ! initialize

      ! Loop over Cholesky vector symmetries.
      ! -------------------------------------

      do iSym=1,nSym

        Nai = LnT1am(iSym,iBatch)
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

          call ChoMP2_OpenB(1,iSym,iBatch)

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
              nVaJi = nVaJi+nVir(iSyma)*NumVec*LnOcc(iSymi,iBatch)
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
              lTot = Nai
              iAdr = Nai*(iVec0+iVec-1)+1
              call ddaFile(lUnit(iSym,iBatch),iOpt,Wrk(kVecai),lTot,iAdr)

              do iSymi=1,nSym
                iSyma = Mul(iSymi,iSym)
                do i=1,LnOcc(iSymi,iBatch)
                  kOff1 = kVecai+LiT1am(iSyma,iSymi,iBatch)+nVir(iSyma)*(i-1)
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

                do j=1,LnOcc(iSymj,iBatch)
                  do i=1,j

                    ij = LiMatij(iSymj,iSymj,iBatch)+iTri(i,j)

                    kOffi = kVec+iVaJi(iSymj)+nVir(iSymb)*NumVec*(i-1)
                    kOffj = kVec+iVaJi(iSymj)+nVir(iSymb)*NumVec*(j-1)
                    kOffM = kMabij+LiT2am(1)+nMatab(1)*(ij-1)+iMatab(iSymb,iSymb)

                    call DGEMM_('N','T',nVir(iSymb),nVir(iSymb),NumVec,One,Wrk(kOffi),nVir(iSymb),Wrk(kOffj),nVir(iSymb),One, &
                                Wrk(kOffM),nVir(iSymb))

                  end do
                end do

                do iSymi=1,iSymj-1

                  iSyma = Mul(iSymi,iSym)
                  iSymab = Mul(iSyma,iSymb)
                  iSymij = iSymab

                  if ((LnOcc(iSymi,iBatch) > 0) .and. (nVir(iSyma) > 0)) then

                    do j=1,LnOcc(iSymj,iBatch)
                      do i=1,LnOcc(iSymi,iBatch)

                        ij = LiMatij(iSymi,iSymj,iBatch)+LnOcc(iSymi,iBatch)*(j-1)+i

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

          ! Close Cholesky vector files.
          ! ----------------------------

          call ChoMP2_OpenB(2,iSym,iBatch)

        end if

      end do ! iSym

    else ! level 2 BLAS for diagonal batches.

      ! Loop over Cholesky vector symmetries.
      ! -------------------------------------

      do iSym=1,nSym

        Nai = LnT1am(iSym,iBatch)
        Nbj = LnT1am(iSym,jBatch)
        if ((Nai > 0) .and. (Nbj > 0) .and. (nEnrVec(iSym) > 0)) then

          ! Setup Cholesky vector batching.
          ! -------------------------------

          if (jBatch == iBatch) then
            MinMem = Nai
          else
            MinMem = Nai+Nbj
          end if
          NumVec = min(lWrk0/MinMem,nEnrVec(iSym))
          if (NumVec < 1) call SysAbendMsg(SecNam,'insufficient memory','[1]')
          nBat = (nEnrVec(iSym)-1)/NumVec+1

          ! Open Cholesky vector files.
          ! ---------------------------

          call ChoMP2_OpenB(1,iSym,iBatch)
          if (jBatch /= iBatch) call ChoMP2_OpenB(1,iSym,jBatch)

          ! Cholesky vector batch loop.
          ! ---------------------------

          do iBat=1,nBat

            if (iBat == nBat) then
              NumV = nEnrVec(iSym)-NumVec*(nBat-1)
            else
              NumV = NumVec
            end if
            iVec1 = NumVec*(iBat-1)+1

            kVai = kEnd0
            kVbj = kVai+Nai*NumV
            if (jBatch == iBatch) then
              kEnd1 = kVbj
              kVbj = kVai
            else
              kEnd1 = kVbj+Nbj*NumV
            end if
            lWrk1 = lWrk-kEnd1+1
            if (lWrk1 < 0) call SysAbendMsg(SecNam,'insufficient memory','[2]') ! this would be a bug...

            ! Read vectors.
            ! -------------

            iOpt = 2
            lTot = Nai*NumV
            iAdr = Nai*(iVec1-1)+1
            call ddaFile(lUnit(iSym,iBatch),iOpt,Wrk(kVai),lTot,iAdr)
            if (jBatch /= iBatch) then
              iOpt = 2
              lTot = Nbj*NumV
              iAdr = Nbj*(iVec1-1)+1
              call ddaFile(lUnit(iSym,jBatch),iOpt,Wrk(kVbj),lTot,iAdr)
            end if

            ! Compute integral contribution.
            ! ------------------------------

            Fac = X(min((iBat-1),1))
            kXint = kXaibj+LiT2am(iSym)
            if (iBatch == jBatch) then
              call dGeMM_Tri('N','T',Nai,Nai,NumV,One,Wrk(kVai),Nai,Wrk(kVai),Nai,Fac,Wrk(kXint),Nai)
            else
              call DGEMM_('N','T',Nai,Nbj,NumV,One,Wrk(kVai),Nai,Wrk(kVbj),Nbj,Fac,Wrk(kXint),Nai)
            end if

          end do ! Cholesky vector batch

          ! Close Cholesky vector files.
          ! ----------------------------

          call ChoMP2_OpenB(2,iSym,iBatch)
          if (jBatch /= iBatch) call ChoMP2_OpenB(2,iSym,jBatch)

        end if

      end do ! iSym

    end if

    ! Compute contribution to MP2 energy correction.
    ! ----------------------------------------------

    call ChoMP2_Energy_Contr(EMP2,EOcc,EVir,Wrk(kXaibj),LnT2am,LiT2am,iBatch,jBatch)

  end do ! jBatch
  if (Verbose) call ChoMP2_Energy_Prt(SecNam,2,iBatch)
end do ! iBatch

! Finish table.
! -------------

if (Verbose) call ChoMP2_Energy_Prt(SecNam,3,iBatch)

! Delete files if requested.
! --------------------------

if (Delete) then
  do iBatch=1,nBatch
    do iSym=1,nSym
      call ChoMP2_OpenB(1,iSym,iBatch)
      call ChoMP2_OpenB(3,iSym,iBatch)
    end do
  end do
end if

! Change sign on energy.
! ----------------------

EMP2 = -EMP2

end subroutine ChoMP2_Energy_Srt
