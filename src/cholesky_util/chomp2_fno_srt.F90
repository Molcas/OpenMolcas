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

subroutine ChoMP2_fno_Srt(irc,Delete,P_ab,P_ii,EOcc,EVir,Wrk,lWrk)
!
!  F. Aquilante, Geneva May 2008  (snick to Pedersen's code)

use ChoMP2, only: LnOcc, LnT1am, LiT1am, LiMatij, lUnit

implicit real*8(a-h,o-z)
logical Delete
real*8 EOcc(*), EVir(*), Wrk(lWrk), P_ab(*), P_ii(*)
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "chfnopt.fh"
character*10 ThisNm
character*17 SecNam
parameter(SecNam='ChoMP2_fno_Srt',ThisNm='fno_Srt')
integer nEnrVec(8), LnT2am, LiT2am(8), kP(8), lP(8)
integer nVaJi, iVaJi(8)
integer iDummy
parameter(iDummy=-999999)
real*8 xsDnom
! Statement functions
MulD2h(i,j) = ieor(i-1,j-1)+1
iTri(i,j) = max(i,j)*(max(i,j)-3)/2+i+j

irc = 0

kP(1) = 1
lP(1) = 0
do iS=2,nSym
  kP(iS) = kP(iS-1)+nVir(iS-1)**2
  lP(iS) = lP(iS-1)+nOcc(iS-1)
end do

! Set number of vectors.
! ----------------------

if (DecoMP2) then
  call iCopy(nSym,nMP2Vec,1,nEnrVec,1)
else
  call iCopy(nSym,NumCho,1,nEnrVec,1)
end if

if (MP2_small) then

  ! Loop over occupied orbital batches.
  ! -----------------------------------

  do iBatch=1,nBatch

    jBatch = iBatch

    call ChoMP2_Energy_GetInd(LnT2am,LiT2am,iBatch,jBatch)

    kXaibj = 1
    kEnd0 = kXaibj+LnT2am
    lWrk0 = lWrk-kEnd0+1
    if (lWrk0 < 1) call ChoMP2_Quit(SecNam,'insufficient memory','[0]')

    if (ChoAlg == 2) then

      kMabij = kXaibj  ! rename pointer
      call FZero(Wrk(kMabij),LnT2am) ! initialize

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

          if (lWrk1 < Nai) call ChoMP2_Quit(SecNam,'Insufficient memory','[ChoAlg.2.1]')

          ! Set up batch over Cholesky vectors.
          ! -----------------------------------

          nVec = min(lWrk1/Nai,nEnrVec(iSym))
          if (nVec < 1) call ChoMP2_Quit(SecNam,'Insufficient memory','[ChoAlg.2.2]') ! should not happen
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
              iSyma = MulD2h(iSymi,iSym)
              iVaJi(iSymi) = nVaJi
              nVaJi = nVaJi+nVir(iSyma)*NumVec*LnOcc(iSymi,iBatch)
            end do

            ! Pointer to reordered vectors: kVec.
            ! -----------------------------------

            kVec = kEnd1
            kEnd2 = kVec+nVaJi
            lWrk2 = lWrk-kEnd2+1
            if (lWrk2 < 0) call ChoMP2_Quit(SecNam,'Insufficient memory','[ChoAlg.2.3]') ! should not happen

            ! Read one vector at a time and reorder.
            ! --------------------------------------

            iVec0 = iVec1-1
            do iVec=1,NumVec

              iOpt = 2
              lTot = Nai
              iAdr = Nai*(iVec0+iVec-1)+1
              call ddaFile(lUnit(iSym,iBatch),iOpt,Wrk(kVecai),lTot,iAdr)

              do iSymi=1,nSym
                iSyma = MulD2h(iSymi,iSym)
                do i=1,LnOcc(iSymi,iBatch)
                  kOff1 = kVecai+LiT1am(iSyma,iSymi,iBatch)+nVir(iSyma)*(i-1)
                  kOff2 = kVec+iVaJi(iSymi)+nVir(iSyma)*NumVec*(i-1)+nVir(iSyma)*(iVec-1)
                  call dCopy_(nVir(iSyma),Wrk(kOff1),1,Wrk(kOff2),1)
                end do
              end do

            end do

            ! Compute M(ab,ii) .
            ! ------------------

            do iSymj=1,nSym

              iSymb = MulD2h(iSymj,iSym)

              if (nVir(iSymb) > 0) then

                do j=1,LnOcc(iSymj,iBatch)

                  i = j

                  ij = LiMatij(iSymj,iSymj,iBatch)+iTri(i,j)

                  kOffi = kVec+iVaJi(iSymj)+nVir(iSymb)*NumVec*(i-1)
                  kOffj = kVec+iVaJi(iSymj)+nVir(iSymb)*NumVec*(j-1)
                  kOffM = kMabij+LiT2am(1)+nMatab(1)*(ij-1)+iMatab(iSymb,iSymb)

                  call dGeMM_('N','T',nVir(iSymb),nVir(iSymb),NumVec,1.0d0,Wrk(kOffi),nVir(iSymb),Wrk(kOffj),nVir(iSymb),1.0d0, &
                              Wrk(kOffM),nVir(iSymb))

                end do

              end if

            end do

          end do

          call Cho_GAdGOp(Wrk(kMabij),LnT2am,'+')

          ! Close Cholesky vector files.
          ! ----------------------------

          call ChoMP2_OpenB(2,iSym,iBatch)

          do iSymj=1,nSym

            iSymb = MulD2h(iSymj,iSym)

            if (nVir(iSymb) > 0) then

              do j=1,LnOcc(iSymj,iBatch)

                i = j

                ij = LiMatij(iSymj,iSymj,iBatch)+iTri(i,j)

                kOffM = kMabij+LiT2am(1)+nMatab(1)*(ij-1)+iMatab(iSymb,iSymb)

                ! Compute energy contribution
                ! -------------------------------------
                do jb=1,nVir(iSymb)
                  do ja=1,nVir(iSymb)
                    Dnom = EVir(iVir(iSymb)+ja)+EVir(iVir(iSymb)+jb)-2.0d0*EOcc(iOcc(iSymj)+j)
                    xsDnom = Dnom/(Dnom**2+shf**2)
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

    end if

  end do ! iBatch

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

  return

end if

! Loop over occupied orbital batches.
! --------------------------------------------------

do iBatch=1,nBatch

  jBatch = iBatch

  call ChoMP2_Energy_GetInd(LnT2am,LiT2am,iBatch,jBatch)

  kXaibj = 1
  kEnd0 = kXaibj+LnT2am
  lWrk0 = lWrk-kEnd0+1
  if (lWrk0 < 1) call ChoMP2_Quit(SecNam,'insufficient memory','[0]')

  if (ChoAlg == 2) then

    kMabij = kXaibj  ! rename pointer
    call FZero(Wrk(kMabij),LnT2am) ! initialize

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

        if (lWrk1 < Nai) call ChoMP2_Quit(SecNam,'Insufficient memory','[ChoAlg.2.1]')

        ! Set up batch over Cholesky vectors.
        ! -----------------------------------

        nVec = min(lWrk1/Nai,nEnrVec(iSym))
        if (nVec < 1) call ChoMP2_Quit(SecNam,'Insufficient memory','[ChoAlg.2.2]') ! should not happen
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
            iSyma = MulD2h(iSymi,iSym)
            iVaJi(iSymi) = nVaJi
            nVaJi = nVaJi+nVir(iSyma)*NumVec*LnOcc(iSymi,iBatch)
          end do

          ! Pointer to reordered vectors: kVec.
          ! -----------------------------------

          kVec = kEnd1
          kEnd2 = kVec+nVaJi
          lWrk2 = lWrk-kEnd2+1
          if (lWrk2 < 0) call ChoMP2_Quit(SecNam,'Insufficient memory','[ChoAlg.2.3]') ! should not happen

          ! Read one vector at a time and reorder.
          ! --------------------------------------

          iVec0 = iVec1-1
          do iVec=1,NumVec

            iOpt = 2
            lTot = Nai
            iAdr = Nai*(iVec0+iVec-1)+1
            call ddaFile(lUnit(iSym,iBatch),iOpt,Wrk(kVecai),lTot,iAdr)

            do iSymi=1,nSym
              iSyma = MulD2h(iSymi,iSym)
              do i=1,LnOcc(iSymi,iBatch)
                kOff1 = kVecai+LiT1am(iSyma,iSymi,iBatch)+nVir(iSyma)*(i-1)
                kOff2 = kVec+iVaJi(iSymi)+nVir(iSyma)*NumVec*(i-1)+nVir(iSyma)*(iVec-1)
                call dCopy_(nVir(iSyma),Wrk(kOff1),1,Wrk(kOff2),1)
              end do
            end do

          end do

          ! Compute M(ab,ii) .
          ! ------------------

          do iSymj=1,nSym

            iSymb = MulD2h(iSymj,iSym)

            if (nVir(iSymb) > 0) then

              do j=1,LnOcc(iSymj,iBatch)

                i = j

                ij = LiMatij(iSymj,iSymj,iBatch)+iTri(i,j)

                kOffi = kVec+iVaJi(iSymj)+nVir(iSymb)*NumVec*(i-1)
                kOffj = kVec+iVaJi(iSymj)+nVir(iSymb)*NumVec*(j-1)
                kOffM = kMabij+LiT2am(1)+nMatab(1)*(ij-1)+iMatab(iSymb,iSymb)

                call dGeMM_('N','T',nVir(iSymb),nVir(iSymb),NumVec,1.0d0,Wrk(kOffi),nVir(iSymb),Wrk(kOffj),nVir(iSymb),1.0d0, &
                            Wrk(kOffM),nVir(iSymb))

              end do

            end if

          end do

          call Cho_GAdGOp(Wrk(kMabij),LnT2am,'+')

          ! Close Cholesky vector files.
          ! ----------------------------

          call ChoMP2_OpenB(2,iSym,iBatch)

          do iSymj=1,nSym

            iSymb = MulD2h(iSymj,iSym)

            if (nVir(iSymb) > 0) then

              do j=1,LnOcc(iSymj,iBatch)

                i = j

                ij = LiMatij(iSymj,iSymj,iBatch)+iTri(i,j)

                kOffM = kMabij+LiT2am(1)+nMatab(1)*(ij-1)+iMatab(iSymb,iSymb)

                ! Compute T(a,b)[i] and energy contrib.
                ! -------------------------------------
                do jb=1,nVir(iSymb)
                  do ja=1,nVir(iSymb)
                    Dnom = EVir(iVir(iSymb)+ja)+EVir(iVir(iSymb)+jb)-2.0d0*EOcc(iOcc(iSymj)+j)
                    xsDnom = Dnom/(Dnom**2+shf**2)
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
                call dGeMM_('N','N',nVir(iSymb),nVir(iSymb),nVir(iSymb),1.0d0,Wrk(kOffM),nVir(iSymb),Wrk(kOffM),nVir(iSymb),1.0d0, &
                            P_ab(kP(iSymb)),nVir(iSymb))

              end do

            end if

          end do

        end do

      end if

    end do ! iSym

  end if

end do ! iBatch

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

end subroutine ChoMP2_fno_Srt
