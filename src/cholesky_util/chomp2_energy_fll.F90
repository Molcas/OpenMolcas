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

use ChoMP2, only: LiMatij

implicit real*8(a-h,o-z)
logical Delete
real*8 EOcc(*), EVir(*), Wrk(lWrk)
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
character*10 ThisNm
character*17 SecNam
parameter(SecNam='ChoMP2_Energy_Fll',ThisNm='Energy_Fll')
integer nEnrVec(8), LnT2am, LiT2am(8)
integer nVaJi, iVaJi(8)
real*8 X(0:1)
data X/0.0d0,1.0d0/
! Statement functions
MulD2h(i,j) = ieor(i-1,j-1)+1
iTri(i,j) = max(i,j)*(max(i,j)-3)/2+i+j

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
  call iCopy(nSym,nMP2Vec,1,nEnrVec,1)
else
  iTyp = 1
  call iCopy(nSym,NumCho,1,nEnrVec,1)
end if

! Set (ai|bj) indices.
! --------------------

call ChoMP2_Energy_GetInd(LnT2am,LiT2am,1,1)

! Allocate memory for integrals.
! ------------------------------

kXaibj = 1
kEnd0 = kXaibj+LnT2am
lWrk0 = lWrk-kEnd0+1
if (lWrk0 < 0) call ChoMP2_Quit(SecNam,'insufficient memory','[0]')

! Initialize MP2 energy correction.
! ---------------------------------

EMP2 = 0.0d0

! Special code for ChoAlg=2:
! compute M(ab,ij) = (ai|bj) with i<=j using level 3 BLAS.
! For ChoAlg=1: use strictly lower triangular storage (=>
! level 2 BLAS).
! --------------------------------------------------------

if (ChoAlg == 2) then ! level 3 BLAS algorithm

  kMabij = kXaibj ! rename pointer
  call FZero(Wrk(kMabij),LnT2am) ! initialize

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

      if (lWrk1 < Nai) call ChoMP2_Quit(SecNam,'Insufficient memory','[ChoAlg.2.1]')

      ! Set up batch over Cholesky vectors.
      ! -----------------------------------

      nVec = min(lWrk1/Nai,nEnrVec(iSym))
      if (nVec < 1) call ChoMP2_Quit(SecNam,'Insufficient memory','[ChoAlg.2.2]') ! should not happen
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
          iSyma = MulD2h(iSymi,iSym)
          iVaJi(iSymi) = nVaJi
          nVaJi = nVaJi+nVir(iSyma)*NumVec*nOcc(iSymi)
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
          lTot = nT1am(iSym)
          iAdr = nT1am(iSym)*(iVec0+iVec-1)+1
          call ddaFile(lUnit_F(iSym,iTyp),iOpt,Wrk(kVecai),lTot,iAdr)

          do iSymi=1,nSym
            iSyma = MulD2h(iSymi,iSym)
            do i=1,nOcc(iSymi)
              kOff1 = kVecai+iT1am(iSyma,iSymi)+nVir(iSyma)*(i-1)
              kOff2 = kVec+iVaJi(iSymi)+nVir(iSyma)*NumVec*(i-1)+nVir(iSyma)*(iVec-1)
              call dCopy_(nVir(iSyma),Wrk(kOff1),1,Wrk(kOff2),1)
            end do
          end do

        end do

        ! Compute M(ab,ij) for i<=j.
        ! First do iSymi=iSymj, then iSymi<iSymj.
        ! ---------------------------------------

        do iSymj=1,nSym

          iSymb = MulD2h(iSymj,iSym)

          if (nVir(iSymb) > 0) then

            do j=1,nOcc(iSymj)
              do i=1,j

                ij = LiMatij(iSymj,iSymj,1)+iTri(i,j)

                kOffi = kVec+iVaJi(iSymj)+nVir(iSymb)*NumVec*(i-1)
                kOffj = kVec+iVaJi(iSymj)+nVir(iSymb)*NumVec*(j-1)
                kOffM = kMabij+LiT2am(1)+nMatab(1)*(ij-1)+iMatab(iSymb,iSymb)

                call DGEMM_('N','T',nVir(iSymb),nVir(iSymb),NumVec,1.0d0,Wrk(kOffi),nVir(iSymb),Wrk(kOffj),nVir(iSymb),1.0d0, &
                            Wrk(kOffM),nVir(iSymb))

              end do
            end do

            do iSymi=1,iSymj-1

              iSyma = MulD2h(iSymi,iSym)
              iSymab = MulD2h(iSyma,iSymb)
              iSymij = iSymab

              if ((nOcc(iSymi) > 0) .and. (nVir(iSyma) > 0)) then

                do j=1,nOcc(iSymj)
                  do i=1,nOcc(iSymi)

                    ij = LiMatij(iSymi,iSymj,1)+nOcc(iSymi)*(j-1)+i

                    kOffi = kVec+iVaJi(iSymi)+nVir(iSyma)*NumVec*(i-1)
                    kOffj = kVec+iVaJi(iSymj)+nVir(iSymb)*NumVec*(j-1)
                    kOffM = kMabij+LiT2am(iSymij)+nMatab(iSymab)*(ij-1)+iMatab(iSyma,iSymb)

                    call DGEMM_('N','T',nVir(iSyma),nVir(iSymb),NumVec,1.0d0,Wrk(kOffi),nVir(iSyma),Wrk(kOffj),nVir(iSymb),1.0d0, &
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
      if (NumVec < 1) call ChoMP2_Quit(SecNam,'insufficient memory','[2]')
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
        call dGeMM_Tri('N','T',nT1am(iSym),nT1am(iSym),NumV,1.0d0,Wrk(kVec),nT1am(iSym),Wrk(kVec),nT1am(iSym),Fac,Wrk(kXint), &
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
