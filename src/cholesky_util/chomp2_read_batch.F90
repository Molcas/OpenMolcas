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
! Copyright (C) 2008, Jonas Bostrom                                    *
!***********************************************************************

subroutine ChoMP2_Read_Batch(LnPQRSprod,LiPQRSprod,Wrk,lWrk,iBatch,jBatch,kXpqrs)
!
! Jonas Bostrom, Aug 2008. (Generalization and modification of some
!                           stuff in ChoMP2_energy_*)
!
! Purpose: Reads parts of cholesky vectors from disk and
!          multiply them into two integral batches (or one batch if
!          jBatch=iBatch).

use Cholesky, only: nSym, NumCho
use ChoMP2, only: ChoAlg, LnPQprod, lUnit_F, nBatch, nPQ_prod
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LnPQRSprod, LiPQRSprod(8), lWrk, iBatch, jBatch
integer(kind=iwp), intent(out) :: kXpqrs
real(kind=wp), intent(out) :: Wrk(lWrk)
integer(kind=iwp) :: iAdr, iBat, iOpt, iSym, iTyp, iVec, iVec1, jVec, kEnd0, kEnd1, kEnd2, kOff, kRead, kVai, kVbj, kXint, lTot, &
                     lWrk0, lWrk1, lWrk2, MinMem, Nai, nBat, Nbj, nEnrVec(8), NumV, NumVec
real(kind=wp) :: Fac
real(kind=wp), parameter :: X(0:1) = [Zero,One]
character(len=*), parameter :: SecNam = 'ChoMP2_Read_Batch'

! Set number and type of vectors.
! -------------------------------

iTyp = 1
nEnrVec(1:nSym) = NumCho(1:nSym)

! Allocate memory for integrals.
! ------------------------------

kXpqrs = 1
kEnd0 = kXpqrs+LnPQRSprod
lWrk0 = lWrk-kEnd0+1
if (lWrk0 < 1) call SysAbendMsg(SecNam,'insufficient memory','[0]')

! Special code for iBatch=jBatch and ChoAlg=2:
! compute M(ab,ij) = (ai|bj) with i<=j using level 3 BLAS.
! For ChoAlg=1: use strictly lower triangular storage (=>
! level 2 BLAS).
! --------------------------------------------------------

if (ChoAlg == 2) then
  write(u6,*) 'No support for Cholesky algorithm 2'
else

  ! Loop over Cholesky vector symmetries.
  ! -------------------------------------

  do iSym=1,nSym

    Nai = LnPQprod(iSym,iBatch)
    Nbj = LnPQprod(iSym,jBatch)
    if ((Nai > 0) .and. (Nbj > 0) .and. (nEnrVec(iSym) > 0)) then

      ! Allocate memory for reading 1 vector.
      ! -------------------------------------

      kRead = kEnd0
      if (nBatch /= 1) then
        kEnd1 = kRead+nPQ_prod(iSym)
        lWrk1 = lWrk-kEnd1+1
        if (lWrk1 < 1) call SysAbendMsg(SecNam,'insufficient memory','[0.1]')
      else
        kEnd1 = kRead
        lWrk1 = lWrk0
      end if

      ! Setup Cholesky vector batching.
      ! -------------------------------

      if (jBatch == iBatch) then
        MinMem = Nai
      else
        MinMem = Nai+Nbj
      end if
      NumVec = min(lWrk1/MinMem,nEnrVec(iSym))
      if (NumVec < 1) call SysAbendMsg(SecNam,'insufficient memory','[1]')
      nBat = (nEnrVec(iSym)-1)/NumVec+1

      ! Open Cholesky vector file.
      ! --------------------------

      call ChoMP2_OpenF(1,iTyp,iSym)

      ! Cholesky vector batch loop.
      ! ---------------------------

      do iBat=1,nBat

        if (iBat == nBat) then
          NumV = nEnrVec(iSym)-NumVec*(nBat-1)
        else
          NumV = NumVec
        end if
        iVec1 = NumVec*(iBat-1)+1

        kVai = kEnd1
        kVbj = 0
        kEnd2 = 0
        lWrk2 = 0
        if (nBatch /= 1) then
          kVbj = kVai+Nai*NumV
          if (jBatch == iBatch) then
            kEnd2 = kVbj
            kVbj = kVai
          else
            kEnd2 = kVbj+Nbj*NumV
          end if
          lWrk2 = lWrk-kEnd2+1
          if (lWrk2 < 0) call SysAbendMsg(SecNam,'insufficient memory','[2]') ! this would be a bug...
        end if

        ! Read vectors, copy out sub-blocks.
        ! ----------------------------------

        if (nBatch == 1) then

          iOpt = 2
          lTot = nPQ_prod(iSym)*NumV
          iAdr = nPQ_prod(iSym)*(iVec1-1)+1
          call ddaFile(lUnit_F(iSym,iTyp),iOpt,Wrk(kVai),lTot,iAdr)
        else
          do iVec=1,NumV

            jVec = iVec1+iVec-1
            iOpt = 2
            lTot = nPQ_prod(iSym)
            iAdr = nPQ_prod(iSym)*(jVec-1)+1
            call ddaFile(lUnit_F(iSym,iTyp),iOpt,Wrk(kRead),lTot,iAdr)
            kOff = kVai+Nai*(iVec-1)
            call ChoMP2_Srt(Wrk(kRead),Wrk(kOff),1,iSym,iBatch)
            if (jBatch /= iBatch) then
              kOff = kVbj+Nbj*(iVec-1)
              call ChoMP2_Srt(Wrk(kRead),Wrk(kOff),1,iSym,jBatch)
            end if
          end do
        end if

        ! Compute integral contribution.
        ! ------------------------------

        Fac = X(min((iBat-1),1))
        kXint = kXpqrs+LiPQRSprod(iSym)
        if (iBatch == jBatch) then
          call dGeMM_Tri('N','T',Nai,Nai,NumV,One,Wrk(kVai),Nai,Wrk(kVai),Nai,Fac,Wrk(kXint),Nai)
        else
          call DGEMM_('N','T',Nai,Nbj,NumV,One,Wrk(kVai),Nai,Wrk(kVbj),Nbj,Fac,Wrk(kXint),Nai)
        end if
      end do           ! Cholesky vector batch

      ! Close Cholesky vector files.
      ! ----------------------------

      call ChoMP2_OpenF(2,iTyp,iSym)

    end if

  end do               ! iSym

end if

return

end subroutine ChoMP2_Read_Batch
