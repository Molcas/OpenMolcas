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
! Copyright (C) 2007, Francesco Aquilante                              *
!***********************************************************************

subroutine Cho_SOSmp2_Energy(irc,EMP2,EOcc,EVir,Delete)
! Francesco Aquilante, May 2007.
!
! Purpose: compute "Scaled Opposite-Spin" MP2 energy correction from
!          MO Cholesky vectors of the matrix M(ai,bj)=(ai|bj)^2.

use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(out) :: EMP2
real(kind=wp), intent(in) :: EOcc(*), EVir(*)
logical(kind=iwp), intent(in) :: Delete
integer(kind=iwp) :: ia, iAdr, iaiS, iaS, iaSoff(8), iaSym, iaT, iBat, ii, iiSoff(8), iiSym, iiT, iSym, iTyp, jSym, jVec, kRead, &
                     lTot, lWrk, MaxNVec, Nai, nAt, nBat, nEnrVec(8), nIt, NKVec, nOV, NumVec, nVec
real(kind=wp) :: Dmax
integer(kind=iwp), allocatable :: iD_bj(:)
real(kind=wp), allocatable :: W(:), Wrk(:), Y(:,:)
character(len=*), parameter :: SecNam = 'Cho_SOSmp2_Energy'
real(kind=wp), external :: ddot_
#include "cholesky.fh"
#include "chomp2.fh"
#include "chomp2_cfg.fh"

irc = 0

iTyp = 2
nEnrVec(:) = nMP2Vec

! Initialize SOS-MP2 energy correction.
! -------------------------------------

EMP2 = Zero

! Some offsets
! ------------
nIt = nOcc(1)
nAt = nVir(1)
iiSoff(1) = 0
iaSoff(1) = 0
do iSym=2,nSym
  iiSoff(iSym) = nIt
  iaSoff(iSym) = nAt
  nIt = nIt+nOcc(iSym)
  nAt = nAt+nVir(iSym)
end do
MaxNVec = nIt*nAt

call mma_allocate(W,MaxNVec,label='Ea-Ei')
call mma_allocate(iD_bj,MaxNVec,label='iD_bj')

! Loop over Cholesky vector symmetries.
! -------------------------------------

do jSym=1,nSym

  Nai = nT1am(jSym)
  if ((Nai > 0) .and. (nEnrVec(jSym) > 0)) then

    nOV = 0
    do iiSym=1,nSym
      iaSym = Mul(iiSym,jSym)
      do ii=1,nOcc(iiSym)
        iiT = iiSoff(iiSym)+ii
        iaS = nOV+nVir(iaSym)*(ii-1)
        do ia=1,nVir(iaSym)
          iaT = iaSoff(iaSym)+ia
          iaiS = iaS+ia
          W(iaiS) = EVir(iaT)-EOcc(iiT)
        end do
      end do
      nOV = nOV+nVir(iaSym)*nOcc(iiSym) ! ... = Nai
    end do

    ! Cholesky decompsition of the Orbital Energy
    ! Denominators (OED)
    ! -------------------------------------------
    call CHO_GET_ORD_bj(nOV,MaxNVec,OED_Thr,W,iD_bj,NKVec,Dmax)

    if (Verbose .or. (NKVec < 1)) then
      write(u6,'(A)') '---------------------------------------'
      write(u6,'(A,I2,A)') 'Orbital energy denominators CD (sym=',jSym,')'
      write(u6,'(A)') '---------------------------------------'
      write(u6,'(1X,A,I3,A,I9,A,1P,D25.16)') 'Number of vectors needed: ',NKVec,'   ( nAocc x nAvir : ',nOV,' ), max residual:',Dmax
      call xFlush(u6)
    end if

    if (NKVec > 0) then

      call mma_allocate(Y,nOV,NKVec,label='Yai_k')
      ! init to one the 1st col
      Y(:,1) = One

      call CHO_GET_OED_cd(.true.,nOV,W,NKVec,iD_bj,1,Y,Y)

      call mma_maxDBLE(lWrk)
      call mma_allocate(Wrk,lWrk,label='GetMax')

      ! Set up batch over Cholesky vectors.
      ! -----------------------------------

      nVec = min(lWrk/(Nai+NKVec),nEnrVec(jSym))
      if (nVec < 1) then
        call ChoMP2_Quit(SecNam,'Insufficient memory','Batch setup')
      end if
      nBat = (nEnrVec(jSym)-1)/nVec+1

      kRead = NKVec*nVec+1

      ! Open Cholesky vector files.
      ! ---------------------------

      call ChoMP2_OpenF(1,iTyp,jSym)

      ! Start vector batch loop.
      ! ------------------------

      do iBat=1,nBat

        if (iBat == nBat) then
          NumVec = nEnrVec(jSym)-nVec*(nBat-1)
        else
          NumVec = nVec
        end if
        jVec = nVec*(iBat-1)+1

        ! Read vectors
        ! ------------
        lTot = nT1am(jSym)*NumVec
        iAdr = nT1am(jSym)*(jVec-1)+1
        call ddaFile(lUnit_F(jSym,iTyp),2,Wrk(kRead),lTot,iAdr)

        ! Compute   E(k,J) = sum_ai Y(ai,k) * R(ai,J)
        ! --------------------------------------------
        call DGEMM_('T','N',NKVec,NumVec,Nai,One,Y,Nai,Wrk(kRead),Nai,Zero,Wrk,NKVec)

        ! Compute (unscaled) SOS-MP2 energy
        ! -----------------------------------

        EMP2 = EMP2+ddot_(NKVec*NumVec,Wrk,1,Wrk,1)

      end do ! Cholesky vector batch

      ! Close Cholesky vector files.
      ! ----------------------------
      call ChoMP2_OpenF(2,iTyp,jSym)

      call mma_deallocate(Wrk)
      call mma_deallocate(Y)

    end if

  end if

end do

call mma_deallocate(W)
call mma_deallocate(iD_bj)

! If requested, delete vector files.
! ----------------------------------

if (Delete) then
  do jSym=1,nSym
    call ChoMP2_OpenF(1,iTyp,jSym)
    call ChoMP2_OpenF(3,iTyp,jSym)
  end do
end if

! Change sign and use proper factor on energy.
! --------------------------------------------

!-tbp, December 2012: removed factor 2
!-tbp  (wrong result for two-electron systems with factor 2)
!-tbp EMP2 = -Two*EMP2
EMP2 = -EMP2

return

end subroutine Cho_SOSmp2_Energy
