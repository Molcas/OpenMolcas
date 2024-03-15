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
! Copyright (C) 2008, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_CheckBackTra(iTyp,COcc,CVir,lU_AO)
!
! Thomas Bondo Pedersen, Jan. 2008.
!
! Purpose: check backtransformation of vectors (MO->AO).
!          A summary is printed at the end of this routine.
!
! The check is simple, comparing the quantities
!
!    X(J) = sum_alpha,beta L(J;alpha,beta)
!    Y(J) = sum_ai L(ai,J)*P(a)*P(i)
!
! where
!
!    P(i) = sum_alpha COcc(i,alpha)
!    P(a) = sum_alpha CVir(alpha,a)
!
! The reported abs. min., abs. max, average, and RMS errors
! are calculated from the vector
!
!    D(J) = X(J) - Y(J)

use Symmetry_Info, only: Mul
use Cholesky, only: nBas, nSym
use ChoMP2, only: iAOVir, iOcc, iT1am, iT1AOT, iVir, lUnit_F, nMP2Vec, nOcc, nOccT, nT1am, nVir, nVirT
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iTyp, lU_AO(*)
real(kind=wp), intent(in) :: COcc(*), CVir(*)
integer(kind=iwp) :: a, Al, AlBe, i, iAdr, iOpt, iSym, iSyma, iSymAl, iSymBe, iSymi, J, kC, kP, kP_, lTot, na, nAlBe, nMP2Vec_Tot, &
                     nOcc_Max
real(kind=wp) :: AbsMaxErr, AbsMinErr, AvgErr, Err(4,8), RMSErr
real(kind=wp), allocatable :: POcc(:), PVir(:), Q(:), X(:), Y(:), D(:), V(:)
real(kind=wp), external :: ddot_

! Initializations.
! ----------------

nMP2Vec_Tot = 0

AvgErr = Zero
RMSErr = Zero

nOcc_Max = nOcc(1)
do iSym=2,nSym
  nOcc_Max = max(nOcc_Max,nOcc(iSym))
end do

call mma_allocate(POcc,nOccT,Label='POcc')
call mma_allocate(PVir,nVirT,Label='PVir')
call mma_allocate(Q,nOcc_Max,Label='Q')

! P(i) = sum_alpha COcc(i,alpha)
! ------------------------------

POcc(:) = Zero
do iSymi=1,nSym
  iSymAl = iSymi
  if ((nOcc(iSymi) > 0) .and. (nBas(iSymAl) > 0)) then
    kP = iOcc(iSymi)
    do Al=1,nBas(iSymAl)
      kC = iT1AOT(iSymi,iSymAl)+nOcc(iSymi)*(Al-1)
      POcc(kP+1:kP+nOcc(iSymi)) = POcc(kP+1:kP+nOcc(iSymi))+COcc(kC+1:kC+nOcc(iSymi))
    end do
  end if
end do

! P(a) = sum_alpha CVir(alpha,a)
! ------------------------------

PVir(:) = Zero
do iSyma=1,nSym
  iSymAl = iSyma
  if ((nVir(iSyma) > 0) .and. (nBas(iSymAl) > 0)) then
    kP_ = iVir(iSyma)
    do a=1,nVir(iSyma)
      kP = kP_+a
      kC = iAOVir(iSymAl,iSyma)+nBas(iSymAl)*(a-1)
      PVir(kP) = PVir(kP)+sum(CVir(kC+1:kC+nBas(iSyma)))
    end do
  end if
end do

! Check each symmetry block.
! --------------------------

do iSym=1,nSym

  if (nMP2Vec(iSym) > 0) then

    ! Allocation.
    ! -----------

    call mma_allocate(X,nMP2Vec(iSym),Label='X')
    call mma_allocate(Y,nMP2Vec(iSym),Label='Y')

    ! Zero result arrays.
    ! -------------------

    X(:) = Zero
    Y(:) = Zero

    ! X(J) = sum_alpha,beta L(J;alpha,beta)
    ! -------------------------------------

    nAlBe = 0
    do iSymBe=1,nSym
      iSymAl = Mul(iSymBe,iSym)
      nAlBe = nAlBe+nBas(iSymAl)*nBas(iSymBe)
    end do

    call mma_allocate(V,nMP2Vec(iSym),Label='V')
    do AlBe=1,nAlBe
      iOpt = 2
      lTot = nMP2Vec(iSym)
      iAdr = nMP2Vec(iSym)*(AlBe-1)+1
      call ddaFile(lU_AO(iSym),iOpt,V,lTot,iAdr)
      X(:) = X(:)+V(:)
    end do
    call mma_deallocate(V)

    ! Y(J) = sum_ai L(ai,J)*P(a)*P(i)
    ! -------------------------------

    call mma_allocate(V,nT1Am(iSym),Label='V')
    iOpt = 1
    call ChoMP2_OpenF(iOpt,iTyp,iSym)
    do J=1,nMP2Vec(iSym)
      iOpt = 2
      lTot = nT1Am(iSym)
      iAdr = nT1Am(iSym)*(J-1)+1
      call ddaFile(lUnit_F(iSym,iTyp),iOpt,V,lTot,iAdr)
      do iSymi=1,nSym
        iSyma = Mul(iSymi,iSym)
        na = max(nVir(iSyma),1)
        call dGeMV_('T',nVir(iSyma),nOcc(iSymi),One,V(1+iT1Am(iSyma,iSymi)),na,PVir(1+iVir(iSyma)),1,Zero,Q,1)
        Y(J) = Y(J)+dDot_(nOcc(iSymi),Q,1,POcc(1+iOcc(iSymi)),1)
      end do
    end do
    iOpt = 2
    call ChoMP2_OpenF(iOpt,iTyp,iSym)
    call mma_deallocate(V)

    ! Calculate errors.
    ! -----------------

    call mma_allocate(D,nMP2Vec(iSym),Label='D')
    D(:) = X(:)-Y(:)
    Err(1,iSym) = abs(D(1))
    Err(2,iSym) = abs(D(1))
    Err(3,iSym) = D(1)
    do J=1,nMP2Vec(iSym)-1
      Err(1,iSym) = min(Err(1,iSym),abs(D(1+J)))
      Err(2,iSym) = max(Err(2,iSym),abs(D(1+J)))
      Err(3,iSym) = Err(3,iSym)+D(1+J)
    end do
    AvgErr = AvgErr+Err(3,iSym)
    Err(3,iSym) = Err(3,iSym)/real(nMP2Vec(iSym),kind=wp)
    Err(4,iSym) = dDot_(nMP2Vec(iSym),D,1,D,1)
    RMSErr = RMSErr+Err(4,iSym)
    Err(4,iSym) = Err(4,iSym)/real(nMP2Vec(iSym),kind=wp)
    Err(4,iSym) = sqrt(Err(4,iSym))
    call mma_deallocate(D)

    ! Deallocation.
    ! -------------

    call mma_deallocate(Y)
    call mma_deallocate(X)

  else

    Err(:,iSym) = Zero

  end if

  if (iSym == 1) then
    AbsMinErr = Err(1,iSym)
    AbsMaxErr = Err(2,iSym)
    nMP2Vec_Tot = max(nMP2Vec(iSym),0)
  else
    AbsMinErr = min(AbsMinErr,Err(1,iSym))
    AbsMaxErr = max(AbsMaxErr,Err(2,iSym))
    nMP2Vec_Tot = nMP2Vec_Tot+max(nMP2Vec(iSym),0)
  end if

end do

AvgErr = AvgErr/real(nMP2Vec_Tot,kind=wp)
RMSErr = RMSErr/real(nMP2Vec_Tot,kind=wp)
RMSErr = sqrt(RMSErr)

! Report.
! -------

call Cho_Head('MO Vector Backtransformation Check','=',80,u6)
write(u6,'(/,2X,A,/,2X,A)') 'Symmetry  Min. Abs. Error  Max. Abs. Error    Average Error       RMS Error', &
                            '----------------------------------------------------------------------------'
do iSym=1,nSym
  write(u6,'(4X,I2,4X,4(3X,ES14.6))') iSym,(Err(i,iSym),i=1,4)
end do
write(u6,'(2X,A)') '----------------------------------------------------------------------------'
write(u6,'(2X,A,2X,4(3X,ES14.6))') 'Total:',AbsMinErr,AbsMaxErr,AvgErr,RMSErr
write(u6,'(2X,A,/)') '----------------------------------------------------------------------------'

! Free memory.
! ------------

call mma_deallocate(Q)
call mma_deallocate(PVir)
call mma_deallocate(POcc)

end subroutine ChoMP2_CheckBackTra
