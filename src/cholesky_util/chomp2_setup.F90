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

subroutine ChoMP2_Setup(irc)
!
! Thomas Bondo Pedersen, Oct. 2004 / Feb. 2005.
!
! Purpose: setup of Cholesky MP2 program.

use ChoMP2, only: ChoMP2_allocated, iFirst, iFirstS, NumOcc
use ChoMP2, only: LnOcc, LnT1am, LiT1am, LnMatij, LiMatij
use ChoMP2, only: lUnit, NumBatOrb, LnBatOrb
use ChoMP2, only: LnPQprod, LiPQprod
use stdalloc

implicit real*8(a-h,o-z)
#include "cholesky.fh"
#include "choorb.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
character*12 SecNam
parameter(SecNam='ChoMP2_Setup')
character*2 Unt
integer NumVec(8), nFrac(8)
logical Accepted
integer bsize, blast
real*8 lX
logical ChoMP2_Setup_MemChk
! Statement function
MulD2h(i,j) = ieor(i-1,j-1)+1

irc = 0

! Setup index arrays and counters.
! --------------------------------

if ((DecoMP2 .or. DoDens) .and. (ThrMP2 <= 0.0d0)) call Get_dScalar('Cholesky Threshold',ThrMP2)

call ChoMP2_GetInf(nOrb,nOcc,nFro,nDel,nVir)
iOcc(1) = 0
iBatOrb(1) = 0
iVir(1) = 0
iFro(1) = 0
iDel(1) = 0
nOccT = nOcc(1)
nVirT = nVir(1)
nFroT = nFro(1)
nDelT = nDel(1)
if (.false.) then
  nBatOrbT = nVir(1)+nOcc(1)+nFro(1)+nDel(1)
else
  nBatOrbT = nOcc(1)
end if
do iSym=2,nSym
  iOcc(iSym) = nOccT
  iVir(iSym) = nVirT
  iFro(iSym) = nFroT
  iDel(iSym) = nDelT
  iBatOrb(iSym) = nBatOrbT
  nOccT = nOccT+nOcc(iSym)
  nVirT = nVirT+nVir(iSym)
  nFroT = nFroT+nFro(iSym)
  nDelT = nDelT+nDel(iSym)
  if (.false.) then
    nBatOrbT = nBatOrbT+nOcc(iSym)+nVir(iSym)+nFro(iSym)+nDel(iSym)
  else
    nBatOrbT = nBatOrbT+nOcc(iSym)
  end if
end do

do iSym=1,nSym
  nT1am(iSym) = 0
  do iSymi=1,nSym
    iSyma = MulD2h(iSymi,iSym)
    iT1am(iSyma,iSymi) = nT1am(iSym)
    nT1am(iSym) = nT1am(iSym)+nVir(iSyma)*nOcc(iSymi)
  end do
end do

if (.false.) then
  do iSym=1,nSym
    nPQ_prod(iSym) = 0
    nPQ_prodij(iSym) = 0
    nPQ_prodia(iSym) = 0
    nPQ_prodab(iSym) = 0
    do iSymQ=1,nSym
      iSymP = MulD2h(iSymQ,iSym)
      iPQ_prod(iSymP,iSymQ) = nPQ_prod(iSym)

      nPQ_prod(iSym) = nPQ_Prod(iSym)+(nOcc(iSymP)+nVir(iSymP)+nFro(iSymP)+nDel(iSymP))*(nOcc(iSymQ)+nVir(iSymQ)+nFro(iSymQ)+ &
                       nDel(iSymQ))
      nPQ_prodij(iSym) = nPQ_Prodij(iSym)+(nFro(iSymP)+nOcc(iSymP))*(nFro(iSymQ)+nOcc(iSymQ))
      nPQ_prodia(iSym) = nPQ_Prodia(iSym)+(nFro(iSymP)+nOcc(iSymP))*(nVir(iSymQ)+nDel(iSymQ))
      nPQ_prodab(iSym) = nPQ_Prodab(iSym)+(nVir(iSymP)+nDel(iSymP))*(nVir(iSymQ)+nDel(iSymQ))
    end do
  end do

end if

do iSym=1,nSym
  nT1AOT(iSym) = 0
  do iSymAl=1,nSym
    iSymi = MulD2h(iSymAl,iSym)
    iT1AOT(iSymi,iSymAl) = nT1AOT(iSym)
    nT1AOT(iSym) = nT1AOT(iSym)+nOcc(iSymi)*nBas(iSymAl)
  end do
end do

do iSym=1,nSym
  nAOVir(iSym) = 0
  do iSyma=1,nSym
    iSymAl = MulD2h(iSyma,iSym)
    iAOVir(iSymAl,iSyma) = nAOVir(iSym)
    nAOVir(iSym) = nAOVir(iSym)+nBas(iSymAl)*nVir(iSyma)
  end do
end do

if (ChoAlg == 2) then
  do iSym=1,nSym
    nMatab(iSym) = 0
    do iSymb=1,nSym
      iSyma = MulD2h(iSymb,iSym)
      iMatab(iSyma,iSymb) = nMatab(iSym)
      nMatab(iSym) = nMatab(iSym)+nVir(iSyma)*nVir(iSymb)
    end do
  end do
else
  call iZero(nMatab,8)
  call iZero(iMatab,64)
end if

! If batching over occuped orbitals is forced by user, there better
! be orbitals to batch over!
! -----------------------------------------------------------------

if (ForceBatch .and. (nBatOrbT == 1)) ForceBatch = .false.

! Setup batches over occupied orbitals.
! -------------------------------------

! nPQ_prodx will be the largest number of products in one symmetry
nPQProdx = 0
if (.false.) then
  nPQProdx = nPQ_Prod(1)
  do iSym=2,nSym
    nPQProdx = max(nPQProdx,nPQ_Prod(iSym))
  end do
end if

nT1amx = nT1am(1)
do iSym=2,nSym
  nT1amx = max(nT1amx,nT1am(iSym))
end do

if (DecoMP2) then
  do iSym=1,nSym
    NumVec(iSym) = min(NumCho(iSym),nT1am(iSym))
  end do
else
  call iCopy(nSym,NumCho,1,NumVec,1)
end if
call Cho_GAiGOp(NumVec,nSym,'max')

! nBatOrbT is the total numbers of orbitals to batch over,
! if densities are not to be computed this is set equal to nOccT

if (nBatOrbT < 6) then
  mBatch = nBatOrbT
else
  mBatch = nBatOrbT/2+1
end if
Accepted = .false.
nBatch = 0
do while ((nBatch < mBatch) .and. (.not. Accepted))

  nBatch = nBatch+1
  if (nBatch == mBatch) then
    nBatch = nBatOrbT
    do iSym=1,nSym
      nFrac(iSym) = max(NumVec(iSym),1)
    end do
  else
    do iSym=1,nSym
      nFrac(iSym) = max(min(NumVec(iSym),20),1)
    end do
  end if

  call ChoMP2_deallocate(irc)
  ChoMP2_allocated = .true.

  call mma_allocate(iFirst,nBatch,Label='iFirst')
  call mma_allocate(iFirstS,nSym,nBatch,Label='iFirstS')
  call mma_allocate(NumOcc,nBatch,Label='NumOcc')
  call mma_allocate(LnOcc,nSym,nBatch,Label='LnOcc')
  call mma_allocate(LnT1am,nSym,nBatch,Label='LnT1am')
  call mma_allocate(LiT1am,nSym,nSym,nBatch,Label='LiT1am')
  call mma_allocate(lUnit,nSym,nBatch,Label='lUnit')

  if (ChoAlg == 2) then
    call mma_allocate(LnMatij,nSym,nBatch,Label='LnMatij')
    call mma_allocate(LiMatij,nSym,nSym,nBatch,Label='LiMatij')
  else
    call mma_allocate(LnMatij,1,1,Label='LnMatij')
    call mma_allocate(LiMatij,1,1,1,Label='LiMatij')
  end if

  ! Generalization of NumOcc for arbitrary quantity to batch over
  ! Would be good to kill NumOcc safely and only use one...
  call mma_allocate(NumBatOrb,nBatch,Label='NumBatOrb')
  ! Generalization of LnOcc for arbitrary quantity to batch over.
  call mma_allocate(LnBatOrb,nSym,nBatch,Label='LnBatOrb')
  if (.false.) then
    call mma_allocate(LnPQprod,nSym,nBatch,Label='LnPQprod')
    call mma_allocate(LiPQprod,nSym,nSym,nBatch,Label='LiPQprod')
  else
    call mma_allocate(LnPQprod,1,1,Label='LnPQprod')
    call mma_allocate(LiPQprod,1,1,1,Label='LiPQprod')
  end if

  call ChoMP2_Setup_Index(iFirst,iFirstS,NumOcc,LnOcc,NumBatOrb,LnBatOrb,LnT1am,LiT1am,LnPQprod,LiPQprod,LnMatij,LiMatij,nSym, &
                          nBatch)

  call mma_maxDBLE(lWork)
  if (.false.) then
    ! All Memory available minus one full vector and some small
    ! vectors for the PCG-algorithm.
    lAvail = lWork-nPQprodx-l_Mp2Lagr*9
  else if (Laplace .and. SOS_MP2) then
    lX = 0.0d0
    do iSym=1,nSym
      if ((nT1am(iSym) > 0) .and. (NumVec(iSym) > 0)) then
        bsize = min(Laplace_BlockSize,NumVec(iSym))
        nBlock = (NumVec(iSym)-1)/bsize+1
        blast = NumVec(iSym)-bsize*(nBlock-1)
        xM = dble(NumVec(iSym))
        xn = dble(nBlock)
        xb = dble(bsize)
        xbp = dble(blast)
        lX = max(lX,0.5d0*(xM*(xM+1.0d0)+(xn-1.0d0)*xb*(xb-1.0d0)+xbp*(xbp-1.0d0)))
      end if
    end do
    l_X = int(lX)
    if (l_X < 0) then
      write(Lupri,'(A,A)') SecNam,': dimension of X matrix is negative!'
      write(Lupri,'(A,I15)') 'l_X=',l_X
      if (lX > 0.0d0) then
        write(LuPri,'(A)') 'This seems to be an integer overflow!'
        call Cho_RWord2Byte(lX,Byte,Unt)
        write(LuPri,'(A,1P,D15.6,A,D15.6,1X,A,A)') 'In double precision, lX=',lX,' words (',Byte,Unt,')'
      end if
      irc = 1
      lAvail = 0
    else
      lAvail = lWork-l_X
    end if
  else
    lAvail = lWork-nT1amx ! all minus one vector (for reading)
  end if
  call GAiGOp_Scal(lAvail,'min')
  ! The argument LnPQprod is only used for the case where full
  ! Lpq-vectors are transformed for densities. Will be a dummy arg
  ! for regular MP2.
  Accepted = ChoMP2_Setup_MemChk(LnT1am,LnPQprod,NumVec,nFrac,nSym,nBatch,lAvail)

  if (ForceBatch .and. (nBatch == 1)) Accepted = .false. ! force batching

  if (.not. Accepted) call ChoMP2_deallocate(irc)

  if (irc == 1) then
    Accepted = .false.
    nBatch = mBatch ! break while loop
  end if

end do

if (.not. Accepted) then
  write(6,*) 'Accepted = false'
  irc = -1
  return
end if

! Initialize file units.
! ----------------------

do iSym=1,nSym
  do iTyp=1,nTypF
    call ChoMP2_OpenF(0,iTyp,iSym)
  end do
end do

end subroutine ChoMP2_Setup
