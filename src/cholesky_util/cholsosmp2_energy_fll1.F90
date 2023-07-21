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
! Copyright (C) 2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoLSOSMP2_Energy_Fll1(N,w,t,EOcc,EVir,Delete,EMP2,irc)
!
! Thomas Bondo Pedersen, December 2012.
!
! Compute Laplace-SOS-MP2 energy correction from full Cholesky
! vectors (i.e., not batched), reading the vectors only once
! at the expense of memory.

use stdalloc

implicit none
integer N
real*8 w(N)
real*8 t(N)
real*8 EOcc(*)
real*8 EVir(*)
logical Delete
real*8 EMP2
integer irc
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "cholesky.fh"
character(len=22), parameter :: SecNam = 'ChoLSOSMP2_Energy_Fll1'
real*8, external :: dDot_
integer nEnrVec(8)
integer iClos
integer iTyp
integer iSym
integer l_X
integer Nai
integer nBlock
integer iOpt, iAddr, l_Tot
integer q
integer iBlock, jBlock
integer iSymi, iSyma
integer i, a
integer ipi, ipj
integer nVeci, nVecj
integer iVec
integer ip0, ip1
real*8 tq, Eq
real*8, allocatable :: X(:), V(:,:)
integer j, k
real*8 epsi, epsa
integer MulD2h
! Statement functions
epsi(j,k) = EOcc(iOcc(k)+j)
epsa(j,k) = EVir(iVir(k)+j)
MulD2h(j,k) = ieor(j-1,k-1)+1

! init return code
irc = 0

! init energy
EMP2 = 0.0d0

! check input (incl. common block variables)
if (nBatch /= 1) then
  irc = -1
  return
end if
if (N /= Laplace_nGridPoints) then
  irc = -2
  return
end if
if (Laplace_BlockSize < 1) then
  irc = -3
  return
end if

! determine if files are to be deleted after use
if (Delete) then
  iClos = 3
else
  iClos = 2
end if

! set number and type of vectors
if (DecoMP2) then
  iTyp = 2
  call iCopy(nSym,nMP2Vec,1,nEnrVec,1)
else
  iTyp = 1
  call iCopy(nSym,NumCho,1,nEnrVec,1)
end if

! compute energy correction
do iSym=1,nSym
  Nai = nT1am(iSym)
  if ((Nai > 0) .and. (nEnrVec(iSym) > 0)) then
    ! set vector block size and allocate X matrix
    l_X = min(Laplace_BlockSize,nEnrVec(iSym))
    l_X = l_X**2
    call mma_allocate(X,l_X,Label='X')
    ! compute number of vector blocks
    nBlock = (nEnrVec(iSym)-1)/Laplace_BlockSize+1
    ! open vector file
    call ChoMP2_OpenF(1,iTyp,iSym)
    ! allocate memory for vectors
    l_Tot = Nai*nEnrVec(iSym)
    call mma_allocate(V,l_Tot,2,Label='V')
    ! Read all vectors
    iOpt = 2
    iAddr = 1
    call dDAFile(lUnit_F(iSym,iTyp),iOpt,V(:,1),l_Tot,iAddr)
    ! loop over Laplace grid
    do q=1,N
      ! init energy for this q
      Eq = 0.0d0
      ! scale grid point by 1/2
      tq = 0.5d0*t(q)
      ! scale vectors
      call dCopy_(l_Tot,V(:,1),1,V(:,2),1)
      do iVec=1,nEnrVec(iSym)
        ip0 = Nai*(iVec-1)
        do iSymi=1,nSym
          iSyma = MulD2h(iSym,iSymi)
          ip1 = ip0+iT1am(iSyma,iSymi)
          do i=1,nOcc(iSymi)
            call dScal_(nVir(iSyma),exp(epsi(i,iSymi)*tq),V(ip1+nVir(iSyma)*(i-1)+1,2),1)
          end do
          do a=1,nVir(iSyma)
            call dScal_(nOcc(iSymi),exp(-epsa(a,iSyma)*tq),V(ip1+a,2),nVir(iSyma))
          end do
        end do
      end do
      ! loop over vector blocks to compute
      ! X(J,K) = sum_ai L(ai,J)*L(ai,K)*exp(-(e(a)-a(i))*t(q)/2)
      ! Eq += w(q)*sum_JK [X(J,K)]**2
      do jBlock=1,nBlock
        ipj = 1+Nai*Laplace_BlockSize*(jBlock-1)
        if (jBlock == nBlock) then
          nVecj = nEnrVec(iSym)-Laplace_BlockSize*(nBlock-1)
        else
          nVecj = Laplace_BlockSize
        end if
        do iBlock=jBlock,nBlock
          ipi = 1+Nai*Laplace_BlockSize*(iBlock-1)
          if (iBlock == nBlock) then
            nVeci = nEnrVec(iSym)-Laplace_BlockSize*(nBlock-1)
          else
            nVeci = Laplace_BlockSize
          end if
          call dGEMM_('T','N',nVeci,nVecj,Nai,1.0d0,V(ipi,2),Nai,V(ipj,2),Nai,0.0d0,X,nVeci)
          if (iBlock == jBlock) then
            Eq = Eq+0.5d0*dDot_(nVeci*nVecj,X,1,X,1)
          else
            Eq = Eq+dDot_(nVeci*nVecj,X,1,X,1)
          end if
        end do
      end do
      ! accumulate in EMP2
      EMP2 = EMP2-w(q)*Eq
    end do
    ! deallocate memory
    call mma_deallocate(V)
    call mma_deallocate(X)
    ! close (and possibly delete) vector file
    call ChoMP2_OpenF(iClos,iTyp,iSym)
  end if
end do

! Scale energy (only half of X computed)
EMP2 = 2.0d0*EMP2

end subroutine ChoLSOSMP2_Energy_Fll1
