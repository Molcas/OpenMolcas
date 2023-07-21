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

subroutine ChoLSOSMP2_Energy_Srt(N,w,t,EOcc,EVir,Delete,EMP2,irc)
!
! Thomas Bondo Pedersen, December 2012.
!
! Compute Laplace-SOS-MP2 energy correction from sorted Cholesky
! vectors (i.e., occupied orbitals processed in batches).

use ChoMP2, only: iFirstS, LnOcc, LnT1am, LiT1am, lUnit
use Constants
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
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
character(len=21), parameter :: SecNam = 'ChoLSOSMP2_Energy_Srt'
real*8, external :: dDot_
#if !defined (_I8_) || defined (_DEBUGPRINT_)
character*2 Unt
real*8 Byte
#endif
integer nEnrVec(8)
integer iSym
integer l_X
integer Nai
integer nBlock
integer iOpt, iAddr, l_Tot
integer q
integer iBlock, jBlock
integer iSymi, iSyma
integer i, a, ii
integer ipi, ipj
integer nVeci, nVecj
integer iVec
integer ip0, ip1
integer ipX, lenX
integer iBatch
integer bsize
integer blast
real*8 lX, xM, xn, xb, xbp
real*8 tq, Eq, wq
real*8, allocatable :: X(:), V(:)
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
EMP2 = Zero

! check input (incl. common block variables)
if (nBatch < 2) then
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

! set number of vectors
if (DecoMP2) then
  call iCopy(nSym,nMP2Vec,1,nEnrVec,1)
else
  call iCopy(nSym,NumCho,1,nEnrVec,1)
end if

! allocate X
lX = Zero
do iSym=1,nSym
  if ((nT1am(iSym) > 0) .and. (nEnrVec(iSym) > 0)) then
    bsize = min(Laplace_BlockSize,nEnrVec(iSym))
    nBlock = (nEnrVec(iSym)-1)/bsize+1
    blast = nEnrVec(iSym)-bsize*(nBlock-1)
    xM = dble(nEnrVec(iSym))
    xn = dble(nBlock)
    xb = dble(bsize)
    xbp = dble(blast)
    lX = max(lX,0.5d0*(xM*(xM+1.0d0)+(xn-1.0d0)*xb*(xb-1.0d0)+xbp*(xbp-1.0d0)))
  end if
end do
l_X = int(lX)
#if !defined (_I8_) || defined (_DEBUGPRINT_)
if (l_X < 0) then
  write(Lupri,'(A,A)') SecNam,': dimension of X matrix is negative!'
  write(Lupri,'(A,I15)') 'l_X=',l_X
  if (lX > 0.0d0) then
    write(LuPri,'(A)') 'This seems to be an integer overflow!'
    call Cho_RWord2Byte(lX,Byte,Unt)
    write(LuPri,'(A,1P,D15.6,A,D15.6,1X,A,A)') 'In double precision, lX=',lX,' words (',Byte,Unt,')'
  end if
  irc = 1
  return
end if
#endif
call mma_allocate(X,l_X,Label='X')

! allocate vector array
Nai = 0
do iBatch=1,nBatch
  do iSym=1,nSym
    Nai = max(Nai,LnT1am(iSym,iBatch)*nEnrVec(iSym))
  end do
end do
call mma_allocate(V,Nai,Label='V')

! compute energy correction
do q=1,N
  ! init energy for this q
  Eq = 0.0d0
  ! scale grid point by 1/2
  tq = 0.5d0*t(q)
  ! scale weight by 2 (only lower blocks of X computed)
  wq = 2.0d0*w(q)
  do iSym=1,nSym
    if (nEnrVec(iSym) > 0) then
      ! init X for this symmetry
      bsize = min(Laplace_BlockSize,nEnrVec(iSym))
      nBlock = (nEnrVec(iSym)-1)/bsize+1
      blast = nEnrVec(iSym)-bsize*(nBlock-1)
      lenX = nEnrVec(iSym)*(nEnrVec(iSym)+1)/2+(nBlock-1)*(bsize*(bsize-1)/2)+blast*(blast-1)/2
#     ifdef _DEBUGPRINT_
      if (lenX > l_X) then
        call WarningMessage(2,SecNam//': insufficient X allocation')
        write(6,'(A,2(1X,I10))') 'lenX,l_X=',lenX,l_X
        call Abend()
      end if
#     endif
      X(1:lenX) = Zero
      do iBatch=1,nBatch
        Nai = LnT1am(iSym,iBatch)
        if (Nai > 0) then
          ! open file, read vectors, close file
          ! do not delete file here - it may be needed later
          ! (because of the loop over q)
          call ChoMP2_OpenB(1,iSym,iBatch)
          iOpt = 2
          l_Tot = Nai*nEnrVec(iSym)
          iAddr = 1
#         ifdef _DEBUGPRINT_
          if (l_Tot > size(V)) then
            call WarningMessage(2,SecNam//': insufficient V allocation')
            write(6,'(A,2(1X,I10))') 'l_Tot,SIZE(V)=',l_Tot,size(V)
            call Abend()
          end if
#         endif
          call dDAFile(lUnit(iSym,iBatch),iOpt,V,l_Tot,iAddr)
          call ChoMP2_OpenB(2,iSym,iBatch)
          ! scale vectors
          do iVec=1,nEnrVec(iSym)
            ip0 = Nai*(iVec-1)
            do iSymi=1,nSym
              if (LnOcc(iSymi,iBatch) > 0) then
                iSyma = MulD2h(iSym,iSymi)
                ip1 = ip0+LiT1am(iSyma,iSymi,iBatch)
                do i=0,LnOcc(iSymi,iBatch)-1
                  ii = iFirstS(iSymi,iBatch)+i
                  call dScal_(nVir(iSyma),exp(epsi(ii,iSymi)*tq),V(ip1+nVir(iSyma)*i+1),1)
                end do
                do a=1,nVir(iSyma)
                  call dScal_(LnOcc(iSymi,iBatch),exp(-epsa(a,iSyma)*tq),V(ip1+a),nVir(iSyma))
                end do
              end if
            end do
          end do
          ! loop over vector blocks to compute
          ! X(J,K) +=
          ! sum_ai L(ai,J)*L(ai,K)*exp(-(e(a)-e(i))*t(q)/2)
          ipX = 1
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
              call dGEMM_('T','N',nVeci,nVecj,Nai,1.0d0,V(ipi),Nai,V(ipj),Nai,1.0d0,X(ipX),nVeci)
              ipX = ipX+nVeci*nVecj
            end do
          end do
#         ifdef _DEBUGPRINT_
          if (lenX /= (ipX-1)) then
            call WarningMessage(2,SecNam//': dimension problem [1]')
            write(6,'(A,I10,A,I10)') 'lenX=',lenX,' ipX-1=',ipX-1
            call Abend()
          end if
#         endif
        end if
      end do
      ! compute energy contribution
      ! Eq += sum_JK [X(J,K)]**2
      ipX = 1
      do jBlock=1,nBlock
        if (jBlock == nBlock) then
          nVecj = nEnrVec(iSym)-Laplace_BlockSize*(nBlock-1)
        else
          nVecj = Laplace_BlockSize
        end if
        do iBlock=jBlock,nBlock
          if (iBlock == nBlock) then
            nVeci = nEnrVec(iSym)-Laplace_BlockSize*(nBlock-1)
          else
            nVeci = Laplace_BlockSize
          end if
          if (iBlock == jBlock) then
            Eq = Eq+0.5d0*dDot_(nVeci*nVecj,X(ipX),1,X(ipX),1)
          else
            Eq = Eq+dDot_(nVeci*nVecj,X(ipX),1,X(ipX),1)
          end if
          ipX = ipX+nVeci*nVecj
        end do
      end do
#     ifdef _DEBUGPRINT_
      if (lenX /= (ipX-1)) then
        call WarningMessage(2,SecNam//': dimension problem [2]')
        write(6,'(A,I10,A,I10)') 'lenX=',lenX,' ipX-1=',ipX-1
        call Abend()
      end if
#     endif
    end if
  end do
  ! scale Eq
  Eq = -wq*Eq
  ! Accumulate in EMP2
  EMP2 = EMP2+Eq
end do

! deallocations
call mma_deallocate(V)
call mma_deallocate(X)

! delete files if requested
if (Delete) then
  do iBatch=1,nBatch
    do iSym=1,nSym
      call ChoMP2_OpenB(1,iSym,iBatch)
      call ChoMP2_OpenB(3,iSym,iBatch)
    end do
  end do
end if

end subroutine ChoLSOSMP2_Energy_Srt
